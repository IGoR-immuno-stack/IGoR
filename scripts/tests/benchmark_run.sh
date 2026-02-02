#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
source $SCRIPT_DIR/config.sh

BENCH_DIR=$(mktemp -d)
LOG_FILE="$BENCH_DIR/benchmark.log"

# Default sizes
# Standalone sampling sizes (Generation only)
GEN_SIZES=(10000 100000)
# Pipeline sizes (Generate -> Read -> Align -> Infer)
# Alignment is computationally expensive (~0.3s/seq), so defaults are kept small.
PIPELINE_SIZES=(100 500 1000)

IGORCALL="$IGORBIN -set_wd $BENCH_DIR"
MODEL_PARMS="$TESTREF/demo_inference/final_parms.txt"
MODEL_MARGINALS="$TESTREF/demo_inference/final_marginals.txt"
GENOMIC_V="$TESTINPUT/genomicVs_with_primers.fasta"
GENOMIC_D="$TESTINPUT/genomicDs.fasta"
GENOMIC_J="$TESTINPUT/genomicJs_all_curated.fasta"

# Results storage
declare -a RESULTS

# Helper for timing
get_time() {
    python3 -c 'import time; print(time.time())'
}

# Helper to run silent and capture time
run_step() {
    local step_name="$1"
    local unique_id="$2"
    local command_str="$3"
    
    printf "%-20s | %-30s | " "$step_name" "$unique_id"
    
    local start=$(get_time)
    # Run command and redirect ALL output to log
    echo "CMD: $command_str" >> "$LOG_FILE"
    if eval "$command_str" >> "$LOG_FILE" 2>&1; then
        local end=$(get_time)
        local duration=$(python3 -c "print(f'{($end - $start):.4f}')")
        echo "${duration}s"
        RESULTS+=("$step_name|$unique_id|$duration")
    else
        echo "FAILED (See log)"
        RESULTS+=("$step_name|$unique_id|FAILED")
        return 1
    fi
}

run_standalone_generation() {
    local n_seqs="$1"
    local batch_name="gen_standalone_${n_seqs}"
    # Use -threads 1 for consistency in simple throughput tests, 
    # or we could expose threads if IGoR supported parallel generation (usually it's per-sequence fast anyway).
    local cmd="$IGORCALL -batch $batch_name -threads 1 -set_custom_model \"$MODEL_PARMS\" \"$MODEL_MARGINALS\" -generate $n_seqs --seed 12345"
    
    run_step "Standalone_Sampling" "N=$n_seqs" "$cmd"
}

run_pipeline_bench() {
    local n_seqs="$1"
    local n_threads="$2"
    local pfx="pipe_${n_seqs}_${n_threads}"
    
    # 1. Pipeline Generation
    # This generation is specifically for the pipeline workflow
    local gen_cmd="$IGORCALL -batch $pfx -threads 1 -set_custom_model \"$MODEL_PARMS\" \"$MODEL_MARGINALS\" -generate $n_seqs --seed 12345"
    run_step "Pipeline_Gen" "N=$n_seqs, T=1" "$gen_cmd"
    
    local seq_file="$BENCH_DIR/${pfx}_generated/generated_seqs_werr.csv"

    # 2. Pipeline Read
    local read_cmd="$IGORCALL -batch $pfx -read_seqs \"$seq_file\""
    run_step "Pipeline_Read" "N=$n_seqs" "$read_cmd"

    # 3. Pipeline Align
    local align_cmd="$IGORCALL -batch $pfx -threads $n_threads -align --V -set_genomic --V \"$GENOMIC_V\" && \
                     $IGORCALL -batch $pfx -threads $n_threads -align --D -set_genomic --D \"$GENOMIC_D\" && \
                     $IGORCALL -batch $pfx -threads $n_threads -align --J -set_genomic --J \"$GENOMIC_J\""
    run_step "Pipeline_Align" "N=$n_seqs, T=$n_threads" "$align_cmd"

    # 4. Pipeline Infer
    local infer_cmd="$IGORCALL -batch $pfx -threads $n_threads -set_custom_model \"$MODEL_PARMS\" \"$MODEL_MARGINALS\" -infer --N_iter 2"
    run_step "Pipeline_Infer" "N=$n_seqs, T=$n_threads" "$infer_cmd"
}

print_summary() {
    echo ""
    echo "=================================================================="
    echo "BENCHMARK SUMMARY"
    echo "=================================================================="
    printf "%-20s | %-25s | %s\n" "Operation" "Configuration" "Time (s)"
    echo "------------------------------------------------------------------"
    for row in "${RESULTS[@]}"; do
        IFS='|' read -r op conf time <<< "$row"
        printf "%-20s | %-25s | %s\n" "$op" "$conf" "$time"
    done
    echo "------------------------------------------------------------------"
    echo "Logs: $LOG_FILE"
}

# ---------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------

# Handle cases where arguments are passed with a leading '--'
if [[ "${1:-}" == "--" ]]; then
    shift
fi

MODE="${1:-all}"
OPTIONAL_ARG="${2:-}"

echo "Starting Benchmarks... (Logs: $LOG_FILE)"
echo "------------------------------------------------------------------"
printf "%-20s | %-30s | %s\n" "Step" "Config" "Time"
echo "------------------------------------------------------------------"

if [[ "$MODE" == "sampling" ]]; then
    if [[ -n "$OPTIONAL_ARG" ]]; then
        # User specified size (e.g., pixi run benchmark -- sampling 1000000)
        run_standalone_generation "$OPTIONAL_ARG"
    else
        # Run default sampling sizes
        for size in "${GEN_SIZES[@]}"; do
            run_standalone_generation $size
        done
    fi
elif [[ "$MODE" == "gen" ]]; then
     # Alias for default sampling
    for size in "${GEN_SIZES[@]}"; do
        run_standalone_generation $size
    done
fi

if [[ "$MODE" == "pipeline" || "$MODE" == "all" ]]; then
    for size in "${PIPELINE_SIZES[@]}"; do
        # Single Thread
        run_pipeline_bench $size 1
        
        # Multi Thread (only for larger sizes)
        if [ "$size" -ge 500 ]; then
            run_pipeline_bench $size 4
        fi
    done
fi

if [[ "$MODE" == "all" ]]; then
   # Also run the standard sampling sizes for 'all'
   for size in "${GEN_SIZES[@]}"; do
       run_standalone_generation $size
   done
fi

print_summary

# Cleanup?
# rm -rf "$BENCH_DIR"
