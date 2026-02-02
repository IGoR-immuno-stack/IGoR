#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
source $SCRIPT_DIR/config.sh

BENCH_DIR=$(mktemp -d)
LOG_FILE="$BENCH_DIR/benchmark.log"

# Default sizes
# Standalone sampling sizes (Generation only)
GEN_SIZES=(1000 10000 100000 1000000)
# Pipeline sizes (Generate -> Read -> Align -> Infer)
PIPELINE_SIZES=(100 500 1000)

IGORCALL="$IGORBIN -set_wd $BENCH_DIR"
MODEL_PARMS="$TESTREF/demo_inference/final_parms.txt"
MODEL_MARGINALS="$TESTREF/demo_inference/final_marginals.txt"
GENOMIC_V="$TESTINPUT/genomicVs_with_primers.fasta"
GENOMIC_D="$TESTINPUT/genomicDs.fasta"
GENOMIC_J="$TESTINPUT/genomicJs_all_curated.fasta"

# Results storage
OPTIONAL_ARGS=()
declare -a RESULTS
CURRENT_CONTEXT=""

VERBOSE=false

# Helper to process args
process_args() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --verbose)
                VERBOSE=true
                shift
                ;;
            --)
                shift
                ;;
            *)
                # Store positional args
                if [[ -z "${MODE:-}" ]]; then
                    MODE="$1"
                else
                    OPTIONAL_ARGS+=("$1")
                fi
                shift
                ;;
        esac
    done
}

process_args "$@"
MODE="${MODE:-all}"

# Helper for timing
get_time() {
    python3 -c 'import time; print(time.time())'
}

# Helper to run silent and capture time
run_step() {
    local step_name="$1"
    local command_str="$2"
    
    # Indented output for steps
    printf "  %-25s : " "$step_name"
    
    local start=$(get_time)
    # Run command and redirect ALL output to log
    echo "CMD: $command_str" >> "$LOG_FILE"
    if eval "$command_str" >> "$LOG_FILE" 2>&1; then
        local end=$(get_time)
        local duration=$(python3 -c "print(f'{($end - $start):.4f}')")
        echo "${duration}s"
        RESULTS+=("$CURRENT_CONTEXT|$step_name|$duration")
    else
        echo "FAILED (See log)"
        RESULTS+=("$CURRENT_CONTEXT|$step_name|FAILED")
        return 1
    fi
}

run_standalone_generation() {
    local n_seqs="$1"
    local batch_name="gen_standalone_${n_seqs}"
    local cmd="$IGORCALL -batch $batch_name -threads 1 -set_custom_model \"$MODEL_PARMS\" \"$MODEL_MARGINALS\" -generate $n_seqs --seed 12345"
    
    echo ""
    echo "------------------------------------------------------------------"
    echo "[Standalone] Generation (N=$n_seqs)"
    CURRENT_CONTEXT="Standalone (N=$n_seqs)"
    
    run_step "Generation" "$cmd"
}

run_pipeline_bench() {
    local n_seqs="$1"
    local n_threads="$2"
    local pfx="pipe_${n_seqs}_${n_threads}"
    
    echo ""
    echo "------------------------------------------------------------------"
    echo "[Pipeline] Full Workflow (N=$n_seqs, Threads=$n_threads)"
    CURRENT_CONTEXT="Pipeline (N=$n_seqs, T=$n_threads)"

    # 1. Pipeline Generation
    local gen_cmd="$IGORCALL -batch $pfx -threads 1 -set_custom_model \"$MODEL_PARMS\" \"$MODEL_MARGINALS\" -generate $n_seqs --seed 12345"
    run_step "Generation" "$gen_cmd"
    
    local seq_file="$BENCH_DIR/${pfx}_generated/generated_seqs_werr.csv"

    # 2. Pipeline Read
    local read_cmd="$IGORCALL -batch $pfx -read_seqs \"$seq_file\""
    run_step "Read Seqs" "$read_cmd"

    # 3. Pipeline Align
    if [ "$VERBOSE" = true ]; then
        local align_v="$IGORCALL -batch $pfx -threads $n_threads -align --V -set_genomic --V \"$GENOMIC_V\""
        run_step "Align (V)" "$align_v"

        local align_d="$IGORCALL -batch $pfx -threads $n_threads -align --D -set_genomic --D \"$GENOMIC_D\""
        run_step "Align (D)" "$align_d"

        local align_j="$IGORCALL -batch $pfx -threads $n_threads -align --J -set_genomic --J \"$GENOMIC_J\""
        run_step "Align (J)" "$align_j"
    else
        local align_cmd="$IGORCALL -batch $pfx -threads $n_threads -align --V -set_genomic --V \"$GENOMIC_V\" && \
                         $IGORCALL -batch $pfx -threads $n_threads -align --D -set_genomic --D \"$GENOMIC_D\" && \
                         $IGORCALL -batch $pfx -threads $n_threads -align --J -set_genomic --J \"$GENOMIC_J\""
        run_step "Alignment" "$align_cmd"
    fi

    # 4. Pipeline Infer
    local infer_cmd="$IGORCALL -batch $pfx -threads $n_threads -set_custom_model \"$MODEL_PARMS\" \"$MODEL_MARGINALS\" -infer --N_iter 2"
    run_step "Inference (2 iter)" "$infer_cmd"
}

print_summary() {
    echo ""
    echo "=================================================================="
    echo "BENCHMARK SUMMARY"
    echo "=================================================================="
    printf "%-30s | %-25s | %s\n" "Context" "Step" "Time (s)"
    echo "------------------------------------------------------------------"
    for row in "${RESULTS[@]}"; do
        IFS='|' read -r context step time <<< "$row"
        printf "%-30s | %-25s | %s\n" "$context" "$step" "$time"
    done
    echo "------------------------------------------------------------------"
    echo "Logs: $LOG_FILE"
}

# ---------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------

echo "Starting Benchmarks... (Logs: $LOG_FILE)"

if [[ "$MODE" == "sampling" ]]; then
    if [[ ${#OPTIONAL_ARGS[@]} -gt 0 ]]; then
        for size in "${OPTIONAL_ARGS[@]}"; do
            run_standalone_generation "$size"
        done
    else
        for size in "${GEN_SIZES[@]}"; do
            run_standalone_generation $size
        done
    fi
elif [[ "$MODE" == "gen" ]]; then
    for size in "${GEN_SIZES[@]}"; do
        run_standalone_generation $size
    done
fi

if [[ "$MODE" == "pipeline" || "$MODE" == "all" ]]; then
    for size in "${PIPELINE_SIZES[@]}"; do
        run_pipeline_bench $size 1
        
        if [ "$size" -ge 500 ]; then
            run_pipeline_bench $size 4
        fi
    done
fi

if [[ "$MODE" == "all" ]]; then
   for size in "${GEN_SIZES[@]}"; do
       run_standalone_generation $size
   done
fi

print_summary

# Cleanup
# rm -rf "$BENCH_DIR"
