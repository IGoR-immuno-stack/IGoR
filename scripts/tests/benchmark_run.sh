#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
source $SCRIPT_DIR/config.sh

# ------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------
BENCH_DIR=$(mktemp -d)
LOG_FILE="$BENCH_DIR/benchmark.log"
REPORT_DIR="benchmark_reports"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
REPORT_BASE_DIR="$REPORT_DIR/benchmark_$TIMESTAMP"
mkdir -p "$REPORT_BASE_DIR"

# Default sizes
GEN_SIZES=(1000 10000 100000 1000000)
PIPELINE_SIZES=(100 500 1000)

IGORCALL="$IGORBIN -set_wd $BENCH_DIR"
MODEL_PARMS="$TESTREF/demo_inference/final_parms.txt"
MODEL_MARGINALS="$TESTREF/demo_inference/final_marginals.txt"
GENOMIC_V="$TESTINPUT/genomicVs_with_primers.fasta"
GENOMIC_D="$TESTINPUT/genomicDs.fasta"
GENOMIC_J="$TESTINPUT/genomicJs_all_curated.fasta"

# Results storage
declare -a RESULTS
VERBOSE=false


SPECIFIC_BENCHMARK=""

# ------------------------------------------------------------------
# Benchmark definitions
# ------------------------------------------------------------------
# Each benchmark: index|name|size|threads|context
declare -a BENCHMARK_DEFINITIONS=(
    "1|gen_1k|1000|1|Standalone (N=1k)"
    "2|gen_10k|10000|1|Standalone (N=10k)"
    "3|gen_100k|100000|1|Standalone (N=100k)"
    "4|gen_1M|1000000|1|Standalone (N=1M)"
    "5|pipe_100_t1|100|1|Pipeline (N=100, T=1)"
    "6|pipe_500_t1|500|1|Pipeline (N=500, T=1)"
    "7|pipe_1000_t1|1000|1|Pipeline (N=1000, T=1)"
    "8|pipe_500_t4|500|4|Pipeline (N=500, T=4)"
    "9|pipe_1000_t4|1000|4|Pipeline (N=1000, T=4)"
)

# ------------------------------------------------------------------
# Helper functions (defined before they are used)
# ------------------------------------------------------------------
list_benchmarks() {
    echo "Available benchmarks:"
    echo "---------------------"
    echo ""
    echo "STANDALONE SAMPLING (N=1k to 1M)"
    printf "%-3s | %-15s | %-12s | %s\n" "Idx" "Name" "Size" "Context"
    echo "-----|-----------------|-------------|------------------------------------"
    for def in "${BENCHMARK_DEFINITIONS[@]:0:4}"; do
        IFS='|' read -r idx name size threads context <<< "$def"
        printf "%-3s | %-15s | %-12s | %s\n" "$idx" "$name" "${size} seqs" "$context"
    done
    echo ""
    echo "FULL PIPELINE (Generate -> Read -> Align -> Infer)"
    printf "%-3s | %-15s | %-12s | %s\n" "Idx" "Name" "Size" "Context"
    echo "-----|-----------------|-------------|------------------------------------"
    for def in "${BENCHMARK_DEFINITIONS[@]:4}"; do
        IFS='|' read -r idx name size threads context <<< "$def"
        printf "%-3s | %-15s | %-12s | %s\n" "$idx" "$name" "${size} seqs" "$context"
    done
    echo ""
}

get_time() {
    python3 -c 'import time; print(time.time())'
}

get_benchmark_def() {
    local idx="$1"
    for def in "${BENCHMARK_DEFINITIONS[@]}"; do
        IFS='|' read -r def_idx name size threads context <<< "$def"
        if [[ "$def_idx" == "$idx" ]]; then
            echo "$def"
            return
        fi
    done
    echo ""
}

parse_benchmark_spec() {
    local spec="$1"
    local benchmarks=()

    if [[ "$spec" == "all" ]] || [[ -z "$spec" ]]; then
        for def in "${BENCHMARK_DEFINITIONS[@]}"; do
            IFS='|' read -r idx name size threads context <<< "$def"
            benchmarks+=("$idx")
        done
        echo "${benchmarks[@]}"
        return
    fi

    # Handle comma-separated values
    IFS=',' read -ra ITEMS <<< "$spec"
    for item in "${ITEMS[@]}"; do
        item=$(echo "$item" | xargs)  # trim whitespace

        # Handle range (e.g., 1-5)
        if [[ "$item" =~ ^[0-9]+-[0-9]+$ ]]; then
            local start=$(echo "$item" | cut -d'-' -f1)
            local end=$(echo "$item" | cut -d'-' -f2)
            for ((i=start; i<=end; i++)); do
                if [[ -n "$(get_benchmark_def "$i")" ]]; then
                    benchmarks+=("$i")
                else
                    echo "Error: Invalid benchmark index '$i'" >&2
                    exit 1
                fi
            done
        # Single index
        elif [[ "$item" =~ ^[0-9]+$ ]]; then
            if [[ -n "$(get_benchmark_def "$item")" ]]; then
                benchmarks+=("$item")
            else
                echo "Error: Invalid benchmark index '$item'" >&2
                exit 1
            fi
        else
            echo "Error: Invalid benchmark specification '$item'" >&2
            echo "Use indices (e.g., 1, 5-8, or see --list)" >&2
            exit 1
        fi
    done

    # Remove duplicates and sort
    printf '%s\n' "${benchmarks[@]}" | sort -nu | tr '\n' ' '
}

# ------------------------------------------------------------------
# Usage
# ------------------------------------------------------------------
usage() {
    cat << EOF
Usage: $0 [OPTIONS] [MODE] [EXTRA_ARGS...]

Run IGoR performance benchmarks and generate detailed reports.

Arguments:
  MODE                Benchmark mode: all, sampling, pipeline, gen
                      Default: all

  EXTRA_ARGS          Additional arguments depending on mode:
                      - For 'sampling'/'gen': custom sizes [size1 size2 ...]
                      - For 'pipeline': (none, uses default sizes)

  BENCHMARK_ID        Specific benchmark index to run (overrides MODE)
                      See --list for available indices

Options:
  -h, --help          Show this help message
  -l, --list          List available benchmarks with indices
  -v, --verbose       Show detailed output for each step
  -i, --id ID         Run specific benchmark by index
  -k, --keep          Keep temporary output files

Benchmark Categories:
  Standalone Sampling (N=1k to 1M):
    1-4: Generation with varying sequence counts

  Full Pipeline (N=100,500,1000; T=1,4):
    5-9: Generate -> Read -> Align -> Infer workflow

Examples (pixi):
  pixi run benchmark                  Run all benchmarks
  pixi run benchmark all              Run all benchmarks (same as above)
  pixi run benchmark sampling         Run standalone sampling benchmarks
  pixi run benchmark pipeline         Run full pipeline benchmarks
  pixi run benchmark sampling 5000 50000  Custom sampling sizes
  pixi run benchmark 1                Run benchmark #1 only
  pixi run benchmark 5-8              Run benchmarks #5 through #8
  pixi run benchmark 1,3,5            Run benchmarks 1, 3, and 5
  pixi run benchmark --list           List available benchmarks
  pixi run benchmark -v --list        Verbose listing

Direct execution:
  $0                           Run all benchmarks
  $0 1                          Run benchmark #1 only
  $0 -l                         List available benchmarks

EOF
}

# ------------------------------------------------------------------
# Parse arguments
# ------------------------------------------------------------------
KEEP_OUTPUT=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            usage
            exit 0
            ;;
        -l|--list)
            list_benchmarks
            exit 0
            ;;
        -v|--verbose)
            VERBOSE=true
            shift
            ;;
        -k|--keep)
            KEEP_OUTPUT=true
            shift
            ;;
        -i|--id)
            if [[ -z "${2:-}" ]]; then
                echo "Error: --id requires an argument" >&2
                exit 1
            fi
            SPECIFIC_BENCHMARK="$2"
            shift 2
            ;;
        --)
            shift
            break
            ;;
        -*)
            echo "Error: Unknown option '$1'" >&2
            usage >&2
            exit 1
            ;;
        *)
            break
            ;;
    esac
done

# ------------------------------------------------------------------
# Benchmark execution functions
# ------------------------------------------------------------------
run_standalone_generation() {
    local idx="$1"
    local n_seqs="$2"
    local batch_name="gen_${n_seqs}"
    local cmd="$IGORCALL -batch $batch_name -threads 1 -set_custom_model \"$MODEL_PARMS\" \"$MODEL_MARGINALS\" -generate $n_seqs --seed 12345"

    echo ""
    echo "------------------------------------------------------------------"
    echo "[${idx}] Standalone Generation (N=$n_seqs)"
    CURRENT_CONTEXT="Standalone (N=$n_seqs)"

    run_step "Generation" "$cmd"
}

run_pipeline_bench() {
    local idx="$1"
    local n_seqs="$2"
    local n_threads="$3"
    local pfx="pipe_${n_seqs}_${n_threads}"

    echo ""
    echo "------------------------------------------------------------------"
    echo "[${idx}] Full Pipeline (N=$n_seqs, Threads=$n_threads)"
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

run_step() {
    local step_name="$1"
    local command_str="$2"

    printf "  %-25s : " "$step_name"

    # Show command if verbose mode
    if $VERBOSE; then
        printf '\n    $ %s\n' "$command_str"
    fi

    local start=$(get_time)
    echo "CMD: $command_str" >> "$LOG_FILE"
    if eval "$command_str" >> "$LOG_FILE" 2>&1; then
        local end=$(get_time)
        local duration=$(python3 -c "print(f'{($end - $start):.4f}')")
        echo "${duration}s"
        RESULTS+=("$CURRENT_CONTEXT|$step_name|$duration|OK")
    else
        echo "FAILED (See log)"
        RESULTS+=("$CURRENT_CONTEXT|$step_name|FAILED|ERROR")
        return 1
    fi
}

print_summary() {
    echo ""
    echo "=================================================================="
    echo "BENCHMARK SUMMARY"
    echo "=================================================================="
    printf "%-30s | %-25s | %-10s | %s\n" "Context" "Step" "Time (s)" "Status"
    echo "------------------------------------------------------------------"
    for row in "${RESULTS[@]}"; do
        IFS='|' read -r context step time status <<< "$row"
        printf "%-30s | %-25s | %-10s | %s\n" "$context" "$step" "$time" "$status"
    done
    echo "------------------------------------------------------------------"
    echo "Full log: $LOG_FILE"
    echo "Report:   $REPORT_BASE_DIR/SUMMARY.txt"
}

print_detailed_summary() {
    local summary_file="$REPORT_BASE_DIR/SUMMARY.txt"

    cat > "$summary_file" << EOF
==============================================================================
BENCHMARK SUMMARY REPORT
==============================================================================
Timestamp:  $(date)
Log file:   $LOG_FILE

------------------------------------------------------------------------------
DETAILED RESULTS
------------------------------------------------------------------------------
EOF

    printf "%-30s | %-25s | %-10s | %s\n" "Context" "Step" "Time (s)" "Status" >> "$summary_file"
    echo "------------------------------------------------------------------" >> "$summary_file"

    for row in "${RESULTS[@]}"; do
        IFS='|' read -r context step time status <<< "$row"
        printf "%-30s | %-25s | %-10s | %s\n" "$context" "$step" "$time" "$status" >> "$summary_file"
    done
    echo "------------------------------------------------------------------" >> "$summary_file"
}

# ------------------------------------------------------------------
# MAIN
# ------------------------------------------------------------------
MODE="${1:-}"
# Track if we're using custom sizes
USE_CUSTOM_SIZES=false
CUSTOM_SIZES=()

if [[ -n "$SPECIFIC_BENCHMARK" ]]; then
    # Specific benchmark ID overrides mode
    BENCHMARKS_TO_RUN=($(parse_benchmark_spec "$SPECIFIC_BENCHMARK"))
else
    # Check if first arg is a number (benchmark ID)
    if [[ "${1:-}" =~ ^[0-9]+$ ]] || [[ "${1:-}" =~ ^[0-9]+-[0-9]+$ ]] || [[ "${1:-}" =~ ^[0-9,]+$ ]]; then
        BENCHMARKS_TO_RUN=($(parse_benchmark_spec "$1"))
        shift
        MODE=""
    else
        MODE="${1:-all}"
        shift

        case "$MODE" in
            sampling|gen)
                # Use custom sizes if provided
                if [[ $# -gt 0 ]]; then
                    CUSTOM_SIZES=("$@")
                    USE_CUSTOM_SIZES=true
                fi
                # Build dynamic benchmark list
                BENCHMARKS_TO_RUN=()
                if $USE_CUSTOM_SIZES; then
                    for ((i=0; i<${#CUSTOM_SIZES[@]}; i++)); do
                        BENCHMARKS_TO_RUN+=("$((i+1))")
                    done
                else
                    for ((i=0; i<${#GEN_SIZES[@]}; i++)); do
                        BENCHMARKS_TO_RUN+=("$((i+1))")
                    done
                fi
                ;;
            pipeline)
                BENCHMARKS_TO_RUN=(5 6 7 8 9)
                ;;
            *)
                BENCHMARKS_TO_RUN=($(parse_benchmark_spec "all"))
                ;;
        esac
    fi
fi

echo "Starting Benchmarks..."
echo "Log file: $LOG_FILE"
echo "Report:   $REPORT_BASE_DIR"
if $USE_CUSTOM_SIZES; then
    echo "Custom sizes: ${CUSTOM_SIZES[*]}"
fi
if $VERBOSE; then
    echo "Verbose: enabled"
fi
echo "Benchmarks: ${BENCHMARKS_TO_RUN[*]}"
echo ""

for bench_idx in "${BENCHMARKS_TO_RUN[@]}"; do
    # For custom sizes in sampling mode, use the custom size instead of definition
    actual_size=""
    array_idx=""
    if $USE_CUSTOM_SIZES && [[ "$MODE" == "sampling" || "$MODE" == "gen" ]]; then
        # Get size from CUSTOM_SIZES array (index is 1-based)
        array_idx=$((bench_idx - 1))
        if [[ $array_idx -lt ${#CUSTOM_SIZES[@]} ]]; then
            actual_size="${CUSTOM_SIZES[$array_idx]}"
        else
            # Fallback to definition if out of bounds
            def=$(get_benchmark_def "$bench_idx")
            IFS='|' read -r idx name actual_size threads context <<< "$def"
        fi
        # Reconstruct definition with custom size
        def="$bench_idx|gen_${actual_size}|${actual_size}|1|Standalone (N=${actual_size})"
    else
        def=$(get_benchmark_def "$bench_idx")
    fi
    
    if [[ -z "$def" ]]; then
        echo "Warning: Skipping invalid benchmark index $bench_idx"
        continue
    fi

    IFS='|' read -r idx name size threads context <<< "$def"

    # Check if it's standalone or pipeline
    if [[ "$name" == gen_* ]]; then
        run_standalone_generation "$idx" "$size"
    elif [[ "$name" == pipe_* ]]; then
        run_pipeline_bench "$idx" "$size" "$threads"
    fi
done

print_summary
print_detailed_summary

# Cleanup
if ! $KEEP_OUTPUT; then
    rm -rf "$BENCH_DIR"
    echo ""
    echo "Cleaned temporary files: $BENCH_DIR"
fi

exit 0
