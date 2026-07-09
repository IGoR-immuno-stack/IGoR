#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
source $SCRIPT_DIR/config.sh

# ------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------
PROF_DIR=$(mktemp -d)
trap 'rm -rf "$PROF_DIR"' EXIT

REPORT_DIR="profile_reports"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

# Fixed target: smallest pipeline size, single thread by default (matches
# benchmark_run.sh's pipe_100_t1 configuration). Thread count can be
# overridden with --threads, e.g. to reproduce false-sharing behaviour.
N_SEQS=100
N_THREADS=1
PFX="profile_run"

# perf record sampling frequency for the general (call-graph) session, and
# for perf c2c (cache-to-cache) sessions -- c2c needs a much higher rate
# since HITM/false-sharing events are comparatively rare.
PERF_FREQ=999
C2C_FREQ=60000

# Exact aggregate hardware-counter totals (perf stat, counting mode) to go
# alongside the sampled perf record data -- gives IPC, cache/branch miss
# rates and frontend/backend stall cycles with no sampling error, at the
# cost of a separate run of the same command.
STAT_EVENTS="cycles,instructions,L1-dcache-loads,L1-dcache-load-misses,LLC-loads,LLC-load-misses,branches,branch-misses,stalled-cycles-backend,stalled-cycles-frontend"

MODEL_PARMS="$TESTREF/demo_inference/final_parms.txt"
MODEL_MARGINALS="$TESTREF/demo_inference/final_marginals.txt"
GENOMIC_V="$TESTINPUT/genomicVs_with_primers.fasta"
GENOMIC_D="$TESTINPUT/genomicDs.fasta"
GENOMIC_J="$TESTINPUT/genomicJs_all_curated.fasta"

IGORCALL="$IGORBIN -w $PROF_DIR"

# ------------------------------------------------------------------
# Usage
# ------------------------------------------------------------------
usage() {
    cat << EOF
Usage: $0 [BACKEND] [MODE] [OPTIONS]

Profile IGoR's generate/align/infer subcommands with perf or
valgrind/callgrind, on a fixed smallest-sample-size (N=$N_SEQS) pipeline
run. Cache-miss statistics are collected by default with both backends.
With the perf backend, each mode also gets an exact "perf stat" pass
(cycles, instructions, L1-dcache and LLC loads/misses, branches,
branch-misses, stalled-cycles-frontend/backend) alongside the sampled
perf record data.

BACKEND   perf | callgrind | auto   (default: auto)
          auto tries perf first, and falls back to callgrind if perf is
          unavailable or blocked by /proc/sys/kernel/perf_event_paranoid.

MODE      generate | align | infer | all   (default: all)
          Each mode runs in its own independent profiling session (a
          separate perf record / valgrind invocation and output file) --
          "all" profiles the three modes one after another and never
          merges sessions.

Options:
  -h, --help              Show this help message
  --threads=N             Threads passed to igor's -j flag (default: 1).
                           Override to reproduce OpenMP scaling / false
                           sharing behaviour, e.g. with --c2c.
  --collect-fn=FUNCTION   callgrind only. Disables collection at process
                           start (--collect-atstart=no) and toggles it on
                           entering/leaving FUNCTION (--toggle-collect),
                           so CLI parsing, config/model loading, and other
                           setup/IO outside FUNCTION isn't counted. Repeat
                           to toggle on multiple functions. Ignored (with
                           a warning) when the backend is perf.
  --c2c                    perf only. Additionally records a "perf c2c"
                           (cache-to-cache) session per mode, to surface
                           false sharing between OpenMP threads. Ignored
                           (with a warning) when the backend is callgrind.
                           Meaningless with the default --threads=1.

Examples (pixi):
  pixi run profile                       Auto-detect backend, profile all modes
  pixi run profile -- perf generate      Force perf, profile generation only
  pixi run profile -- callgrind infer    Force callgrind, profile inference only
  pixi run profile -- callgrind align --collect-fn=Aligner::align_seqs
                                          Force callgrind, only collect inside
                                          align_seqs (skip setup/IO)
  pixi run profile -- perf align --threads=4 --c2c
                                          Profile 4-threaded alignment and
                                          look for false sharing

Output:
  Reports are written to profile_reports/profile_<timestamp>_<backend>/<mode>/,
  with the raw profiler data at its tool-default filename (perf.data or
  callgrind.out) so it can be opened without needing to pass -i/--in.
EOF
}

# ------------------------------------------------------------------
# Parse arguments
# ------------------------------------------------------------------
BACKEND="auto"
MODE="all"
COLLECT_FNS=()
C2C=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            usage
            exit 0
            ;;
        --threads=*)
            N_THREADS="${1#--threads=}"
            shift
            ;;
        --threads)
            if [[ -z "${2:-}" ]]; then
                echo "Error: --threads requires an argument" >&2
                exit 1
            fi
            N_THREADS="$2"
            shift 2
            ;;
        --collect-fn=*)
            COLLECT_FNS+=("${1#--collect-fn=}")
            shift
            ;;
        --collect-fn)
            if [[ -z "${2:-}" ]]; then
                echo "Error: --collect-fn requires an argument" >&2
                exit 1
            fi
            COLLECT_FNS+=("$2")
            shift 2
            ;;
        --c2c)
            C2C=true
            shift
            ;;
        perf|callgrind|auto)
            BACKEND="$1"
            shift
            ;;
        generate|align|infer|all)
            MODE="$1"
            shift
            ;;
        *)
            echo "Error: Unknown argument '$1'" >&2
            usage >&2
            exit 1
            ;;
    esac
done

# ------------------------------------------------------------------
# Backend selection
# ------------------------------------------------------------------
perf_usable() {
    if ! command -v perf &> /dev/null; then
        echo "perf: not found on PATH" >&2
        return 1
    fi
    local paranoid
    paranoid=$(cat /proc/sys/kernel/perf_event_paranoid 2>/dev/null || echo 999)
    if (( paranoid > 1 )); then
        echo "perf: unprivileged sampling blocked (kernel.perf_event_paranoid=$paranoid)" >&2
        echo "      Fix with one of:" >&2
        echo "        sudo sysctl -w kernel.perf_event_paranoid=1" >&2
        echo "        sudo -E pixi run profile -- perf $MODE" >&2
        return 1
    fi
    return 0
}

callgrind_usable() {
    command -v valgrind &> /dev/null
}

case "$BACKEND" in
    auto)
        if perf_usable; then
            BACKEND="perf"
        else
            echo "Falling back to callgrind." >&2
            if ! callgrind_usable; then
                echo "Error: valgrind not found on PATH either. Cannot profile." >&2
                exit 1
            fi
            BACKEND="callgrind"
        fi
        ;;
    perf)
        perf_usable || exit 1
        ;;
    callgrind)
        if ! callgrind_usable; then
            echo "Error: valgrind not found on PATH." >&2
            exit 1
        fi
        ;;
esac

if [[ ${#COLLECT_FNS[@]} -gt 0 && "$BACKEND" != "callgrind" ]]; then
    echo "Warning: --collect-fn is only supported with the callgrind backend; ignoring it for backend '$BACKEND'." >&2
    COLLECT_FNS=()
fi

if $C2C && [[ "$BACKEND" != "perf" ]]; then
    echo "Warning: --c2c is only supported with the perf backend; ignoring it for backend '$BACKEND'." >&2
    C2C=false
fi

REPORT_BASE_DIR="$REPORT_DIR/profile_${TIMESTAMP}_${BACKEND}"
mkdir -p "$REPORT_BASE_DIR"
LOG_FILE="$REPORT_BASE_DIR/profile.log"

echo "Backend: $BACKEND"
echo "Mode(s): $MODE"
echo "Threads: $N_THREADS"
echo "Report:  $REPORT_BASE_DIR"
echo ""

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------
run_setup() {
    local cmd="$1"
    echo "CMD: $cmd" >> "$LOG_FILE"
    eval "$cmd" >> "$LOG_FILE" 2>&1
}

profile_step() {
    local step_name="$1"
    local cmd="$2"
    local outdir="$REPORT_BASE_DIR/$step_name"
    mkdir -p "$outdir"

    printf "  %-12s : " "$step_name"
    echo "CMD: $cmd" >> "$LOG_FILE"

    local collect_flags=""
    if [[ ${#COLLECT_FNS[@]} -gt 0 ]]; then
        collect_flags="--collect-atstart=no"
        for fn in "${COLLECT_FNS[@]}"; do
            collect_flags="$collect_flags --toggle-collect=$fn"
        done
    fi

    local profile_cmd
    case "$BACKEND" in
        perf)
            # cycles for the call-graph, cache-references/cache-misses for
            # cache-miss statistics, collected in the same session.
            profile_cmd="perf record -F $PERF_FREQ -g --call-graph dwarf \
                         -e cycles:pp,cache-references,cache-misses \
                         -o ${outdir}/perf.data -- $cmd"
            ;;
        callgrind)
            # --cache-sim=yes adds simulated L1/LL cache-miss counts to the
            # same callgrind output (no extra file/session needed).
            profile_cmd="valgrind --tool=callgrind --cache-sim=yes $collect_flags \
                         --callgrind-out-file=${outdir}/callgrind.out -- $cmd"
            ;;
    esac

    if eval "$profile_cmd" >> "$LOG_FILE" 2>&1; then
        echo "OK"
    else
        echo "FAILED (see $LOG_FILE)"
        return 1
    fi

    case "$BACKEND" in
        perf)
            perf report --stdio -i "${outdir}/perf.data" > "${outdir}/report.txt" 2>> "$LOG_FILE"
            ;;
        callgrind)
            callgrind_annotate "${outdir}/callgrind.out" > "${outdir}/report.txt" 2>> "$LOG_FILE"
            ;;
    esac

    if [[ "$BACKEND" == "perf" ]]; then
        printf "  %-12s : " "$step_name (stat)"
        local stat_cmd="perf stat -e $STAT_EVENTS -o ${outdir}/perf_stat.txt -- $cmd"
        echo "CMD: $stat_cmd" >> "$LOG_FILE"
        if eval "$stat_cmd" >> "$LOG_FILE" 2>&1; then
            echo "OK"
        else
            echo "FAILED (see $LOG_FILE)"
            return 1
        fi
    fi

    if $C2C; then
        printf "  %-12s : " "$step_name (c2c)"
        mkdir -p "${outdir}/c2c"
        local c2c_cmd="perf c2c record -F $C2C_FREQ -o ${outdir}/c2c/perf.data -- $cmd"
        echo "CMD: $c2c_cmd" >> "$LOG_FILE"
        if eval "$c2c_cmd" >> "$LOG_FILE" 2>&1; then
            echo "OK"
        else
            echo "FAILED (see $LOG_FILE)"
            return 1
        fi
        perf c2c report --stdio -i "${outdir}/c2c/perf.data" > "${outdir}/c2c/report.txt" 2>> "$LOG_FILE"
    fi
}

# ------------------------------------------------------------------
# Setup (unprofiled): prepare model/genomic config and pipeline data so
# that "infer" always has alignment results to work with, regardless of
# which mode(s) are actually profiled.
# ------------------------------------------------------------------
run_setup "$IGORCALL init"
run_setup "$IGORCALL config set model.source custom"
run_setup "$IGORCALL config set model.parms \"$MODEL_PARMS\""
run_setup "$IGORCALL config set model.marginals \"$MODEL_MARGINALS\""
run_setup "$IGORCALL config set genomic.V \"$GENOMIC_V\""
run_setup "$IGORCALL config set genomic.D \"$GENOMIC_D\""
run_setup "$IGORCALL config set genomic.J \"$GENOMIC_J\""
run_setup "$IGORCALL config set generate.seed 12345"
run_setup "$IGORCALL config set infer.iterations 2"
run_setup "$IGORCALL -b $PFX -j 1 generate $N_SEQS"

SEQ_FILE="$PROF_DIR/${PFX}_generated/generated_seqs_werr.csv"
run_setup "$IGORCALL -b $PFX import-seqs \"$SEQ_FILE\""
run_setup "$IGORCALL -b $PFX -j $N_THREADS align --gene V"
run_setup "$IGORCALL -b $PFX -j $N_THREADS align --gene D"
run_setup "$IGORCALL -b $PFX -j $N_THREADS align --gene J"

# ------------------------------------------------------------------
# Profile each requested mode, each in its own session
# ------------------------------------------------------------------
run_generate() {
    # Separate batch so this doesn't disturb the pipeline data used by
    # the align/infer sessions.
    profile_step "generate" "$IGORCALL -b ${PFX}_gen -j $N_THREADS generate $N_SEQS"
}

run_align() {
    profile_step "align" "$IGORCALL -b $PFX -j $N_THREADS align --gene V && \
                           $IGORCALL -b $PFX -j $N_THREADS align --gene D && \
                           $IGORCALL -b $PFX -j $N_THREADS align --gene J"
}

run_infer() {
    profile_step "infer" "$IGORCALL -b $PFX -j $N_THREADS infer"
}

case "$MODE" in
    all)
        run_generate
        run_align
        run_infer
        ;;
    generate)
        run_generate
        ;;
    align)
        run_align
        ;;
    infer)
        run_infer
        ;;
esac

# ------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------
echo ""
echo "=================================================================="
echo "PROFILING COMPLETE"
echo "=================================================================="
echo "Report: $REPORT_BASE_DIR"
echo "Log:    $LOG_FILE"
echo ""
case "$BACKEND" in
    perf)
        echo "Explore interactively: cd $REPORT_BASE_DIR/<mode> && perf report"
        echo "(cache-references/cache-misses were recorded alongside cycles)"
        echo "Exact aggregate counters: $REPORT_BASE_DIR/<mode>/perf_stat.txt"
        if $C2C; then
            echo "False sharing: cd $REPORT_BASE_DIR/<mode>/c2c && perf c2c report"
        fi
        echo "Flamegraph (needs https://github.com/brendangregg/FlameGraph on PATH):"
        echo "  cd $REPORT_BASE_DIR/<mode> && perf script | stackcollapse-perf.pl | flamegraph.pl > flame.svg"
        ;;
    callgrind)
        echo "Explore interactively: kcachegrind $REPORT_BASE_DIR/<mode>/callgrind.out"
        echo "Or read the text summary: $REPORT_BASE_DIR/<mode>/report.txt"
        echo "(L1/LL cache-miss simulation was included via --cache-sim=yes)"
        ;;
esac

exit 0
