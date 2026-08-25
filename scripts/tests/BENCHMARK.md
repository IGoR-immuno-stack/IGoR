# IGoR Scalability Benchmarks

This document describes the scalability benchmarks implemented for IGoR.

## Benchmark Script

The benchmark logic is encapsulated in `scripts/tests/benchmark_run.sh`.
Results are logged to a temporary file (path printed at start) and a summary table is printed to stdout.

## Usages

### 1. Default (All Standard Tests)
Runs both the standard pipeline benchmarks and the default standalone sampling tests.
```bash
pixi run benchmark
```

### 2. Standalone Sampling (Generation)
Tests the throughput of the generator.

**Default Sizes**:
```bash
pixi run benchmark sampling
```

**Custom / Stress Test**:
Run with specific sequence counts. You can pass multiple values to see scaling behavior.
```bash
# Example: Scale from 1k to 1 million
pixi run benchmark -- sampling 1000 10000 100000 1000000
```

### 3. Pipeline Benchmarks
Runs the full workflow (Generate $\to$ Read $\to$ Align $\to$ Infer) to test end-to-end performance and alignment/inference scaling.
```bash
pixi run benchmark pipeline
```

## Methodology

### Metrics Measured
We measure **Wall Clock Time** (in seconds).

#### Standalone
*   **Generation (Standalone)**: Pure synthetic sequence generation throughput. Isolated from pipeline steps.

#### Pipeline Steps
*   **Pipeline: Generation**: Generation step specifically to prepare input for the pipeline tests.
*   **Pipeline: Read**: Parsing sequences.
*   **Pipeline: Alignment**: Alignment against V, D, and J genomic templates.
    *   *Verbose Mode*: Splits into **Align V**, **Align D**, **Align J**.
*   **Pipeline: Inference**: EM algorithm inference.

### Test Cases

| Category | Sequences ($N$) | Details |
| :--- | :--- | :--- |
| **Standalone Sampling** | 10k, 100k, (Custom) | Throughput testing. Can be scaled to billions via arguments. |
| **Full Pipeline** | 100, 500, 1000 | Comprehensive test. restricted sizes due to Alignment intensity. |

## Profiling

Beyond timing, `scripts/tests/profile_run.sh` profiles the `generate`,
`align`, and `infer` subcommands individually (each in its own profiling
session) on a fixed single-threaded, N=100 pipeline run, to find hotspots
worth optimizing.

```bash
pixi run profile                       # auto-detect backend, profile all 3 modes
pixi run profile -- perf generate      # force perf, profile generation only
pixi run profile -- callgrind infer    # force callgrind, profile inference only
```

Two backends are supported:
*   **`perf`**: low-overhead sampling, best signal, but needs the system
    `perf` binary and unprivileged `perf_event_open` access. If
    `/proc/sys/kernel/perf_event_paranoid` is greater than `1`, unprivileged
    sampling is blocked; either run
    `sudo sysctl -w kernel.perf_event_paranoid=1` once, or run the task
    under `sudo -E`.
*   **`callgrind`** (valgrind): no special privileges needed, works
    identically on any Linux box — the portable fallback.

With `BACKEND=auto` (the default), the script uses `perf` if it's usable and
transparently falls back to `callgrind` otherwise.

Reports are written to `profile_reports/profile_<timestamp>_<backend>/<mode>/`
— one directory per mode, containing the raw profile at its tool's default
filename (`perf.data` or `callgrind.out`, so it can be opened without `-i`)
and a text summary (`report.txt`). Explore interactively with
`cd .../<mode> && perf report` or `kcachegrind callgrind.out`.

Cache-miss statistics are collected by default with both backends: `perf`
records `cache-references`/`cache-misses` alongside `cycles` in the same
session, and `callgrind` runs with `--cache-sim=yes` to simulate L1/LL
cache misses in the same output.

### `perf record` vs `perf stat`

`perf record` (used for `<mode>/perf.data`) is *sampling*: it periodically
interrupts the program and tags the current instruction pointer/call stack,
which is what gives per-function attribution, but the totals are a
statistical subsample, not exact counts. `perf stat` is *counting*: it
reads the hardware counters directly before/after the whole run, giving
exact totals with no per-function breakdown. The two aren't
interchangeable — you can't recover `perf stat`'s exact numbers from a
`perf.data` sample capture after the fact.

For that reason, with the `perf` backend each mode also gets a `perf stat`
pass at `<mode>/perf_stat.txt`, covering `cycles`, `instructions`,
`L1-dcache-loads`/`L1-dcache-load-misses`, `LLC-loads`/`LLC-load-misses`,
`branches`/`branch-misses`, and `stalled-cycles-frontend`/`backend` — enough
to compute IPC, cache/branch miss rates, and see whether stalls are
frontend- (fetch/decode) or backend- (execution/memory) bound. Hardware
typically exposes only 4-8 general-purpose counters per core, so requesting
this many events causes `perf stat` to time-multiplex them; it reports
scaled estimates in that case (noted in the output). On hybrid CPUs (see
below) `perf stat` also splits per-PMU like `perf report` does, and some
events may show `<not supported>` on the `cpu_atom` (E-core) PMU.

### Reading `perf report`'s "Available samples" list

On hybrid Intel CPUs (12th‑gen+, e.g. Alder Lake/Raptor Lake and newer —
anything with performance + efficiency cores), `perf report` will present
each event twice, once per PMU, and ask which one to display:

```
926 cpu_atom/cycles/
10K cpu_core/cycles/
322 cpu_atom/cache-references/
9K  cpu_core/cache-references/
164 cpu_atom/cache-misses/
2K  cpu_core/cache-misses/
```

*   **`cpu_core`** is the PMU for the performance ("P") cores; **`cpu_atom`**
    is the PMU for the efficiency ("E") cores (named `atom` because their
    microarchitecture descends from Intel's Atom line). The kernel schedules
    threads across both core types, so a single run's samples get split
    between the two PMUs depending on which physical core each sample landed
    on — it's not two different metrics, just the same metric bucketed by
    where it was recorded. If a run's samples are lopsided toward
    `cpu_atom` (e.g. the process got scheduled onto E-cores), timings won't
    be comparable to a run that landed mostly on `cpu_core`; pin with
    `taskset` to one core type for consistent, repeatable profiles if this
    matters.
*   **`cycles`**: CPU clock cycles consumed while running — the base signal
    for "where is time being spent," independent of clock speed variation.
*   **`cache-references`**: memory accesses that queried the cache
    hierarchy (the PMU's generic last-level-cache reference event) — the
    denominator for a miss rate.
*   **`cache-misses`**: the subset of `cache-references` that weren't
    satisfied by cache and had to go to slower memory. A high
    `cache-misses`/`cache-references` ratio in a function flags poor data
    locality (bad access patterns, or false sharing across threads) — a
    common, distinct-from-`cycles` reason a function is slow.

### Scoping callgrind to a specific function

Each profiled subcommand invocation includes CLI parsing and
config/model/genomic-template loading, which can dominate a callgrind
profile of a small (N=100) run. Use `--collect-fn=FUNCTION` (repeatable) to
disable collection at process start and only toggle it on while inside
`FUNCTION` (and its callees) — i.e. `--collect-atstart=no` +
`--toggle-collect=FUNCTION` under the hood. Only supported with the
`callgrind` backend; ignored (with a warning) otherwise.

```bash
pixi run profile -- callgrind align --collect-fn=Aligner::align_seqs
```

### Finding OpenMP false sharing with perf c2c

`--c2c` (perf backend only) additionally records a `perf c2c`
(cache-to-cache) session per mode, written to `<mode>/c2c/perf.data` with
its report at `<mode>/c2c/report.txt` — useful for spotting false sharing
between OpenMP threads. It's meaningless with the default single-threaded
target, so pair it with `--threads`:

```bash
pixi run profile -- perf align --threads=4 --c2c
```
