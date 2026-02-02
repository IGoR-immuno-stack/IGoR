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
Run with a specific sequence count (e.g., 1,000,000). To go huge (e.g., 1,000,000,000), typical limits apply (disk space, time).
```bash
# Example: 1 million sequences
pixi run benchmark -- sampling 1000000
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
*   **Standalone_Sampling**: Pure synthetic sequence generation throughput. Isolated from pipeline steps.

#### Pipeline Steps
*   **Pipeline_Gen**: Generation step specifically to prepare input for the pipeline tests.
*   **Pipeline_Read**: Parsing sequences.
*   **Pipeline_Align**: Alignment against V, D, and J genomic templates.
*   **Pipeline_Infer**: EM algorithm inference.

### Test Cases

| Category | Sequences ($N$) | Details |
| :--- | :--- | :--- |
| **Standalone Sampling** | 10k, 100k, (Custom) | Throughput testing. Can be scaled to billions via arguments. |
| **Full Pipeline** | 100, 500, 1000 | Comprehensive test. restricted sizes due to Alignment intensity. |
