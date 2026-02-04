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
