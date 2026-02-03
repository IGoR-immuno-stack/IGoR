# IGoR Fast Generation Benchmark Results

## Summary

The optimized fast generation implementation provides significant speedups over the original method, especially for large batch sizes.

## Performance Comparison (Human TCR Beta)

| Sequences | Original | Fast | Speedup | Fast Rate (seq/s) |
|-----------|----------|------|---------|-------------------|
| 1,000 | 0.053s | 0.034s | 1.56x | 29,221 |
| 10,000 | 0.173s | 0.038s | 4.55x | 262,049 |
| 100,000 | 1.459s | 0.127s | 11.49x | 787,843 |
| 1,000,000 | 14.503s | 0.727s | 19.95x | 1,375,708 |
| 10,000,000 | 151.366s | 5.746s | **26.34x** | **1,740,360** |

## Key Optimizations

1. **O(log n) Sampling**: Binary search on precomputed CDFs instead of linear search through realization maps
2. **Cached Data Structures**: Gene sequences and deletion/insertion values are cached in vectors for O(1) indexed access
3. **Parallel Generation**: Multiple threads generate sequences concurrently with thread-local RNG and buffers
4. **Buffered I/O**: Output is batched and written in chunks to reduce I/O overhead
5. **Reduced Allocations**: Reusable RealizationBuffer objects avoid per-sequence memory allocations

## Usage

### Command Line
```bash
# Fast generation with auto-detected threads
igor -set_wd /output -batch test -set_custom_model model_parms.txt model_marginals.txt \
    -set_genomic --V genomicVs.fasta --D genomicDs.fasta --J genomicJs.fasta \
    -generate 10000000 --fast

# Specify thread count
igor ... -generate 10000000 --fast --threads 8

# Sequences only (no realizations output)
igor ... -generate 10000000 --fast --sequences-only

# Reproducible with seed
igor ... -generate 10000000 --fast --seed 42
```

### Programmatic (C++)
```cpp
igor::GenerationConfig config;
config.num_threads = 8;
config.generate_errors = true;
config.seed = 42;

genmodel.generate_sequences_fast(10000000, "sequences.csv", "realizations.csv", config);

auto stats = genmodel.get_last_generation_stats();
std::cout << "Rate: " << stats.sequences_per_second << " seq/s" << std::endl;
```

## Running Benchmarks

```bash
# Full benchmark suite
./scripts/benchmark_generation.sh

# Sampling microbenchmarks
./scripts/benchmark_sampling_only.sh
```

## Test Results

All unit tests pass (2172 assertions):
- CategoricalSampler correctness
- ConditionalSampler correctness
- DinucleotideMarkovSampler correctness
- RealizationBuffer operations
- CDF normalization
- Binary search correctness
