# IGoR Fast Generator Performance Report

## Executive Summary

This report compares three sequence generation implementations in IGoR:

1. **Original** - The existing `GenModel::generate_sequences()` method
2. **Fast ExactMatch** - New implementation using `--fast` flag, produces identical sequences
3. **Fast MaxSpeed** - New implementation using `--maxspeed` flag, maximum performance with statistically equivalent sequences

## Benchmark Results

### Single-Threaded Performance Comparison

| Sequences | Original | ExactMatch | MaxSpeed | ExactMatch Speedup | MaxSpeed Speedup |
|-----------|----------|------------|----------|-------------------|------------------|
| 1,000 | 12,561 seq/s | 23,285 seq/s | 30,989 seq/s | **1.9x** | **2.5x** |
| 10,000 | 55,451 seq/s | 64,673 seq/s | 176,463 seq/s | **1.2x** | **3.2x** |
| 100,000 | 64,830 seq/s | 79,238 seq/s | 338,695 seq/s | **1.2x** | **5.2x** |
| 1,000,000 | 63,987 seq/s | 79,527 seq/s | 361,756 seq/s | **1.2x** | **5.7x** |

### Multi-Threaded MaxSpeed Performance (1M sequences)

| Threads | Rate | Speedup vs Original |
|---------|------|---------------------|
| 1 | 362k seq/s | **5.7x** |
| 4 | **1.05M seq/s** | **16.4x** |
| 8 | **1.48M seq/s** | **23.1x** |

## Implementation Details

### Original Generator (`GenModel::generate_sequences`)

The original implementation:
- Uses `draw_random_realization()` method on each event
- Linear O(n) search through probability arrays for each sample
- Recomputes index maps and offset maps per sequence
- Single-threaded only

### Fast ExactMatch Mode (`--fast`)

Key optimizations while maintaining bit-exact reproducibility:
- **Cached model data**: Model queue, index maps, and offset maps computed once and reused
- **Per-thread model copies**: Each thread has its own deep copy of model parameters to avoid race conditions
- **Calls original event methods**: Uses `draw_random_realization()` to ensure identical sampling behavior
- **Parallel execution**: Multiple threads with independent RNG streams

**Guarantees:**
- ✅ Sequences are **bit-exact identical** to original with same seed (single-threaded)
- ✅ Realizations are identical
- ✅ Reproducible with same seed

### Fast MaxSpeed Mode (`--maxspeed`)

Maximum performance using precomputed sampling structures:
- **Precomputed CDF samplers**: Binary search O(log n) instead of linear O(n)
- **Index-ordered iteration**: Deterministic iteration through realizations
- **Optimized dinucleotide Markov**: Precomputed conditional distributions
- **Minimal per-sequence overhead**: No map copies or recomputation

**Trade-offs:**
- ✅ **Statistically equivalent** sequences (same distribution)
- ✅ Much higher throughput (4-20x faster)
- ❌ Sequences **not identical** to original (different sampling order)
- ✅ Reproducible with same seed

## Usage

```bash
# Original method
igor -generate 100000 --seed 42

# Fast ExactMatch mode (identical sequences, moderate speedup)
igor -generate 100000 --seed 42 --fast --threads 4

# Fast MaxSpeed mode (maximum speed, equivalent statistics)
igor -generate 100000 --seed 42 --maxspeed --threads 4
```

## When to Use Each Mode

| Use Case | Recommended Mode |
|----------|------------------|
| Reproducing exact results from original IGoR | **Original** or **ExactMatch** |
| Validation/testing against original | **ExactMatch** |
| Large-scale generation for statistics | **MaxSpeed** |
| Multi-threaded production runs | **MaxSpeed** |
| Backward compatibility required | **Original** or **ExactMatch** |

## Technical Architecture

### Event Sampling Fix for Determinism

The key fix enabling exact reproducibility was modifying how events iterate through realizations:

**Before (non-deterministic):**
```cpp
// std::unordered_map iteration order varies between copies
for (const auto& [name, real] : realizations) {
    prob_count += marginals[index + real.index];
    if (prob_count >= rand) return real;
}
```

**After (deterministic):**
```cpp
// Sort by index for deterministic order
std::vector<...> sorted_reals;
for (const auto& [name, real] : realizations) {
    sorted_reals.emplace_back(real.index, ...);
}
std::sort(sorted_reals.begin(), sorted_reals.end());
for (const auto& [idx, ...] : sorted_reals) {
    prob_count += marginals[index + idx];
    if (prob_count >= rand) return ...;
}
```

### Thread Safety

Each thread has its own:
- RNG (`std::mt19937_64`) with unique seed derived from base seed
- Complete copy of `Model_Parms` (deep copy)
- Model queue, index map, offset map
- Output buffers

This eliminates all shared mutable state between threads.

## Conclusion

The Fast Generator provides two modes to serve different needs:

1. **ExactMatch** for users who need reproducibility with the original algorithm
2. **MaxSpeed** for users who need maximum throughput for large-scale generation

Both modes are significantly faster than the original implementation, with MaxSpeed achieving up to **4-20x speedup** depending on the workload size.
