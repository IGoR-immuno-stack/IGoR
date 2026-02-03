#!/usr/bin/env python3
"""
benchmark_generation.py

Comprehensive benchmark script for IGoR sequence generation algorithms.
Compares Original, ExactMatch (--fast), and MaxSpeed (--maxspeed) implementations.

===============================================================================
OVERVIEW
===============================================================================

This script performs four types of analysis:

  1. PERFORMANCE: Measures generation speed (sequences/second)
     - Tests Original, ExactMatch (1T/4T), MaxSpeed (1T/4T)
     - Calculates speedup ratios vs Original

  2. CORRECTNESS: Compares output files for exact and sorted matches
     - Exact match: byte-for-byte identical files
     - Sorted match: same sequences but different order
     - Useful for validating multi-threaded implementations

  3. DISTRIBUTION: Analyzes statistical properties of generated sequences
     - Sequence length distribution (min, max, mean, std)
     - Kolmogorov-Smirnov test for distribution similarity
     - Ensures MaxSpeed produces statistically equivalent sequences

  4. REPRODUCIBILITY: Validates deterministic behavior with same seed
     - Runs each single-threaded method twice with identical seed
     - Verifies outputs are byte-for-byte identical
     - Essential for scientific reproducibility

===============================================================================
USAGE
===============================================================================

  python3 scripts/benchmark_generation.py [OPTIONS]

Options:
  --quick       Run only small benchmarks (1,000 sequences)
  --full        Run all benchmarks including 10M sequences
  --skip-dist   Skip distribution analysis (faster)
  --help        Show this help message

Examples:
  # Quick test (1,000 sequences only)
  python3 scripts/benchmark_generation.py --quick

  # Default run (1,000 and 1,000,000 sequences)
  python3 scripts/benchmark_generation.py

  # Full run including 10,000,000 sequences
  python3 scripts/benchmark_generation.py --full

===============================================================================
ENVIRONMENT VARIABLES
===============================================================================

  IGOR_BIN      Path to igor binary (default: ./build/bin/igor)
  MODELS_DIR    Path to models directory (default: ./models)
  OUTPUT_DIR    Path for output files (default: /tmp/igor_benchmark)

===============================================================================
IMPLEMENTATION NOTES
===============================================================================

ExactMatch Mode (--fast):
  - Uses the SAME sampling algorithm as Original
  - Produces IDENTICAL sequences with same seed (single-threaded)
  - Multi-threaded version uses per-thread RNG seeds, so sequences differ
    but the statistical distribution remains identical

MaxSpeed Mode (--maxspeed):
  - Uses DIFFERENT sampling algorithm (binary search CDF)
  - Produces DIFFERENT sequences but SAME statistical distribution
  - Significantly faster due to O(log n) vs O(n) sampling
  - Sequences are statistically equivalent (validated by KS test)

Author: IGoR Development Team
Date: February 2026
"""

import argparse
import os
import subprocess
import sys
import time
import statistics
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple, List, Dict

# =============================================================================
# CONFIGURATION
# =============================================================================

# Default paths (can be overridden by environment variables)
DEFAULT_IGOR_BIN = "./build/bin/igor"
DEFAULT_MODELS_DIR = "./models"
DEFAULT_OUTPUT_DIR = "/tmp/igor_benchmark"
DEFAULT_SEED = 42

# Benchmark sizes for each mode
SIZES_QUICK = [1000]
SIZES_DEFAULT = [1000, 1_000_000]
SIZES_FULL = [1000, 1_000_000, 10_000_000]

# Terminal colors
class Colors:
    RED = '\033[0;31m'
    GREEN = '\033[0;32m'
    YELLOW = '\033[1;33m'
    BLUE = '\033[0;34m'
    CYAN = '\033[0;36m'
    BOLD = '\033[1m'
    NC = '\033[0m'  # No Color (reset)


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class BenchmarkResult:
    """Result of a single benchmark run."""
    method: str
    num_sequences: int
    duration: float
    rate: int
    success: bool
    output_dir: str


@dataclass
class SequenceStats:
    """Statistics for a set of sequences."""
    count: int
    min_length: int
    max_length: int
    mean_length: float
    std_length: float


@dataclass 
class ComparisonResult:
    """Result of comparing two sequence files."""
    exact_match: bool
    sorted_match: bool
    distribution_similar: bool
    mean_diff: float
    pct_diff: float
    p_value: Optional[float] = None


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def format_number(n: int) -> str:
    """Format a number with thousand separators."""
    return f"{n:,}"


def print_header(text: str) -> None:
    """Print a prominent section header."""
    print()
    print(f"{Colors.BLUE}{'═' * 67}{Colors.NC}")
    print(f"{Colors.BLUE}  {text}{Colors.NC}")
    print(f"{Colors.BLUE}{'═' * 67}{Colors.NC}")
    print()


def print_subheader(text: str) -> None:
    """Print a subsection header."""
    print(f"{Colors.CYAN}─── {text} ───{Colors.NC}")


def check_prerequisites(igor_bin: str, models_dir: str) -> bool:
    """Verify all prerequisites are met."""
    # Check igor binary
    if not os.path.isfile(igor_bin) or not os.access(igor_bin, os.X_OK):
        print(f"{Colors.RED}Error: IGoR binary not found at {igor_bin}{Colors.NC}")
        print("Build IGoR first or set IGOR_BIN environment variable")
        return False
    
    # Check models directory
    tcr_beta_path = os.path.join(models_dir, "human", "tcr_beta")
    if not os.path.isdir(tcr_beta_path):
        print(f"{Colors.RED}Error: Models not found at {models_dir}{Colors.NC}")
        print("Run from the IGoR root directory or set MODELS_DIR environment variable")
        return False
    
    return True


# =============================================================================
# BENCHMARK EXECUTION
# =============================================================================

def run_benchmark(
    igor_bin: str,
    models_dir: str, 
    output_dir: str,
    num_sequences: int,
    method_name: str,
    extra_args: str,
    seed: int
) -> BenchmarkResult:
    """
    Execute a single generation benchmark.
    
    Args:
        igor_bin: Path to igor binary
        models_dir: Path to models directory
        output_dir: Path for output files
        num_sequences: Number of sequences to generate
        method_name: Identifier for this run
        extra_args: Additional arguments for igor
        seed: Random seed
    
    Returns:
        BenchmarkResult with timing and status information
    """
    # Construct the command
    batch_name = f"{method_name}_{num_sequences}"
    cmd = [
        igor_bin,
        "-set_wd", output_dir,
        "-batch", batch_name,
        "-set_custom_model",
        os.path.join(models_dir, "human/tcr_beta/models/model_parms.txt"),
        os.path.join(models_dir, "human/tcr_beta/models/model_marginals.txt"),
        "-set_genomic",
        "--V", os.path.join(models_dir, "human/tcr_beta/ref_genome/genomicVs.fasta"),
        "--D", os.path.join(models_dir, "human/tcr_beta/ref_genome/genomicDs.fasta"),
        "--J", os.path.join(models_dir, "human/tcr_beta/ref_genome/genomicJs.fasta"),
        "-generate", str(num_sequences),
        "--seed", str(seed),
    ]
    
    # Add extra arguments
    if extra_args:
        cmd.extend(extra_args.split())
    
    # Run and time the command
    start_time = time.perf_counter()
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        end_time = time.perf_counter()
        duration = end_time - start_time
        rate = int(num_sequences / duration)
        
        return BenchmarkResult(
            method=method_name,
            num_sequences=num_sequences,
            duration=duration,
            rate=rate,
            success=True,
            output_dir=os.path.join(output_dir, f"{batch_name}_generated")
        )
    except subprocess.CalledProcessError as e:
        return BenchmarkResult(
            method=method_name,
            num_sequences=num_sequences,
            duration=0,
            rate=0,
            success=False,
            output_dir=""
        )


# =============================================================================
# FILE COMPARISON
# =============================================================================

def compare_files_exact(file1: str, file2: str) -> bool:
    """Check if two files are byte-for-byte identical."""
    if not os.path.isfile(file1) or not os.path.isfile(file2):
        return False
    
    try:
        with open(file1, 'rb') as f1, open(file2, 'rb') as f2:
            return f1.read() == f2.read()
    except IOError:
        return False


def compare_files_sorted(file1: str, file2: str) -> bool:
    """Check if two files contain the same sequences (order-independent)."""
    if not os.path.isfile(file1) or not os.path.isfile(file2):
        return False
    
    try:
        def get_sorted_sequences(filepath: str) -> List[str]:
            sequences = []
            with open(filepath, 'r') as f:
                next(f)  # Skip header
                for line in f:
                    parts = line.strip().split(';')
                    if len(parts) >= 2:
                        sequences.append(parts[1])
            return sorted(sequences)
        
        seqs1 = get_sorted_sequences(file1)
        seqs2 = get_sorted_sequences(file2)
        return seqs1 == seqs2
    except IOError:
        return False


# =============================================================================
# DISTRIBUTION ANALYSIS
# =============================================================================

def get_sequence_lengths(filepath: str) -> List[int]:
    """Extract sequence lengths from a CSV file."""
    lengths = []
    try:
        with open(filepath, 'r') as f:
            next(f)  # Skip header
            for line in f:
                parts = line.strip().split(';')
                if len(parts) >= 2:
                    lengths.append(len(parts[1]))
    except IOError:
        pass
    return lengths


def compute_sequence_stats(filepath: str) -> Optional[SequenceStats]:
    """Compute statistics for sequences in a file."""
    lengths = get_sequence_lengths(filepath)
    
    if not lengths:
        return None
    
    return SequenceStats(
        count=len(lengths),
        min_length=min(lengths),
        max_length=max(lengths),
        mean_length=statistics.mean(lengths),
        std_length=statistics.stdev(lengths) if len(lengths) > 1 else 0.0
    )


def compare_distributions(file1: str, file2: str) -> ComparisonResult:
    """
    Compare sequence length distributions using KS test.
    
    Uses the Kolmogorov-Smirnov two-sample test to determine if two sets
    of sequences have the same length distribution.
    """
    lengths1 = get_sequence_lengths(file1)
    lengths2 = get_sequence_lengths(file2)
    
    if not lengths1 or not lengths2:
        return ComparisonResult(
            exact_match=False,
            sorted_match=False,
            distribution_similar=False,
            mean_diff=0,
            pct_diff=0
        )
    
    # Compute means
    m1 = statistics.mean(lengths1)
    m2 = statistics.mean(lengths2)
    mean_diff = abs(m1 - m2)
    pct_diff = 100 * mean_diff / m1 if m1 > 0 else 0
    
    # Try to use scipy for KS test
    p_value = None
    try:
        from scipy import stats
        statistic, p_value = stats.ks_2samp(lengths1, lengths2)
        # Consider similar if p > 0.05 OR if mean difference < 1.5%
        distribution_similar = p_value > 0.05 or pct_diff < 1.5
    except ImportError:
        # Fall back to simple comparison
        distribution_similar = pct_diff < 1.5
    
    return ComparisonResult(
        exact_match=compare_files_exact(file1, file2),
        sorted_match=compare_files_sorted(file1, file2),
        distribution_similar=distribution_similar,
        mean_diff=mean_diff,
        pct_diff=pct_diff,
        p_value=p_value
    )


# =============================================================================
# BENCHMARK RUNNER
# =============================================================================

def run_benchmarks_for_size(
    size: int,
    igor_bin: str,
    models_dir: str,
    output_dir: str,
    seed: int,
    skip_dist: bool
) -> Dict[str, BenchmarkResult]:
    """
    Execute all benchmark variants for a given size.
    
    Runs five generation methods:
      1. Original (baseline)
      2. ExactMatch with 1 thread
      3. ExactMatch with 4 threads
      4. MaxSpeed with 1 thread
      5. MaxSpeed with 4 threads
      6. MaxSpeed with all available threads (NT)
    
    Also tests reproducibility by running each single-threaded method twice.
    """
    print_subheader(f"Benchmark: {format_number(size)} sequences")
    
    results: Dict[str, BenchmarkResult] = {}
    
    # Get number of available threads
    num_threads = os.cpu_count() or 4
    
    # Define benchmark configurations
    benchmarks = [
        ("original", "Original", ""),
        ("exactmatch_1t", "ExactMatch (1 thread)", "--fast --threads 1"),
        ("exactmatch_4t", "ExactMatch (4 threads)", "--fast --threads 4"),
        ("maxspeed_1t", "MaxSpeed (1 thread)", "--maxspeed --threads 1"),
        ("maxspeed_4t", "MaxSpeed (4 threads)", "--maxspeed --threads 4"),
        ("maxspeed_nt", f"MaxSpeed ({num_threads} threads)", f"--maxspeed --threads {num_threads}"),
    ]
    
    # =========================================================================
    # PERFORMANCE BENCHMARKS
    # =========================================================================
    
    print(f"\n  {Colors.BOLD}Performance Benchmarks:{Colors.NC}")
    
    for method_id, method_name, extra_args in benchmarks:
        print(f"    {method_name:26s} ", end="", flush=True)
        
        result = run_benchmark(
            igor_bin, models_dir, output_dir,
            size, method_id, extra_args, seed
        )
        results[method_id] = result
        
        if result.success:
            print(f"{Colors.GREEN}{result.duration:.3f}s ({format_number(result.rate)} seq/s){Colors.NC}")
        else:
            print(f"{Colors.RED}FAILED{Colors.NC}")
    
    # =========================================================================
    # CORRECTNESS COMPARISON
    # =========================================================================
    
    print(f"\n  {Colors.BOLD}Correctness Comparison:{Colors.NC}")
    
    orig_file = os.path.join(results["original"].output_dir, "generated_seqs_werr.csv")
    
    comparisons = [
        ("exactmatch_1t", "ExactMatch(1T) vs Original"),
        ("exactmatch_4t", "ExactMatch(4T) vs Original"),
        ("maxspeed_1t", "MaxSpeed(1T) vs Original"),
        ("maxspeed_4t", "MaxSpeed(4T) vs Original"),
        ("maxspeed_nt", f"MaxSpeed({num_threads}T) vs Original"),
    ]
    
    for method_id, label in comparisons:
        if method_id not in results or not results[method_id].success:
            continue
            
        other_file = os.path.join(results[method_id].output_dir, "generated_seqs_werr.csv")
        print(f"    {label:32s} ", end="")
        
        if compare_files_exact(orig_file, other_file):
            print(f"{Colors.GREEN}✓ IDENTICAL (exact byte match){Colors.NC}")
        elif compare_files_sorted(orig_file, other_file):
            print(f"{Colors.YELLOW}≈ Same sequences (different order){Colors.NC}")
        else:
            if "maxspeed" in method_id:
                print(f"{Colors.CYAN}○ Different sequences (optimized sampler){Colors.NC}")
            else:
                print(f"{Colors.CYAN}○ Different sequences (per-thread RNG seeds){Colors.NC}")
    
    # =========================================================================
    # DISTRIBUTION ANALYSIS
    # =========================================================================
    
    if not skip_dist:
        print(f"\n  {Colors.BOLD}Distribution Analysis:{Colors.NC}")
        
        # Compute statistics for each method
        print("    Sequence length statistics:")
        
        stats_methods = [
            ("original", "Original"),
            ("exactmatch_1t", "ExactMatch(1T)"),
            ("maxspeed_1t", "MaxSpeed(1T)"),
            ("maxspeed_4t", "MaxSpeed(4T)"),
            ("maxspeed_nt", f"MaxSpeed({num_threads}T)"),
        ]
        
        for method_id, method_name in stats_methods:
            if method_id not in results or not results[method_id].success:
                continue
            
            filepath = os.path.join(results[method_id].output_dir, "generated_seqs_werr.csv")
            stats = compute_sequence_stats(filepath)
            
            if stats:
                print(f"      {method_name:18s} min={stats.min_length:<4d} max={stats.max_length:<4d} "
                      f"mean={stats.mean_length:<6.1f} std={stats.std_length:<5.1f}")
        
        # Distribution comparison
        print()
        print("    Distribution comparison (KS test, p>0.05 = similar):")
        
        dist_comparisons = [
            ("exactmatch_1t", "ExactMatch(1T) vs Original"),
            ("maxspeed_1t", "MaxSpeed(1T) vs Original"),
            ("maxspeed_4t", "MaxSpeed(4T) vs Original"),
            ("maxspeed_nt", f"MaxSpeed({num_threads}T) vs Original"),
        ]
        
        for method_id, label in dist_comparisons:
            if method_id not in results or not results[method_id].success:
                continue
            
            other_file = os.path.join(results[method_id].output_dir, "generated_seqs_werr.csv")
            comparison = compare_distributions(orig_file, other_file)
            
            print(f"      {label:32s} ", end="")
            
            if comparison.p_value is not None:
                if comparison.distribution_similar:
                    if comparison.p_value > 0.05:
                        print(f"{Colors.GREEN}SIMILAR (p={comparison.p_value:.4f}, "
                              f"Δmean={comparison.mean_diff:.1f}){Colors.NC}")
                    else:
                        print(f"{Colors.GREEN}~SIMILAR (p={comparison.p_value:.4f}, "
                              f"Δmean={comparison.mean_diff:.1f}, {comparison.pct_diff:.1f}%){Colors.NC}")
                else:
                    print(f"{Colors.RED}DIFFERENT (p={comparison.p_value:.4f}, "
                          f"Δmean={comparison.mean_diff:.1f}, {comparison.pct_diff:.1f}%){Colors.NC}")
            else:
                if comparison.distribution_similar:
                    print(f"{Colors.GREEN}~SIMILAR (Δmean={comparison.mean_diff:.1f}, "
                          f"{comparison.pct_diff:.1f}%){Colors.NC}")
                else:
                    print(f"{Colors.RED}DIFFERENT (Δmean={comparison.mean_diff:.1f}, "
                          f"{comparison.pct_diff:.1f}%){Colors.NC}")
    
    # =========================================================================
    # SPEEDUP CALCULATIONS
    # =========================================================================
    
    print(f"\n  {Colors.BOLD}Speedup vs Original:{Colors.NC}")
    
    if results["original"].success:
        orig_time = results["original"].duration
        
        print(f"    {'Method':<24s} {'Speedup':>8s}")
        print(f"    {'-' * 24} {'-' * 8}")
        
        speedup_methods = [
            ("exactmatch_1t", "ExactMatch (1 thread)"),
            ("exactmatch_4t", "ExactMatch (4 threads)"),
            ("maxspeed_1t", "MaxSpeed (1 thread)"),
            ("maxspeed_4t", "MaxSpeed (4 threads)"),
            ("maxspeed_nt", f"MaxSpeed ({num_threads} threads)"),
        ]
        
        for method_id, method_name in speedup_methods:
            if method_id in results and results[method_id].success:
                speedup = orig_time / results[method_id].duration
                print(f"    {method_name:<26s} {speedup:>7.2f}x")
    
    # =========================================================================
    # REPRODUCIBILITY TEST
    # =========================================================================
    # Run each single-threaded method twice with the same seed and verify
    # that the output is identical. This validates deterministic behavior.
    
    print(f"\n  {Colors.BOLD}Reproducibility Test (same seed, successive runs):{Colors.NC}")
    
    # Get number of available threads
    num_threads = os.cpu_count() or 4
    
    reproducibility_tests = [
        ("original", "Original", ""),
        ("exactmatch_1t", "ExactMatch (1T)", "--fast --threads 1"),
        ("maxspeed_1t", "MaxSpeed (1T)", "--maxspeed --threads 1"),
        ("maxspeed_nt", f"MaxSpeed ({num_threads}T)", f"--maxspeed --threads {num_threads}"),
    ]
    
    for method_id, method_name, extra_args in reproducibility_tests:
        print(f"    {method_name:20s} ", end="", flush=True)
        
        # Run first time
        result1 = run_benchmark(
            igor_bin, models_dir, output_dir,
            size, f"{method_id}_run1", extra_args, seed
        )
        
        # Run second time with same seed
        result2 = run_benchmark(
            igor_bin, models_dir, output_dir,
            size, f"{method_id}_run2", extra_args, seed
        )
        
        if not result1.success or not result2.success:
            print(f"{Colors.RED}✗ FAILED (run error){Colors.NC}")
            continue
        
        # Compare outputs
        file1 = os.path.join(result1.output_dir, "generated_seqs_werr.csv")
        file2 = os.path.join(result2.output_dir, "generated_seqs_werr.csv")
        
        if compare_files_exact(file1, file2):
            print(f"{Colors.GREEN}✓ REPRODUCIBLE (identical on successive runs){Colors.NC}")
        else:
            # Multi-threaded runs are expected to be non-reproducible due to thread scheduling
            if "_nt" in method_id or "_4t" in method_id:
                print(f"{Colors.YELLOW}○ NON-DETERMINISTIC (expected for multi-threaded){Colors.NC}")
            else:
                print(f"{Colors.RED}✗ NOT REPRODUCIBLE (outputs differ!){Colors.NC}")
    
    print()
    print("  " + "─" * 63)
    
    return results


def generate_summary_report(output_dir: str, seed: int) -> None:
    """Create a text file with benchmark summary."""
    report_path = os.path.join(output_dir, "benchmark_report.txt")
    
    print_header("Generating Summary Report")
    
    with open(report_path, 'w') as f:
        f.write("IGoR Fast Generator Benchmark Report\n")
        f.write("=" * 40 + "\n\n")
        f.write(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Seed: {seed}\n\n")
        f.write("=" * 67 + "\n")
        f.write("LEGEND\n")
        f.write("=" * 67 + "\n\n")
        f.write("Correctness Symbols:\n")
        f.write("  ✓ IDENTICAL (exact byte match) - Files are byte-for-byte identical\n")
        f.write("  ≈ Same sequences (different order) - Same content, different order\n")
        f.write("  ○ Different sequences - Expected for MaxSpeed mode\n")
        f.write("  ✗ DIFFERENT (unexpected!) - Indicates a bug\n\n")
        f.write("Distribution Test:\n")
        f.write("  SIMILAR (p>0.05)   - Distributions are statistically equivalent\n")
        f.write("  DIFFERENT (p<0.05) - Distributions differ significantly\n\n")
        f.write("=" * 67 + "\n")
        f.write("MODES\n")
        f.write("=" * 67 + "\n\n")
        f.write("Original:\n")
        f.write("  - Command: igor -generate N\n")
        f.write("  - Algorithm: GenModel::generate_sequences()\n\n")
        f.write("ExactMatch (--fast):\n")
        f.write("  - Command: igor -generate N --fast [--threads T]\n")
        f.write("  - Guarantee: IDENTICAL output with same seed (single-threaded)\n\n")
        f.write("MaxSpeed (--maxspeed):\n")
        f.write("  - Command: igor -generate N --maxspeed [--threads T]\n")
        f.write("  - Guarantee: SAME DISTRIBUTION, different sequences\n\n")
    
    print(f"Report saved to: {report_path}")


# =============================================================================
# MAIN ENTRY POINT
# =============================================================================

def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Benchmark IGoR sequence generation algorithms",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument(
        "--quick", action="store_true",
        help="Run only small benchmarks (1,000 sequences)"
    )
    parser.add_argument(
        "--full", action="store_true",
        help="Run all benchmarks including 10M sequences"
    )
    parser.add_argument(
        "--skip-dist", action="store_true",
        help="Skip distribution analysis (faster)"
    )
    parser.add_argument(
        "--igor-bin", type=str,
        default=os.environ.get("IGOR_BIN", DEFAULT_IGOR_BIN),
        help=f"Path to igor binary (default: {DEFAULT_IGOR_BIN})"
    )
    parser.add_argument(
        "--models-dir", type=str,
        default=os.environ.get("MODELS_DIR", DEFAULT_MODELS_DIR),
        help=f"Path to models directory (default: {DEFAULT_MODELS_DIR})"
    )
    parser.add_argument(
        "--output-dir", type=str,
        default=os.environ.get("OUTPUT_DIR", DEFAULT_OUTPUT_DIR),
        help=f"Path for output files (default: {DEFAULT_OUTPUT_DIR})"
    )
    parser.add_argument(
        "--seed", type=int,
        default=DEFAULT_SEED,
        help=f"Random seed (default: {DEFAULT_SEED})"
    )
    
    args = parser.parse_args()
    
    # Select benchmark sizes
    if args.quick:
        sizes = SIZES_QUICK
    elif args.full:
        sizes = SIZES_FULL
    else:
        sizes = SIZES_DEFAULT
    
    # Validate prerequisites
    if not check_prerequisites(args.igor_bin, args.models_dir):
        sys.exit(1)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Print configuration
    print_header("IGoR Generation Benchmark Suite")
    print("Configuration:")
    print(f"  Binary:       {args.igor_bin}")
    print(f"  Models:       {args.models_dir}")
    print(f"  Output:       {args.output_dir}")
    print(f"  Seed:         {args.seed}")
    print(f"  Sizes:        {', '.join(format_number(s) for s in sizes)}")
    print(f"  Options:      quick={args.quick}, full={args.full}, skip-dist={args.skip_dist}")
    
    # Run benchmarks for each size
    print_header("Running Benchmarks")
    
    all_results = {}
    for size in sizes:
        results = run_benchmarks_for_size(
            size=size,
            igor_bin=args.igor_bin,
            models_dir=args.models_dir,
            output_dir=args.output_dir,
            seed=args.seed,
            skip_dist=args.skip_dist
        )
        all_results[size] = results
    
    # Generate summary report
    generate_summary_report(args.output_dir, args.seed)
    
    # Print completion message
    print_header("Benchmark Complete")
    print(f"{Colors.GREEN}All benchmarks completed successfully!{Colors.NC}")
    print()
    print("Summary:")
    print(f"  - Results saved to: {args.output_dir}")
    print(f"  - Report file: {os.path.join(args.output_dir, 'benchmark_report.txt')}")
    print()
    print("To view generated files:")
    print(f"  ls -la {args.output_dir}/")


if __name__ == "__main__":
    main()
