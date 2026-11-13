# scripts/tests

Regression tests that drive IGoR through its command line, and the tools around them.

## Regression tracks

Each track runs `igor` on the inputs in `data/input` and compares what it writes with the
goldens in `data/reference`, through [assert_regression.sh](assert_regression.sh).

- [test_run.sh](test_run.sh) runs the tracks and writes one report per track under
  `regression_reports/`: `pixi run test_regression`, `--list` to name the tracks, `-k` to keep
  their outputs.
- The tracks: [test_cli.sh](test_cli.sh), [test_align.sh](test_align.sh),
  [test_inference.sh](test_inference.sh), [test_no_d_align.sh](test_no_d_align.sh),
  [test_generate.sh](test_generate.sh). Each also runs alone, as `pixi run test_<name>` where
  pixi.toml has the task.
- [config.sh](config.sh) names the binary they run, `build/bin/igor` unless `IGORBIN` names
  another, and the data directories.
- [test_run_conda.sh](test_run_conda.sh) runs them with the igor of a built conda package, installed
  in an environment of its own: `pixi build --output-dir=bundle`, then
  `pixi run test_regression_conda`. It tests the binary that ships, built with the package's
  compiler and flags, which the development build does not share.

## In CI

[ci.yml](../../.github/workflows/ci.yml) runs:

- on Linux, every track on the conda package the job has just built (`test_regression_conda`);
- on macOS and Windows, the `align` track on the development build;
- nightly, every track on a Debug build, where assertions and the tensor bounds checks are on.

## Other tools

- [test_package.sh](test_package.sh): the installed CMake package, from a consumer project
  (`pixi run test_package`).
- [benchmark_run.sh](benchmark_run.sh) and [BENCHMARK.md](BENCHMARK.md): scalability benchmarks
  (`pixi run benchmark`).
- [profile_run.sh](profile_run.sh): perf profiles of the aligner, on Linux (`pixi run profile`).
- [coverage_report.py](coverage_report.py): per-function coverage of the event code under
  migration (`pixi run coverage`), from a GCC coverage build.
- [test_igor_convert.sh](test_igor_convert.sh): the conversions of `igor-convert`; no pixi task
  runs it.
- `models/`: two monoscenario model parameter files that no script reads.
