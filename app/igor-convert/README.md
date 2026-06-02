# igor-convert

Command-line tool for converting between IGoR data formats.

## Overview

`igor-convert` enables seamless conversion between:
- **Parquet** (high-performance binary format)
- **AIRR Rearrangement** (TSV/CSV, one row per sequence)
- **AIRR Alignment** (TSV/CSV, one row per alignment)

## Usage

```bash
igor-convert [options] <input> <output>
```

### Options

- `--compression <type>` — Parquet compression (NONE, SNAPPY, GZIP, ZSTD, LZ4)  
  Default: SNAPPY

- `--delimiter <type>` — AIRR delimiter (TAB, COMMA)  
  Default: TAB for .tsv, COMMA for .csv

- `--help` — Show help message

## Examples

### Parquet ↔ AIRR Rearrangement

```bash
# Parquet to AIRR Rearrangement
igor-convert sequences.parquet output.tsv

# AIRR Rearrangement to Parquet with ZSTD compression
igor-convert --compression ZSTD input.tsv output.parquet
```

### Parquet ↔ AIRR Alignment

```bash
# Parquet to AIRR Alignment
igor-convert sequences.parquet alignments.tsv

# AIRR Alignment to Parquet
igor-convert alignments.tsv sequences.parquet
```

### AIRR Rearrangement ↔ AIRR Alignment

```bash
# Rearrangement to Alignment (one row per alignment)
igor-convert rearrangement.tsv alignment.tsv

# Alignment to Rearrangement (one row per sequence)
igor-convert alignment.tsv rearrangement.tsv
```

## Supported Conversions

| From | To | Notes |
|------|-----|-------|
| Parquet | AIRR Rearrangement | Full round-trip |
| Parquet | AIRR Alignment | Full round-trip |
| AIRR Rearrangement | AIRR Alignment | Schema transformation |
| AIRR Alignment | AIRR Rearrangement | Schema transformation |

## Format Auto-Detection

The tool automatically detects input and output formats based on file extensions and content:

- `.parquet` → Parquet format
- `.tsv` / `.csv` → Auto-detect AIRR schema from header:
  - Contains `v_call`/`d_call`/`j_call` → AIRR Rearrangement
  - Contains `segment` → AIRR Alignment

## Building

The tool is built automatically with IGoR:

```bash
pixi run cmake --build build
```

The executable will be in `build/bin/igor-convert`.

## Use Cases

### 1. Data Migration
Convert legacy IGoR data to modern Parquet format for better performance:
```bash
igor-convert legacy_sequences.tsv sequences.parquet --compression ZSTD
```

### 2. Interoperability
Share IGoR results with AIRR-compliant tools (immunarch, scirpy):
```bash
igor-convert sequences.parquet airr_output.tsv
```

### 3. Format Selection
- Use **Parquet** for large-scale internal analysis (fast I/O)
- Use **AIRR TSV** for data sharing and publication

### 4. Schema Transformation
Convert between AIRR Rearrangement (one row per sequence) and AIRR Alignment (one row per alignment):
```bash
igor-convert rearrangement.tsv alignment.tsv
```

## Known Limitations

See `STREAMING_TASK_PLAN.md` for details on:
- CIGAR string simplification in AIRR Alignment format
- Missing `sequence` field in AIRR Alignment schema
- Optional AIRR fields not implemented

## Dependencies

- IGoR Streaming module
- IGoR Core module
- Apache Arrow/Parquet (via Streaming)
- Sparrow (header-only)
