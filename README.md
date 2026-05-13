# ARGs-Gene-Profiling

# 🧬 DIAMOND + CARD Paired-End ARGs Detection Pipeline

A high-confidence antibiotic resistance gene (ARG) detection pipeline using DIAMOND BLASTX against the CARD (Comprehensive Antibiotic Resistance Database), following ARGs-OAP recommendations.

---

## Overview

This pipeline processes **paired-end unmapped reads** (from a prior genome alignment step) and screens them for antibiotic resistance genes using translated nucleotide search against the CARD protein database. Both mates are analyzed independently and their results are combined per sample.

---

## Features

- ✅ Paired-end read support (mate1 + mate2)
- ✅ High-confidence filtering (≥80% identity, ≥25 aa alignment, E-value ≤ 1e-7)
- ✅ Follows ARGs-OAP threshold recommendations
- ✅ Automatic dependency checking and installation via `mamba`
- ✅ Per-sample and aggregate summary reports (CSV)
- ✅ Timestamped log files
- ✅ Modular, readable Bash architecture

---

## Requirements

| Tool | Version | Source |
|------|---------|--------|
| `diamond` | ≥2.0 | [bioconda](https://anaconda.org/bioconda/diamond) |
| `bash` | ≥4.0 | system |
| `awk`, `wc` | any | system (coreutils) |

> If `diamond` is not found, the script will attempt to install it automatically via `mamba`.

---

## Input

| Item | Description |
|------|-------------|
| Paired-end FASTQ files | Unmapped reads from STAR genome alignment |
| File naming pattern | `<sample>_genome_Unmapped.out.mate1` / `...mate2` |
| DIAMOND database | Pre-built `.dmnd` file from the CARD protein database |

Reads are expected to be stored in:
```
/data/work/alignment_results/genome/
```

The DIAMOND database (`.dmnd`) must exist at:
```
/data/work/diamond_db/card_diamond.dmnd
```

---

## Configuration

Edit the constants at the top of the script to match your environment:

```bash
readonly INPUT_DIR="/data/work/alignment_results/genome"
readonly DIAMOND_DB="/data/work/diamond_db/card_diamond"
readonly OUTPUT_DIR="/data/work/diamond_results"
readonly THREADS=8
readonly MIN_IDENTITY=80          # Minimum % identity
readonly MIN_AMINO_ACIDS=25       # Minimum alignment length in amino acids (= 75 nt)
readonly EVALUE="1e-7"            # E-value cutoff
```

---

## Usage

```bash
# Make executable
chmod +x diamond_card_pipeline.sh

# Run the pipeline
bash diamond_card_pipeline.sh
```

No arguments are required. All configuration is handled via the variables at the top of the script.

---

## Output

All results are written to `OUTPUT_DIR` (default: `/data/work/diamond_results/`).

| File | Description |
|------|-------------|
| `<sample>_diamond_card.tsv` | DIAMOND BLASTX output (format 6) per sample |
| `analysis_summary_YYYYMMDD.csv` | Aggregate CSV with hits and quality metrics for all samples |
| `analysis_log_YYYYMMDDHHMMSS.txt` | Full run log with timestamps |

### TSV Columns (DIAMOND format 6)

Standard BLAST tabular format:

```
qseqid  sseqid  pident  length  mismatch  gapopen  qstart  qend  sstart  send  evalue  bitscore
```

### Summary CSV Columns

```
Sample, Total_Hits, High_Confidence_Hits, Avg_Length, Avg_Identity
```

---

## Thresholds & Rationale

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| `MIN_IDENTITY` | 80% | Standard cutoff for confident ARG identification |
| `MIN_AMINO_ACIDS` | 25 aa (75 nt) | ARGs-OAP minimum alignment recommendation |
| `EVALUE` | 1e-7 | Stringent threshold to reduce false positives |
| `--max-target-seqs` | 1 | Retain best hit per query only |

These thresholds follow the [ARGs-OAP v3](https://github.com/biofuture/Ublastx_stageone) recommendations for metagenomic ARG profiling.

---

## Pipeline Steps

```
1. init_analysis()         → Create output directory, start logging
2. check_dependencies()    → Verify DIAMOND is installed
3. validate_diamond_database() → Confirm .dmnd file exists
4. For each mate1 file:
   └─ process_sample()
      ├─ analyze_mate(mate1) → DIAMOND BLASTX
      ├─ analyze_mate(mate2) → DIAMOND BLASTX
      ├─ cat mate1 + mate2  → combined TSV
      └─ validate_results() → print per-sample stats
5. generate_summary()      → print totals + write CSV
```

---

## Example Output (Console)

```
=== 🚀 DIAMOND + CARD PAIRED-END ANALYSIS STARTED ===
📅 Timestamp: 2025-07-01 09:14:32
...
🧬 Processing sample: SRR12345678
   🔬 Running DIAMOND on mate1...
   🔬 Running DIAMOND on mate2...
✅ SUCCESS: SRR12345678
   ⏱️ Time: 43s
   📊 Total hits: 1284
   🎯 High-confidence hits: 1102
   📏 Avg alignment length: 87.3
   🎯 Avg identity: 91.4
```

---

## Database Preparation

If you need to build the DIAMOND database from CARD:

```bash
# Download CARD protein sequences
wget https://card.mcmaster.ca/latest/data -O card_data.tar.bz2
tar -xjf card_data.tar.bz2

# Build DIAMOND database
diamond makedb \
  --in protein_fasta_protein_homolog_model.fasta \
  --db /data/work/diamond_db/card_diamond \
  --threads 8
```

---

## Citation

If you use this pipeline, please cite:

- **DIAMOND**: Buchfink B, Reuter K, Drost HG. *Sensitive protein alignments at tree-of-life scale using DIAMOND.* Nature Methods, 2021. https://doi.org/10.1038/s41592-021-01101-x
- **CARD**: Alcock et al. *CARD 2023: expanded curation, support for machine learning, and resistome prediction at the CARD.* Nucleic Acids Research, 2023. https://doi.org/10.1093/nar/gkac920
- **ARGs-OAP**: Yin et al. *ARGs-OAP v3.0: antibiotic resistance gene database and analysis platform.* iMeta, 2023.

---

## License

MIT License. See `LICENSE` for details.
