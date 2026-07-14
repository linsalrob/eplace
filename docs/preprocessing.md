# Amplicon preprocessing

`eplace-preprocess` is a standalone pre-annotation command. It converts multiplexed paired-end amplicon reads into an ePLACE-ready ASV FASTA and count tables. It does not run taxonomic annotation.

## Workflow

1. Read synchronized R1, R2, I1 and I2 FASTQ files.
2. Assign each read pair to a sample using the dual-index sample sheet.
3. Trim the forward and reverse locus primers with Cutadapt.
4. Filter reads, learn errors, infer ASVs and merge paired reads with DADA2.
5. Remove chimeras and export ASV sequences and sample counts.

## Requirements

The following programs must be available on `PATH`:

- Cutadapt
- R and `Rscript`
- the R package `dada2`

For example, these may be installed in a Pixi or Conda environment.

## Input files

The command currently expects a standard paired-end four-read set:

- one multiplexed R1 FASTQ or FASTQ.GZ file;
- one multiplexed R2 FASTQ or FASTQ.GZ file;
- one I1 index FASTQ or FASTQ.GZ file;
- one I2 index FASTQ or FASTQ.GZ file;
- a tab-separated sample sheet.

All four FASTQ files must contain the same reads in the same order.

### Sample sheet

The sample sheet must contain exactly the following required columns:

```text
sample_id	index1	index2
Sample01	CTAGCGAA	CTAGTATG
Sample02	GCTCATGA	CGTCTAAT
```

Index matching is exact by default. Use `--index-mismatches 1` to permit one mismatch independently in I1 and I2. Reads with no match or an equally good tie are not assigned.

## Example

```bash
eplace-preprocess \
  --reads-r1 run_R1.fastq.gz \
  --reads-r2 run_R2.fastq.gz \
  --index-i1 run_I1.fastq.gz \
  --index-i2 run_I2.fastq.gz \
  --sample-sheet indexes.tsv \
  --forward-primer AGAGTTTGATCMTGGCTCAG \
  --reverse-primer GWATTACCGCGGCKGCTG \
  --output-dir preprocessing_results \
  --threads 16 \
  --trunc-len-f 255 \
  --trunc-len-r 250 \
  --max-ee-f 2 \
  --max-ee-r 6 \
  --pool pseudo
```

Set both truncation lengths to `0` to disable fixed-length truncation and retain only DADA2 expected-error, `truncQ` and minimum-length filtering.

## Main outputs

The final files are written under `03_dada2/`:

- `asv_sequences.fasta`: ASV IDs and sequences for direct use with ePLACE;
- `asv_table.tsv`: one row per ASV, containing `asv_id`, `sequence`, and counts for every sample;
- `asv_counts_samples_by_rows.tsv`: samples as rows and ASVs as columns;
- `read_tracking.tsv`: input, filtered, denoised, merged and non-chimeric read counts;
- `filtering_summary.tsv`: DADA2 filtering statistics;
- `sequence_table_asv.rds`: the final DADA2 sequence table;
- `sequence_table_prechimera.rds`: the merged table before chimera removal.

Additional outputs include demultiplexed FASTQ files, primer-trimmed FASTQ files, Cutadapt logs and DADA2 error models.

## ePLACE handoff

The FASTA can be passed directly to ePLACE:

```bash
eplace grouped \
  preprocessing_results/03_dada2/asv_sequences.fasta \
  eplace_results \
  --output-classification eplace_results/classifications.tsv
```

The `asv_id` values in the classification output can then be joined back to `asv_table.tsv`.

## Scope

This first implementation demultiplexes reads from separate I1/I2 FASTQ files. It does not yet handle BCL files, single-end reads, inline barcodes embedded in R1/R2, or already-demultiplexed per-sample FASTQ directories.
