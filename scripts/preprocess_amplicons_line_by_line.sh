#!/usr/bin/env bash

# ==============================================================================
# Simple amplicon preprocessing script
# ==============================================================================
# Purpose:
#   1. Point this script at already-demultiplexed FASTQ files.
#   2. Remove primers with Cutadapt.
#   3. Run DADA2.
#   4. Export ASV IDs, sequences, and counts per sample.
#
# This script is deliberately explicit rather than wrapped as a polished CLI.
# Edit the paths and parameters below, then run sections line by line while
# troubleshooting.
#
# Required software:
#   - cutadapt
#   - Rscript
#   - R packages: dada2, ShortRead, Biostrings
#
# Supported input layouts:
#   INPUT_TYPE="paired"  : one R1 and one R2 FASTQ per sample
#   INPUT_TYPE="single"  : one R1 FASTQ per sample
#
# This script assumes the sequencing facility has already demultiplexed samples.
# If you receive one multiplexed FASTQ containing inline sample barcodes, perform
# that facility/library-specific demultiplexing first, then point this script at
# the resulting per-sample FASTQ files.
# ============================================================================== 

set -euo pipefail

# ==============================================================================
# 1. USER SETTINGS
# ==============================================================================

# Choose "paired" or "single".
INPUT_TYPE="paired"

# Folder containing the input FASTQ files.
RAW_DIR="/path/to/raw_fastq"

# Main output folder.
OUT_DIR="/path/to/amplicon_preprocessing"

# File naming patterns.
# Examples matched by these defaults:
#   Sample01_R1.fastq.gz
#   Sample01_R2.fastq.gz
R1_SUFFIX="_R1.fastq.gz"
R2_SUFFIX="_R2.fastq.gz"

# Primers. Degenerate IUPAC bases are accepted by Cutadapt.
FORWARD_PRIMER="AGAGTTTGATCMTGGCTCAG"
REVERSE_PRIMER="GWATTACCGCGGCKGCTG"

# Cutadapt settings.
CUTADAPT_THREADS=8
CUTADAPT_ERROR_RATE=0.15
CUTADAPT_MIN_LENGTH=50
CUTADAPT_ROUNDS=2
DISCARD_UNTRIMMED=true

# DADA2 settings.
DADA2_THREADS=8
TRUNC_LEN_F=255
TRUNC_LEN_R=250
MAX_EE_F=2
MAX_EE_R=6
TRUNC_Q=6
MIN_LEN=50
POOL_MODE="pseudo"          # "pseudo", "true", or "false"
MIN_OVERLAP=10
MAX_MISMATCH=1
CHIMERA_METHOD="consensus"
MIN_FOLD_PARENT=4

# ==============================================================================
# 2. OUTPUT FOLDERS
# ==============================================================================

TRIM_DIR="${OUT_DIR}/01_cutadapt"
DADA2_DIR="${OUT_DIR}/02_dada2"
LOG_DIR="${OUT_DIR}/logs"

mkdir -p "${TRIM_DIR}" "${DADA2_DIR}" "${LOG_DIR}"

# ==============================================================================
# 3. CHECK SOFTWARE
# ==============================================================================

command -v cutadapt >/dev/null 2>&1 || {
    echo "ERROR: cutadapt was not found in PATH" >&2
    exit 1
}

command -v Rscript >/dev/null 2>&1 || {
    echo "ERROR: Rscript was not found in PATH" >&2
    exit 1
}

# ==============================================================================
# 4. INSPECT INPUT FILES
# ==============================================================================

find "${RAW_DIR}" -maxdepth 1 -type f | sort

# Optional integrity check for compressed input files.
find "${RAW_DIR}" -maxdepth 1 -type f -name "*.fastq.gz" -print0 |
while IFS= read -r -d '' file; do
    echo "Checking ${file}"
    gzip -t "${file}"
done

# ==============================================================================
# 5. CUTADAPT PRIMER REMOVAL
# ==============================================================================

if [[ "${INPUT_TYPE}" == "paired" ]]; then

    # --------------------------------------------------------------------------
    # Paired-end input
    # --------------------------------------------------------------------------

    shopt -s nullglob
    r1_files=("${RAW_DIR}"/*"${R1_SUFFIX}")
    shopt -u nullglob

    if [[ ${#r1_files[@]} -eq 0 ]]; then
        echo "ERROR: No R1 files matched ${RAW_DIR}/*${R1_SUFFIX}" >&2
        exit 1
    fi

    for r1 in "${r1_files[@]}"; do

        filename=$(basename "${r1}")
        sample_id=${filename%${R1_SUFFIX}}
        r2="${RAW_DIR}/${sample_id}${R2_SUFFIX}"

        if [[ ! -f "${r2}" ]]; then
            echo "ERROR: Missing R2 file for ${sample_id}: ${r2}" >&2
            exit 1
        fi

        trimmed_r1="${TRIM_DIR}/${sample_id}_R1_trimmed.fastq.gz"
        trimmed_r2="${TRIM_DIR}/${sample_id}_R2_trimmed.fastq.gz"
        cutadapt_log="${LOG_DIR}/${sample_id}_cutadapt.log"

        echo "Trimming paired sample: ${sample_id}"

        cutadapt \
            -g "^${FORWARD_PRIMER}" \
            -G "^${REVERSE_PRIMER}" \
            -e "${CUTADAPT_ERROR_RATE}" \
            -n "${CUTADAPT_ROUNDS}" \
            -m "${CUTADAPT_MIN_LENGTH}" \
            -j "${CUTADAPT_THREADS}" \
            $([[ "${DISCARD_UNTRIMMED}" == "true" ]] && echo "--discard-untrimmed") \
            -o "${trimmed_r1}" \
            -p "${trimmed_r2}" \
            "${r1}" \
            "${r2}" \
            > "${cutadapt_log}" 2>&1

    done

elif [[ "${INPUT_TYPE}" == "single" ]]; then

    # --------------------------------------------------------------------------
    # Single-end input
    # --------------------------------------------------------------------------

    shopt -s nullglob
    r1_files=("${RAW_DIR}"/*"${R1_SUFFIX}")
    shopt -u nullglob

    if [[ ${#r1_files[@]} -eq 0 ]]; then
        echo "ERROR: No single-end files matched ${RAW_DIR}/*${R1_SUFFIX}" >&2
        exit 1
    fi

    for r1 in "${r1_files[@]}"; do

        filename=$(basename "${r1}")
        sample_id=${filename%${R1_SUFFIX}}

        trimmed_r1="${TRIM_DIR}/${sample_id}_R1_trimmed.fastq.gz"
        cutadapt_log="${LOG_DIR}/${sample_id}_cutadapt.log"

        echo "Trimming single-end sample: ${sample_id}"

        # -a removes the reverse-complemented reverse primer if it occurs near
        # the 3' end of a single read spanning the complete short amplicon.
        cutadapt \
            -g "^${FORWARD_PRIMER}" \
            -a "${REVERSE_PRIMER}" \
            -e "${CUTADAPT_ERROR_RATE}" \
            -n "${CUTADAPT_ROUNDS}" \
            -m "${CUTADAPT_MIN_LENGTH}" \
            -j "${CUTADAPT_THREADS}" \
            $([[ "${DISCARD_UNTRIMMED}" == "true" ]] && echo "--discard-untrimmed") \
            -o "${trimmed_r1}" \
            "${r1}" \
            > "${cutadapt_log}" 2>&1

    done

else
    echo "ERROR: INPUT_TYPE must be paired or single" >&2
    exit 1
fi

# Review primer-removal summaries before continuing.
grep -H -E "Total read|Reads with adapters|Pairs written|Reads written" \
    "${LOG_DIR}"/*_cutadapt.log || true

# ==============================================================================
# 6. PASS SETTINGS TO R
# ==============================================================================

export INPUT_TYPE
export TRIM_DIR
export DADA2_DIR
export DADA2_THREADS
export TRUNC_LEN_F
export TRUNC_LEN_R
export MAX_EE_F
export MAX_EE_R
export TRUNC_Q
export MIN_LEN
export POOL_MODE
export MIN_OVERLAP
export MAX_MISMATCH
export CHIMERA_METHOD
export MIN_FOLD_PARENT

# ==============================================================================
# 7. RUN DADA2
# ==============================================================================

Rscript --vanilla - <<'RSCRIPT'

suppressPackageStartupMessages({
    library(dada2)
    library(ShortRead)
    library(Biostrings)
})

input_type <- Sys.getenv("INPUT_TYPE")
trim_dir <- Sys.getenv("TRIM_DIR")
out_dir <- Sys.getenv("DADA2_DIR")

threads <- as.integer(Sys.getenv("DADA2_THREADS"))
trunc_len_f <- as.integer(Sys.getenv("TRUNC_LEN_F"))
trunc_len_r <- as.integer(Sys.getenv("TRUNC_LEN_R"))
max_ee_f <- as.numeric(Sys.getenv("MAX_EE_F"))
max_ee_r <- as.numeric(Sys.getenv("MAX_EE_R"))
trunc_q <- as.integer(Sys.getenv("TRUNC_Q"))
min_len <- as.integer(Sys.getenv("MIN_LEN"))
pool_mode_text <- tolower(Sys.getenv("POOL_MODE"))
min_overlap <- as.integer(Sys.getenv("MIN_OVERLAP"))
max_mismatch <- as.integer(Sys.getenv("MAX_MISMATCH"))
chimera_method <- Sys.getenv("CHIMERA_METHOD")
min_fold_parent <- as.numeric(Sys.getenv("MIN_FOLD_PARENT"))

if (pool_mode_text == "pseudo") {
    pool_mode <- "pseudo"
} else if (pool_mode_text == "true") {
    pool_mode <- TRUE
} else if (pool_mode_text == "false") {
    pool_mode <- FALSE
} else {
    stop("POOL_MODE must be pseudo, true, or false")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
filtered_dir <- file.path(out_dir, "filtered")
dir.create(filtered_dir, recursive = TRUE, showWarnings = FALSE)

get_sample_names <- function(files, suffix) {
    sub(paste0(suffix, "$"), "", basename(files))
}

if (input_type == "paired") {

    fnFs <- sort(list.files(
        trim_dir,
        pattern = "_R1_trimmed\\.fastq\\.gz$",
        full.names = TRUE
    ))

    fnRs <- sort(list.files(
        trim_dir,
        pattern = "_R2_trimmed\\.fastq\\.gz$",
        full.names = TRUE
    ))

    if (length(fnFs) == 0 || length(fnRs) == 0) {
        stop("No paired trimmed FASTQ files were found")
    }

    sample_names_f <- get_sample_names(fnFs, "_R1_trimmed.fastq.gz")
    sample_names_r <- get_sample_names(fnRs, "_R2_trimmed.fastq.gz")

    if (!identical(sample_names_f, sample_names_r)) {
        stop("Forward and reverse sample names do not match")
    }

    sample.names <- sample_names_f

    filtFs <- file.path(filtered_dir, paste0(sample.names, "_F_filt.fastq.gz"))
    filtRs <- file.path(filtered_dir, paste0(sample.names, "_R_filt.fastq.gz"))

    names(fnFs) <- sample.names
    names(fnRs) <- sample.names
    names(filtFs) <- sample.names
    names(filtRs) <- sample.names

    trunc_lengths <- c(trunc_len_f, trunc_len_r)

    filter_out <- filterAndTrim(
        fnFs,
        filtFs,
        fnRs,
        filtRs,
        truncLen = trunc_lengths,
        maxN = 0,
        maxEE = c(max_ee_f, max_ee_r),
        truncQ = trunc_q,
        rm.phix = TRUE,
        compress = TRUE,
        multithread = threads,
        minLen = min_len,
        matchIDs = TRUE
    )

    write.table(
        cbind(sample_id = rownames(filter_out), filter_out),
        file.path(out_dir, "filtering_summary.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    errF <- learnErrors(filtFs, multithread = threads, randomize = TRUE)
    errR <- learnErrors(filtRs, multithread = threads, randomize = TRUE)

    saveRDS(errF, file.path(out_dir, "error_model_forward.rds"))
    saveRDS(errR, file.path(out_dir, "error_model_reverse.rds"))

    derepFs <- derepFastq(filtFs, verbose = TRUE)
    derepRs <- derepFastq(filtRs, verbose = TRUE)
    names(derepFs) <- sample.names
    names(derepRs) <- sample.names

    dadaFs <- dada(
        derepFs,
        err = errF,
        pool = pool_mode,
        multithread = threads
    )

    dadaRs <- dada(
        derepRs,
        err = errR,
        pool = pool_mode,
        multithread = threads
    )

    mergers <- mergePairs(
        dadaFs,
        derepFs,
        dadaRs,
        derepRs,
        minOverlap = min_overlap,
        maxMismatch = max_mismatch,
        verbose = TRUE
    )

    seqtab <- makeSequenceTable(mergers)

    getN <- function(x) sum(getUniques(x))

    tracking <- cbind(
        input = filter_out[, "reads.in"],
        filtered = filter_out[, "reads.out"],
        denoisedF = sapply(dadaFs, getN),
        denoisedR = sapply(dadaRs, getN),
        merged = sapply(mergers, getN)
    )

} else if (input_type == "single") {

    fnFs <- sort(list.files(
        trim_dir,
        pattern = "_R1_trimmed\\.fastq\\.gz$",
        full.names = TRUE
    ))

    if (length(fnFs) == 0) {
        stop("No single-end trimmed FASTQ files were found")
    }

    sample.names <- get_sample_names(fnFs, "_R1_trimmed.fastq.gz")
    filtFs <- file.path(filtered_dir, paste0(sample.names, "_F_filt.fastq.gz"))

    names(fnFs) <- sample.names
    names(filtFs) <- sample.names

    filter_out <- filterAndTrim(
        fnFs,
        filtFs,
        truncLen = trunc_len_f,
        maxN = 0,
        maxEE = max_ee_f,
        truncQ = trunc_q,
        rm.phix = TRUE,
        compress = TRUE,
        multithread = threads,
        minLen = min_len
    )

    write.table(
        cbind(sample_id = rownames(filter_out), filter_out),
        file.path(out_dir, "filtering_summary.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    errF <- learnErrors(filtFs, multithread = threads, randomize = TRUE)
    saveRDS(errF, file.path(out_dir, "error_model_forward.rds"))

    derepFs <- derepFastq(filtFs, verbose = TRUE)
    names(derepFs) <- sample.names

    dadaFs <- dada(
        derepFs,
        err = errF,
        pool = pool_mode,
        multithread = threads
    )

    seqtab <- makeSequenceTable(dadaFs)

    getN <- function(x) sum(getUniques(x))

    tracking <- cbind(
        input = filter_out[, "reads.in"],
        filtered = filter_out[, "reads.out"],
        denoised = sapply(dadaFs, getN)
    )

} else {
    stop("INPUT_TYPE must be paired or single")
}

saveRDS(seqtab, file.path(out_dir, "sequence_table_prechimera.rds"))

seqtab.nochim <- removeBimeraDenovo(
    seqtab,
    method = chimera_method,
    multithread = threads,
    verbose = TRUE,
    minFoldParentOverAbundance = min_fold_parent
)

saveRDS(seqtab.nochim, file.path(out_dir, "sequence_table_asv.rds"))

tracking <- cbind(
    tracking,
    nonchim = rowSums(seqtab.nochim)
)

write.table(
    cbind(sample_id = rownames(tracking), tracking),
    file.path(out_dir, "read_tracking.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# Assign stable ASV labels within this run.
# ------------------------------------------------------------------------------

asv_sequences <- colnames(seqtab.nochim)
asv_ids <- sprintf("ASV%06d", seq_along(asv_sequences))

colnames(seqtab.nochim) <- asv_ids

asv_lookup <- data.frame(
    asv_id = asv_ids,
    sequence = asv_sequences,
    length = nchar(asv_sequences),
    total_reads = colSums(seqtab.nochim),
    prevalence = colSums(seqtab.nochim > 0),
    stringsAsFactors = FALSE
)

write.table(
    asv_lookup,
    file.path(out_dir, "asv_lookup.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# Sample rows, ASV columns.
asv_counts_samples <- data.frame(
    sample_id = rownames(seqtab.nochim),
    seqtab.nochim,
    check.names = FALSE
)

write.table(
    asv_counts_samples,
    file.path(out_dir, "asv_counts_samples_by_rows.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ASV rows with sequence and one count column per sample.
asv_counts_by_asv <- data.frame(
    asv_id = asv_ids,
    sequence = asv_sequences,
    t(seqtab.nochim),
    check.names = FALSE
)

write.table(
    asv_counts_by_asv,
    file.path(out_dir, "asv_id_sequence_counts.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ePLACE-ready FASTA.
asv_dna <- DNAStringSet(asv_sequences)
names(asv_dna) <- asv_ids
writeXStringSet(
    asv_dna,
    file.path(out_dir, "asv_sequences.fasta")
)

cat("\nDADA2 processing complete.\n")
cat("ASVs retained:", ncol(seqtab.nochim), "\n")
cat("Main table:", file.path(out_dir, "asv_id_sequence_counts.tsv"), "\n")
cat("FASTA:", file.path(out_dir, "asv_sequences.fasta"), "\n")

RSCRIPT

# ==============================================================================
# 8. FINAL OUTPUTS
# ==============================================================================

printf '\nMain outputs:\n'
printf '  %s\n' "${DADA2_DIR}/asv_id_sequence_counts.tsv"
printf '  %s\n' "${DADA2_DIR}/asv_sequences.fasta"
printf '  %s\n' "${DADA2_DIR}/asv_counts_samples_by_rows.tsv"
printf '  %s\n' "${DADA2_DIR}/read_tracking.tsv"
printf '  %s\n' "${DADA2_DIR}/filtering_summary.tsv"
