#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(dada2))

parse_args <- function(args) {
  result <- list()
  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[[i]])
    if (i == length(args)) stop("Missing value for --", key)
    result[[key]] <- args[[i + 1]]
    i <- i + 2
  }
  result
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("input-dir", "output-dir", "threads", "trunc-len-f", "trunc-len-r",
              "max-ee-f", "max-ee-r", "trunc-q", "min-length", "min-overlap",
              "max-mismatch", "pool")
missing <- required[!required %in% names(opt)]
if (length(missing)) stop("Missing arguments: ", paste(missing, collapse = ", "))

input_dir <- normalizePath(opt[["input-dir"]], mustWork = TRUE)
output_dir <- opt[["output-dir"]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
filtered_dir <- file.path(output_dir, "filtered")
dir.create(filtered_dir, recursive = TRUE, showWarnings = FALSE)

fnFs <- sort(list.files(input_dir, pattern = "\\.trimmed_R1\\.fastq\\.gz$", full.names = TRUE))
fnRs <- sort(list.files(input_dir, pattern = "\\.trimmed_R2\\.fastq\\.gz$", full.names = TRUE))
if (!length(fnFs) || length(fnFs) != length(fnRs)) stop("Could not find matched trimmed R1/R2 files")

sample_names <- sub("\\.trimmed_R1\\.fastq\\.gz$", "", basename(fnFs))
sample_names_r <- sub("\\.trimmed_R2\\.fastq\\.gz$", "", basename(fnRs))
if (!identical(sample_names, sample_names_r)) stop("Trimmed R1/R2 sample names do not match")

filtFs <- file.path(filtered_dir, paste0(sample_names, ".filtered_R1.fastq.gz"))
filtRs <- file.path(filtered_dir, paste0(sample_names, ".filtered_R2.fastq.gz"))
trunc_len <- c(as.integer(opt[["trunc-len-f"]]), as.integer(opt[["trunc-len-r"]]))

filter_args <- list(
  fwd = fnFs,
  filt = filtFs,
  rev = fnRs,
  filt.rev = filtRs,
  maxN = 0,
  maxEE = c(as.numeric(opt[["max-ee-f"]]), as.numeric(opt[["max-ee-r"]])),
  truncQ = as.integer(opt[["trunc-q"]]),
  rm.phix = TRUE,
  compress = TRUE,
  multithread = as.integer(opt[["threads"]]),
  minLen = as.integer(opt[["min-length"]]),
  matchIDs = TRUE,
  verbose = TRUE
)
if (all(trunc_len > 0)) filter_args$truncLen <- trunc_len
out <- do.call(filterAndTrim, filter_args)
write.table(cbind(sample_id = rownames(out), out), file.path(output_dir, "filtering_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

errF <- learnErrors(filtFs, multithread = as.integer(opt[["threads"]]), randomize = TRUE)
errR <- learnErrors(filtRs, multithread = as.integer(opt[["threads"]]), randomize = TRUE)
saveRDS(errF, file.path(output_dir, "error_model_forward.rds"))
saveRDS(errR, file.path(output_dir, "error_model_reverse.rds"))

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
names(derepFs) <- sample_names
names(derepRs) <- sample_names

pool_value <- switch(opt[["pool"]], "false" = FALSE, "true" = TRUE, "pseudo" = "pseudo")
dadaFs <- dada(derepFs, err = errF, multithread = as.integer(opt[["threads"]]), pool = pool_value)
dadaRs <- dada(derepRs, err = errR, multithread = as.integer(opt[["threads"]]), pool = pool_value)

mergers <- mergePairs(
  dadaFs, derepFs, dadaRs, derepRs,
  minOverlap = as.integer(opt[["min-overlap"]]),
  maxMismatch = as.integer(opt[["max-mismatch"]]),
  verbose = TRUE
)
seqtab <- makeSequenceTable(mergers)
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = as.integer(opt[["threads"]]), verbose = TRUE)
saveRDS(seqtab, file.path(output_dir, "sequence_table_prechimera.rds"))
saveRDS(seqtab_nochim, file.path(output_dir, "sequence_table_asv.rds"))

getN <- function(x) sum(getUniques(x))
track <- cbind(
  input = out[, "reads.in"],
  filtered = out[, "reads.out"],
  denoisedF = sapply(dadaFs, getN),
  denoisedR = sapply(dadaRs, getN),
  merged = sapply(mergers, getN),
  nonchim = rowSums(seqtab_nochim)
)
write.table(cbind(sample_id = rownames(track), track), file.path(output_dir, "read_tracking.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

sequences <- colnames(seqtab_nochim)
totals <- colSums(seqtab_nochim)
order_idx <- order(totals, decreasing = TRUE)
sequences <- sequences[order_idx]
seqtab_nochim <- seqtab_nochim[, order_idx, drop = FALSE]
asv_ids <- sprintf("ASV%06d", seq_along(sequences))

fasta <- file(file.path(output_dir, "asv_sequences.fasta"), "wt")
for (i in seq_along(sequences)) {
  writeLines(c(paste0(">", asv_ids[[i]]), sequences[[i]]), fasta)
}
close(fasta)

asv_table <- data.frame(
  asv_id = asv_ids,
  sequence = sequences,
  t(seqtab_nochim),
  check.names = FALSE
)
write.table(asv_table, file.path(output_dir, "asv_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

counts_samples_rows <- data.frame(sample_id = rownames(seqtab_nochim), seqtab_nochim, check.names = FALSE)
colnames(counts_samples_rows)[-1] <- asv_ids
write.table(counts_samples_rows, file.path(output_dir, "asv_counts_samples_by_rows.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message("DADA2 preprocessing complete: ", length(sequences), " ASVs")
