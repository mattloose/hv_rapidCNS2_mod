#!/usr/bin/env Rscript

# =========================
# Robust methylation classifier script (fixed)
# - Avoids saveRDS() on GRanges (prevents containsOutOfMemoryData/elementType slot errors)
# - Upgrades loaded Bioconductor S4 objects (updateObject) defensively
# - Never installs Bioconductor packages from CRAN at runtime (fails fast instead)
# =========================

suppressPackageStartupMessages({
  library(optparse)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(IRanges)
  library(S4Vectors)
  library(ranger)
  library(matrixStats)
  library(data.table)
  library(glmnet)
})

# ---- Helper: fail fast if packages are missing (no runtime installs) ----
# (In conda envs, installing at runtime can create mismatched Bioc stacks.)
assert_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing R packages: ", paste(missing, collapse = ", "),
      "\nInstall them in your conda env (recommended) and re-run.\n",
      call. = FALSE
    )
  }
}

assert_packages(c(
  "optparse", "GenomicRanges", "GenomeInfoDb", "IRanges", "S4Vectors",
  "ranger", "matrixStats", "data.table", "glmnet"
))

# ---- Helper: upgrade old serialized Bioconductor S4 objects (defensive) ----
safe_updateObject <- function(x) {
  tryCatch(updateObject(x), error = function(e) x)
}

# Function to process methylation data in chunks
process_methylation_chunk <- function(file_path, chunk_size = 100000) {
  processed_chunks <- list()
  chunk_count <- 0

  # Count lines (excluding header)
  con <- file(file_path, "r")
  header <- readLines(con, n = 1)
  total_lines <- 0
  while (length(line <- readLines(con, n = 1)) > 0) {
    total_lines <- total_lines + 1
  }
  close(con)

  message(sprintf("Total lines in file (excluding header): %d", total_lines))

  # Read file in chunks
  while (TRUE) {
    chunk_count <- chunk_count + 1
    skip_rows <- 1 + ((chunk_count - 1) * chunk_size) # +1 for header

    if (skip_rows > total_lines) break

    rows_to_read <- min(chunk_size, total_lines - skip_rows + 1)
    message(sprintf(
      "Processing chunk %d: reading %d rows starting at line %d",
      chunk_count, rows_to_read, skip_rows
    ))

    chunk <- tryCatch({
      fread(
        file_path,
        select = c(1:3, 6, 10:11),
        col.names = c("chr", "start", "end", "strand", "cov", "methylation_percent"),
        skip = skip_rows,
        nrows = rows_to_read,
        data.table = FALSE
      )
    }, error = function(e) {
      message(sprintf("Error reading chunk %d: %s", chunk_count, e$message))
      return(NULL)
    })

    if (is.null(chunk) || nrow(chunk) == 0) break

    chunk <- subset(chunk, !(chr %in% c("chrX", "chrY")))
    chunk <- chunk[complete.cases(chunk), , drop = FALSE]
    chunk$cov <- as.numeric(chunk$cov)
    chunk$methylation_percent <- as.numeric(chunk$methylation_percent)
    chunk$methylation <- chunk$methylation_percent / 100

    chunk_gr <- makeGRangesFromDataFrame(chunk, keep.extra.columns = TRUE)
    processed_chunks[[chunk_count]] <- chunk_gr

    rm(chunk, chunk_gr)
    gc()
  }

  message(sprintf("Processed %d chunks", chunk_count - 1))

  if (length(processed_chunks) == 0) {
    stop("No valid data was processed from the input file")
  }

  combined_gr <- do.call(c, processed_chunks)
  rm(processed_chunks)
  gc()

  combined_gr
}

# Parse arguments
option_list <- list(
  make_option(c("-s", "--sample"), type = "character", default = NULL,
              help = "sample", metavar = "character"),
  make_option(c("-o", "--out_dir"), type = "character", default = NULL,
              help = "output directory", metavar = "character"),
  make_option(c("-i", "--in_file"), type = "character", default = NULL,
              help = "path to modified base called file", metavar = "character"),
  make_option(c("-p", "--probes_file"), type = "character",
              default = "top_probes_hm450.Rdata"),
  make_option(c("-a", "--array_file"), type = "character",
              default = "top_probes_hm450.Rdata"),
  make_option(c("-d", "--training_data"), type = "character",
              default = "capper_betas.RData",
              help = "capper dataset", metavar = "character"),
  make_option(c("-t", "--threads"), type = "numeric", default = 16,
              help = "number of threads", metavar = "character"),
  make_option(c("-c", "--chunk_size"), type = "numeric", default = 100000,
              help = "chunk size for processing methylation data", metavar = "numeric")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$sample) || is.null(opt$out_dir) || is.null(opt$in_file)) {
  stop("Missing required arguments: --sample, --out_dir, --in_file", call. = FALSE)
}

# Create output directory
dir.create(file.path(opt$out_dir), showWarnings = FALSE, recursive = TRUE)

# Process methylation data in chunks
message("Processing methylation data in chunks...")
cur38 <- process_methylation_chunk(opt$in_file, chunk_size = opt$chunk_size)
GenomeInfoDb::genome(cur38) <- "hg38"

# Load array and probes data
message("Loading array data...")
load(opt$array_file)   # expected to define hm450 (GRanges)
load(opt$probes_file)  # expected to define top_probes (character vector)

# Defensive upgrades for S4 objects loaded from Rdata
if (exists("hm450")) hm450 <- safe_updateObject(hm450)

# Basic validations
if (!exists("hm450")) stop("Expected object 'hm450' was not found after loading --array_file", call. = FALSE)
if (!is(hm450, "GRanges")) stop("Object 'hm450' is not a GRanges after loading --array_file", call. = FALSE)
if (!exists("top_probes")) stop("Expected object 'top_probes' was not found after loading --probes_file", call. = FALSE)

# Find overlaps
message("Finding overlaps...")

# NOTE: meth_38_450k was unused downstream; keep minimal work
index <- findOverlaps(hm450, cur38)

hm450.matched <- hm450[queryHits(index)]
mcols(hm450.matched) <- cbind.data.frame(
  mcols(hm450.matched),
  mcols(cur38[subjectHits(index)])
)

# ---- FIX: do NOT keep GRanges for downstream / saveRDS ----
# Convert to a simple data.frame with explicit CpG IDs.
methylation_sample <- as.data.frame(
  mcols(hm450.matched)[, c("cov", "methylation"), drop = FALSE]
)
methylation_sample$cpg <- names(hm450.matched)

# Clear unnecessary objects
rm(cur38, hm450, index, hm450.matched)
gc()

# Save intermediate results safely (data.frame, not GRanges)
fileout <- paste0(opt$out_dir, "/", opt$sample, "_methylation_hg38_HM450.RDS")
message("Saving intermediate results: ", fileout)
saveRDS(methylation_sample, fileout)

# Classification preparation
message("Starting methylation classification...")
case <- methylation_sample

probes <- intersect(unique(case$cpg), top_probes)
message(paste(length(probes), " overlapping CpG sites between sample and reference set.", sep = ""))

# Save probe information
write.table(
  probes,
  file = paste0(opt$out_dir, "/", opt$sample, "_probes_for_training.csv"),
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

# Load and process training data efficiently
message("Loading training data...")

process_training_chunk <- function(file_path, probes, chunk_size = 1000) {
  load(file_path) # expected to define betas (data.frame/matrix-like) with "Dx" column

  if (!exists("betas")) stop("Expected object 'betas' was not found after loading training_data", call. = FALSE)
  if (!("Dx" %in% colnames(betas))) stop("Training data 'betas' must contain a 'Dx' column", call. = FALSE)

  all_cols <- colnames(betas)
  total_rows <- nrow(betas)

  needed_cols <- c(probes, "Dx")
  col_indices <- which(all_cols %in% needed_cols)

  rm(all_cols)
  gc()

  # Calculate chunk size based on available memory (best-effort)
  mem_info <- gc()
  available_mem <- mem_info[2, 2] # MB (heuristic)
  if (is.na(available_mem) || available_mem <= 0) {
    message("Warning: Could not determine available memory, using default chunk size")
    chunk_size <- 1000L
  } else {
    available_mem <- max(available_mem * 0.5, 100)
    estimated_row_size <- 8 * length(col_indices)
    optimal_chunk_size <- floor(available_mem * 1024 * 1024 / estimated_row_size)
    if (is.na(optimal_chunk_size) || optimal_chunk_size <= 0) {
      message("Warning: Invalid optimal chunk size calculated, using default")
      chunk_size <- 1000L
    } else {
      chunk_size <- as.integer(min(chunk_size, optimal_chunk_size))
    }
  }

  chunk_size <- max(1L, min(chunk_size, total_rows))
  chunks <- as.integer(ceiling(total_rows / chunk_size))
  if (is.na(chunks) || chunks <= 0) {
    stop("Invalid number of chunks calculated. Please check input data size.")
  }

  sds <- numeric(length(probes))
  names(sds) <- probes

  message(sprintf("Total rows: %i", total_rows))
  message(sprintf("Number of columns: %i", length(col_indices)))
  message(sprintf("Processing training data in %i chunks (chunk size: %i)...", chunks, chunk_size))

  for (i in 1:chunks) {
    start_row <- as.integer(((i - 1) * chunk_size) + 1)
    end_row <- as.integer(min(i * chunk_size, total_rows))
    if (is.na(start_row) || is.na(end_row) || start_row > end_row) {
      stop(sprintf("Invalid row indices for chunk %i: start=%i, end=%i", i, start_row, end_row))
    }

    message(sprintf("Processing chunk %i/%i (rows %i-%i)", i, chunks, start_row, end_row))

    chunk <- betas[start_row:end_row, col_indices, drop = FALSE]

    # If Dx is included, it's presumably last among selected; drop it for SD calc
    chunk_matrix <- as.matrix(chunk[, -ncol(chunk), drop = FALSE])
    chunk_sds <- matrixStats::colSds(chunk_matrix, na.rm = TRUE)

    # Keep max SD observed per probe across chunks
    sds[names(chunk_sds)] <- pmax(sds[names(chunk_sds)], chunk_sds, na.rm = TRUE)

    rm(chunk, chunk_matrix, chunk_sds)
    gc()
  }

  max_CpG <- 10000L
  maxSDs <- head(order(sds, decreasing = TRUE), n = min(length(sds), max_CpG))
  selected_probes <- names(sds)[maxSDs]

  ts <- betas[, c(selected_probes, "Dx"), drop = FALSE]

  rm(betas, sds, maxSDs)
  gc()

  ts
}

ts <- process_training_chunk(opt$training_data, probes)
cols <- colnames(ts[, -ncol(ts), drop = FALSE])

# Calculate class weights
ts$Dx <- as.factor(ts$Dx)
Dx_counts <- summary(ts$Dx, maxsum = 10000)
Dx_fractions <- min(Dx_counts) / Dx_counts

# Train random forest model
message("Training random forest model...")
rf <- ranger(
  dependent.variable.name = "Dx",
  data = ts[, c(cols, "Dx"), drop = FALSE],
  num.trees = 1000,
  probability = TRUE,
  sample.fraction = Dx_fractions,
  write.forest = TRUE,
  verbose = TRUE,
  num.threads = opt$threads
)

# Calculate predictions and scores
message("Making predictions...")
probs <- predict(rf, ts, predict.all = FALSE)$predictions
scores <- vapply(seq_len(nrow(probs)), function(i) probs[i, which.max(probs[i, ])], numeric(1))
pred <- colnames(probs)[apply(probs, 1, which.max)]
classes <- levels(ts$Dx)

# Train calibration models
message("Training calibration models...")
glm_models <- lapply(classes, function(type) {
  sc <- probs[, type]
  scale_this <- data.frame(
    class = ifelse(pred == type & pred == ts$Dx, 1, 0),
    score = sc
  )
  glm(class ~ score, data = scale_this, family = binomial)
})

rm(ts, probs, scores, pred)
gc()

# Process case data
message("Processing case data...")
case <- as.data.frame(unique(case))
case <- case[match(cols, case$cpg), , drop = FALSE]
case <- case[, c("methylation", "cpg"), drop = FALSE]
case$methylation <- (case$methylation >= 0.5) + 0

df <- t(subset(case, select = methylation))

# Make predictions
message("Making final predictions...")
x_probs <- predict(rf, rbind(df[, cols, drop = FALSE], df[, cols, drop = FALSE]), predict.all = FALSE)$predictions[1, ]
x_score <- x_probs[which.max(x_probs)]
x_pred <- names(x_score)

num_variables <- rf$num.independent.variables

votes <- data.frame(x_probs)
colnames(votes) <- "Freq"
votes$Freq <- votes$Freq / sum(votes$Freq) * 100

# Apply calibration
x_scaled <- lapply(seq_along(classes), function(i) {
  type <- classes[i]
  x_scores <- x_probs[type]
  scaled_scores <- glm_models[[i]]
  predict(scaled_scores, newdata = data.frame(score = x_scores), type = "response")
})
x_calibrated_scores <- unlist(x_scaled) / sum(unlist(x_scaled))
x_calibrated_score <- x_calibrated_scores[x_pred]

rm(rf, glm_models)
gc()

votes$cal_Freq <- x_calibrated_scores
votes$cal_Freq <- votes$cal_Freq / sum(votes$cal_Freq) * 100
votes <- votes[order(votes$Freq), , drop = FALSE]

# Save results
message("Saving results...")
report <- paste(
  paste0("Number of features: ", num_variables),
  paste0("Predicted Class: ", x_pred),
  paste0("Initial Score: ", x_score),
  paste0("Calibrated Score: ", x_calibrated_score),
  sep = "\n"
)

write.table(
  report,
  file = paste0(opt$out_dir, "/", opt$sample, "_calibrated_classification.tsv"),
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

write.table(
  votes,
  file = paste0(opt$out_dir, "/", opt$sample, "_votes.tsv"),
  quote = FALSE
)

message("Done.")