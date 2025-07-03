#!/usr/bin/env Rscript

# Load required packages efficiently
required_packages <- c('optparse', 'GenomicRanges', 'ranger', 'matrixStats', 'data.table', 'glmnet')
for (package in required_packages) {
  if (!require(package, character.only=TRUE, quietly=TRUE)) {
    install.packages(package, repos = "http://cran.us.r-project.org")
    library(package, character.only=TRUE)
  }
}

# Function to process methylation data in chunks
process_methylation_chunk <- function(file_path, chunk_size = 100000) {
  # Initialize empty list to store processed chunks
  processed_chunks <- list()
  chunk_count <- 0
  
  # Get total number of lines in file using R's native capabilities
  con <- file(file_path, "r")
  header <- readLines(con, n = 1)  # Read header
  total_lines <- 0
  while(length(line <- readLines(con, n = 1)) > 0) {
    total_lines <- total_lines + 1
  }
  close(con)
  
  message(sprintf("Total lines in file (excluding header): %d", total_lines))
  
  # Read file in chunks using data.table for efficiency
  while(TRUE) {
    chunk_count <- chunk_count + 1
    skip_rows <- 1 + ((chunk_count - 1) * chunk_size)  # +1 for header
    
    # Check if we've reached the end of the file
    if (skip_rows > total_lines) {
      break
    }
    
    # Calculate rows to read for this chunk
    rows_to_read <- min(chunk_size, total_lines - skip_rows + 1)
    
    message(sprintf("Processing chunk %d: reading %d rows starting at line %d", 
                   chunk_count, rows_to_read, skip_rows))
    
    # Read chunk
    chunk <- tryCatch({
      fread(file_path, 
            select = c(1:3, 6, 10:11),  # Only read needed columns
            col.names = c("chr", "start", "end", "strand", "cov", "methylation_percent"),
            skip = skip_rows,
            nrows = rows_to_read,
            data.table = FALSE)
    }, error = function(e) {
      message(sprintf("Error reading chunk %d: %s", chunk_count, e$message))
      return(NULL)
    })
    
    # Break if no data was read
    if (is.null(chunk) || nrow(chunk) == 0) {
      break
    }
    
    # Process chunk
    chunk <- subset(chunk, !(chr %in% c("chrX", "chrY")))
    chunk <- chunk[complete.cases(chunk), ]
    chunk$cov <- as.numeric(chunk$cov)
    chunk$methylation_percent <- as.numeric(chunk$methylation_percent)
    chunk$methylation <- chunk$methylation_percent/100
    
    # Convert to GRanges and store
    chunk_gr <- makeGRangesFromDataFrame(chunk, keep.extra.columns = TRUE)
    processed_chunks[[chunk_count]] <- chunk_gr
    
    # Clear memory
    rm(chunk, chunk_gr)
    gc()
  }
  
  message(sprintf("Processed %d chunks", chunk_count - 1))
  
  # Combine all chunks
  if (length(processed_chunks) == 0) {
    stop("No valid data was processed from the input file")
  }
  
  combined_gr <- do.call(c, processed_chunks)
  rm(processed_chunks)
  gc()
  
  return(combined_gr)
}

# Parse arguments
option_list = list(
  make_option(c("-s", "--sample"), type="character", default=NULL, 
              help="sample", metavar="character"),
  make_option(c("-o", "--out_dir"), type="character", default=NULL, 
              help="output directory", metavar="character"),
  make_option(c("-i", "--in_file"), type="character", default=NULL,
              help="path to modified base called file", metavar="character"),
  make_option(c("-p", "--probes_file"), type="character", default="top_probes_hm450.Rdata"),
  make_option(c("-a", "--array_file"), type="character", default="top_probes_hm450.Rdata"),
  make_option(c("-d", "--training_data"), type="character", default="capper_betas.RData",
              help="capper dataset", metavar="character"),
  make_option(c("-t", "--threads"), type="numeric", default=16,
              help="number of threads", metavar="character"),
  make_option(c("-c", "--chunk_size"), type="numeric", default=100000,
              help="chunk size for processing methylation data", metavar="numeric")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Create output directory
dir.create(file.path(opt$out_dir), showWarnings = FALSE)

# Process methylation data in chunks
message("Processing methylation data in chunks...")
cur38 <- process_methylation_chunk(opt$in_file, chunk_size = opt$chunk_size)
genome(cur38) <- "hg38"

# Load array data efficiently
message("Loading array data...")
load(opt$array_file)
load(opt$probes_file)

# Find overlaps efficiently
message("Finding overlaps...")
meth_38_450k <- subsetByOverlaps(hm450, cur38)
index <- findOverlaps(hm450, cur38)
hm450.matched <- hm450[queryHits(index)]
mcols(hm450.matched) <- cbind.data.frame(mcols(hm450.matched), 
                                        mcols(cur38[subjectHits(index)]))
methylation_sample <- hm450.matched[, c("cov", "methylation")]

# Clear unnecessary objects
rm(cur38, hm450, meth_38_450k, index, hm450.matched)
gc()

# Save intermediate results
fileout <- paste0(opt$out_dir, "/", opt$sample, "_methylation_hg38_HM450.RDS")
saveRDS(methylation_sample, fileout)

# Classification preparation
message("Starting methylation classification...")
case <- methylation_sample
case$cpg <- names(case)
probes <- intersect(unique(names(case)), top_probes)
message(paste(length(probes), " overlapping CpG sites between sample and reference set.", sep=""))

# Save probe information - fix the write.csv warning
write.table(probes, 
          file = paste0(opt$out_dir, "/", opt$sample, "_probes_for_training.csv"),
          row.names = FALSE, 
          col.names = FALSE,
          quote = FALSE)

# Load and process training data efficiently
message("Loading training data...")

# Function to process training data in chunks
process_training_chunk <- function(file_path, probes, chunk_size = 1000) {
  # Load the data once to get structure
  load(file_path)
  all_cols <- colnames(betas)
  total_rows <- nrow(betas)
  
  # Find which columns we need - pre-allocate memory
  needed_cols <- c(probes, "Dx")
  col_indices <- which(all_cols %in% needed_cols)
  rm(all_cols)  # Clear unused data
  gc()
  
  # Calculate optimal chunk size
  mem_info <- gc()
  available_mem <- mem_info[2, 2]  # Get available memory in MB
  
  # Validate available memory
  if (is.na(available_mem) || available_mem <= 0) {
    message("Warning: Could not determine available memory, using default chunk size")
    chunk_size <- 1000L
  } else {
    # Use 50% of available memory, with a minimum of 100MB
    available_mem <- max(available_mem * 0.5, 100)
    estimated_row_size <- 8 * length(col_indices)  # 8 bytes per numeric value
    
    # Calculate chunk size with validation
    optimal_chunk_size <- floor(available_mem * 1024 * 1024 / estimated_row_size)  # Convert MB to bytes
    if (is.na(optimal_chunk_size) || optimal_chunk_size <= 0) {
      message("Warning: Invalid optimal chunk size calculated, using default")
      chunk_size <- 1000L
    } else {
      chunk_size <- as.integer(min(chunk_size, optimal_chunk_size))
    }
  }
  
  # Ensure chunk_size is valid
  chunk_size <- max(1L, min(chunk_size, total_rows))
  chunks <- as.integer(ceiling(total_rows / chunk_size))
  
  # Validate chunks
  if (is.na(chunks) || chunks <= 0) {
    stop("Invalid number of chunks calculated. Please check input data size.")
  }
  
  # Initialize storage for variable probes
  sds <- numeric(length(probes))
  names(sds) <- probes
  
  # Log the chunking parameters
  message(sprintf("Total rows: %i", total_rows))
  message(sprintf("Number of columns: %i", length(col_indices)))
  message(sprintf("Processing training data in %i chunks (chunk size: %i)...", 
                 chunks, chunk_size))
  
  # Process chunks from the already loaded data
  for (i in 1:chunks) {
    start_row <- as.integer(((i-1) * chunk_size) + 1)
    end_row <- as.integer(min(i * chunk_size, total_rows))
    
    # Validate row indices
    if (is.na(start_row) || is.na(end_row) || start_row > end_row) {
      stop(sprintf("Invalid row indices calculated for chunk %i: start=%i, end=%i", 
                  i, start_row, end_row))
    }
    
    message(sprintf("Processing chunk %i/%i (rows %i-%i)", 
                   i, chunks, start_row, end_row))
    
    # Process chunk from the loaded data
    chunk <- betas[start_row:end_row, col_indices, drop = FALSE]
    
    # Calculate SDs more efficiently using matrixStats
    chunk_matrix <- as.matrix(chunk[, -ncol(chunk), drop = FALSE])
    chunk_sds <- matrixStats::colSds(chunk_matrix, na.rm = TRUE)
    sds[names(chunk_sds)] <- pmax(sds[names(chunk_sds)], chunk_sds, na.rm = TRUE)
    
    # Clear chunk data immediately
    rm(chunk, chunk_matrix, chunk_sds)
    gc()
  }
  
  # Select most variable probes more efficiently
  max_CpG <- 10000L  # Use integer literal
  # Use partial sort for better performance
  maxSDs <- head(order(sds, decreasing = TRUE), n = min(length(sds), max_CpG))
  selected_probes <- names(sds)[maxSDs]
  
  # Create final dataset with only selected probes
  ts <- betas[, c(selected_probes, "Dx")]
  
  # Clear all temporary objects
  rm(betas, sds, maxSDs)
  gc()
  
  return(ts)
}

# Process training data in chunks
ts <- process_training_chunk(opt$training_data, probes)
cols <- colnames(ts[, -ncol(ts)])

# Calculate class weights
ts$Dx <- as.factor(ts$Dx)
Dx_fractions <- min(summary(ts$Dx, maxsum = 10000)) / summary(ts$Dx, maxsum = 10000)

# Train random forest model
message("Training random forest model...")
rf <- ranger(dependent.variable.name = "Dx", 
             data = ts[, c(cols, "Dx")], 
             num.trees = 1000, 
             probability = TRUE, 
             sample.fraction = Dx_fractions,
             write.forest = TRUE,  # We need this for predictions
             #save.memory = TRUE,   # Use memory-efficient storage
             verbose = TRUE,
             num.threads = opt$threads)

# Calculate predictions and scores
message("Making predictions...")
probs <- predict(rf, ts, predict.all = FALSE)$predictions
scores <- unlist(lapply(1:nrow(probs), function(i) {
  probs[i, which.max(probs[i, ])]
}))
pred <- colnames(probs)[apply(probs, 1, which.max)]
classes <- levels(ts$Dx)

# Train calibration models
message("Training calibration models...")
glm_models <- lapply(classes, function(type) {
  scores <- probs[, type]
  scale_this <- data.frame(
    class = ifelse(pred == type & pred == ts$Dx, 1, 0),
    score = scores
  )
  glm(class ~ score, data = scale_this, family = binomial)
})

# Clear training data and intermediate objects
rm(ts, probs, scores, pred)
gc()

# Process case data
message("Processing case data...")
case <- as.data.frame(unique(case))
case <- case[match(cols, case$cpg), ]
case <- case[, c("methylation", "cpg")]
case$methylation <- (case$methylation >= 0.5) + 0
df <- subset(case, select = methylation)
df <- t(df)

# Make predictions
message("Making final predictions...")
x_probs <- predict(rf, rbind(df[, cols], df[, cols]), predict.all = FALSE)$predictions[1, ]
x_score <- x_probs[which.max(x_probs)]
x_pred <- names(x_score)

# Store number of variables before cleaning up
num_variables <- rf$num.independent.variables

# Calculate votes and calibration
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

# Clean up model after predictions
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

write.table(report, 
            file = paste0(opt$out_dir, "/", opt$sample, "_calibrated_classification.tsv"),
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

write.table(votes,
            file = paste0(opt$out_dir, "/", opt$sample, "_votes.tsv"),
            quote = FALSE) 