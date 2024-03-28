# Rscript your_script.R /path/to/input /path/to/output 2000
# Load required libraries
library(Seurat)
library(Matrix)
library(rhdf5)

# Command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Three arguments must be supplied: input_dir, output_dir, and num_hvg.", call. = FALSE)
}

print("Starting sctransform R script.")

input_dir <- args[1]
output_dir <- args[2]
num_hvg <- as.integer(args[3])

# Read data
seurat_data <- Read10X(data.dir = input_dir)

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = seurat_data)

# SCTransform normalization
seurat_obj <- SCTransform(seurat_obj, variable.features.n = num_hvg, verbose = FALSE)
print("Successfully completed SCTransform normalization.")

# Get the count matrix, feature, and cell barcode data
rna_counts <- GetAssayData(seurat_obj, assay = "SCT", slot = "scale.data")
barcodes <- colnames(rna_counts)
features <- rownames(rna_counts)

print("Trying to get rna_counts.")
# Transpose rna_counts
rna_counts <- t(rna_counts)
print("Successfully got RNA counts trying to write file.")

# Ensure the output directory exists, and if not, create it
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Function to write the file with retry mechanism
write_file_with_retry <- function(file_path, max_attempts, wait_time) {
  attempt <- 1
  while (attempt <= max_attempts) {
    tryCatch({
      # Attempt to write the file
      # Replace this with your actual file writing code
      # For example: h5write(..., file = matrix_file_path)
      message("Attempting to write the file.")
      h5createFile(file_path)
      h5write(rna_counts, file_path, "count_matrix")
      # If file is written successfully, break out of the loop
      message("File written successfully")
      return(TRUE)
    }, error = function(e) {
      # If an error occurs, print the error message and retry after waiting
      message(paste("Attempt", attempt, "failed:", e$message))
      Sys.sleep(wait_time)
    })
    attempt <- attempt + 1
  }
  message("Maximum attempts reached. File write failed.")
  return(FALSE)
}

# Define the file path with .gz extension for gzipped output
matrix_file_path <- file.path(output_dir, "matrix.h5")

# Try writing the file with up to 10 attempts, waiting 60 seconds between attempts
# Total wait time will be at most 10 minutes
result <- write_file_with_retry(matrix_file_path, max_attempts = 1000, wait_time = 10)

# Check if writing was successful
if (result) {
  # Check if the file actually exists
  if (file.exists(matrix_file_path)) {
    print("File exists: Writing was successful.")
    # Proceed with the rest of your script
    print("Continuing with the rest of the script.")
  } else {
    print("File does not exist: Writing was not successful despite the result.")
  }
} else {
  # Print a failure message and exit the script with a non-zero status
  print("Failed to write the file after multiple attempts.")
  quit(status = 1, save = "no")
}



# Write features and barcodes to TSV files and then gzip them
features_file_path <- file.path(output_dir, "features.tsv")
write.table(features, features_file_path, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
system(paste("gzip -f", features_file_path))

barcodes_file_path <- file.path(output_dir, "barcodes.tsv")
write.table(barcodes, barcodes_file_path, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
system(paste("gzip -f", barcodes_file_path))

print("Completed sctransform R script.")