#!/usr/bin/env Rscript
library(rmarkdown)

working_dir = getwd()
message(paste("Working directory:", working_dir))

args = commandArgs(trailingOnly = TRUE)

input_params = list(
  path_prefix = args[1],
  vector_type = args[2],
  target_start = as.integer(args[3]),
  target_end = as.integer(args[4]),
  target_start_repcap = as.integer(args[5]),
  target_end_repcap = as.integer(args[6])
)
message("Parameters:")
print(input_params)

# Find the report template relative to this script's filesystem location
rmd_dir = dirname(sub("--file=", "", commandArgs(trailingOnly = FALSE)[4]))
rmd_path = paste0(rmd_dir, "/report.Rmd")
message(paste("Report template location:", rmd_path, sep = " "))

#path_prefix = fs::path_ext_remove(input_params$path_prefix)
out_path = paste0(input_params$path_prefix, "_AAV_report") # Adds .html and .pdf automatically
out_dir = dirname(out_path)
out_filename = basename(out_path)
message(paste("Output location:", out_path, sep = " "))

# Create a temporary directory with appropriate permissions
tmp_dir <- "/tmp/rmarkdown"
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
Sys.setenv(TMPDIR = tmp_dir)
# Set R temporary directory
options(tinytex.verbose = TRUE)
tempdir <- tmp_dir
message(paste("Tempdir:", tmp_dir, sep = " "))
message(paste("Tempdir:", tmp_dir, sep = " "))

listOutTmp <- system(paste("ls -al", shQuote(tmp_dir)), intern = TRUE)
# Display the output
cat(listOutTmp, sep = "\n")

# Copy all files from the Rmd directory to the temporary directory
files_to_copy <- list.files(rmd_dir, full.names = TRUE)
message(paste("Copying files from", rmd_dir, "to", tmp_dir))
copied_successfully <- file.copy(from = files_to_copy, to = tmp_dir, overwrite = TRUE, recursive = FALSE) # recursive=FALSE to copy files, not the dir itself if it were passed

if (!all(copied_successfully)) {
  failed_files <- files_to_copy[!copied_successfully]
  stop(paste("Failed to copy the following files to", tmp_dir, ":", paste(basename(failed_files), collapse=", ")))
}

# Define the path to the Rmd file within the temporary directory
tmp_rmd_path <- file.path(tmp_dir, basename(rmd_path))

rmarkdown::render(
  tmp_rmd_path, # Use the copied Rmd file in the temp directory
  output_format = "all",
  output_file = out_filename,
  output_dir = out_dir,
  knit_root_dir = working_dir,
  intermediates_dir = working_dir,

  params = input_params
)
