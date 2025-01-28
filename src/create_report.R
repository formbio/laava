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

rmarkdown::render(
  rmd_path,
  output_format = "all",
  output_file = out_filename,
  output_dir = out_dir,
  knit_root_dir = working_dir,
  params = input_params
)
