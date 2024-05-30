#!/usr/bin/env Rscript
library(rmarkdown)

working_dir = getwd()
message(paste("Working directory:", working_dir))

args = commandArgs(trailingOnly = TRUE)

rdata_path = args[1]

rmd_dir = dirname(sub("--file=", "", commandArgs(trailingOnly = FALSE)[4]))
rmd_path = paste0(rmd_dir, "/report.Rmd")
message(paste("Report template location:", rmd_path, sep = " "))

input_prefix = fs::path_ext_remove(rdata_path)
out_path = paste0(input_prefix, "_AAV_report") # Adds .html and .pdf automatically
out_dir = dirname(out_path)
out_filename = basename(out_path)
message(paste("Output location:", out_path, sep = " "))

input_params = list(
  rdata_path = rdata_path
)
message("Parameters:")
print(input_params)

rmarkdown::render(
  rmd_path,
  output_format = "all",
  output_file = out_filename,
  output_dir = out_dir,
  knit_root_dir = working_dir,
  params = input_params
)
