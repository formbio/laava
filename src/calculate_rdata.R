#!/usr/bin/env Rscript
library(tidyverse)

# Hush a benign message about "grouped output"
options(dplyr.summarise.inform = FALSE)


working_dir = getwd()
message(paste("Working directory:", working_dir))

args = commandArgs(trailingOnly = TRUE)

r_params = list()
r_params$input_prefix = args[1]
r_params$annot_filename = args[2] # ex: annotation.txt
r_params$sample_id = args[3]
r_params$vector_type = args[4]
r_params$flipflop_summary = ''
if (length(args) > 4) {
  r_params$flipflop_summary = args[5]
}


# Sort order
valid_types <- c('ssAAV', 'scAAV', 'host', 'repcap', 'helper', 'lambda', 'unmapped', 'chimeric')
valid_subtypes <- c('full', 'full-gap', 'left-partial', 'right-partial', 'wtITR-partial', 'mITR-partial', 'partial', 'backbone', 'vector+backbone')


# ----------------------------------------------------------
## Read the annotation file to find the vector target region
# ----------------------------------------------------------
# e.g.
#   NAME=myHost;TYPE=host;
#   NAME=myVector;TYPE=vector;REGION=1795-6553;
#   NAME=mRepCap;TYPE=repcap;REGION=1895-5987;

TARGET_REGION_START <- 0
TARGET_REGION_END <- 0
TARGET_REGION_START_REPCAP <- 0
TARGET_REGION_END_REPCAP <- 0

annot <- read.table(r_params$annot_filename)
for (i in 1:dim(annot)[1]) {
  if (unlist(strsplit(annot[i, ], ';'))[2] == 'TYPE=vector') {
    p <- unlist(strsplit(annot[i, ], ';'))[3];
    s_e <- as.integer(unlist(strsplit(unlist(strsplit(p, '='))[2], '-')));
    TARGET_REGION_START <- s_e[1];
    TARGET_REGION_END <- s_e[2];
  } else if (unlist(strsplit(annot[i, ], ';'))[2] == 'TYPE=repcap') {
    p <- unlist(strsplit(annot[i, ], ';'))[3];
    s_e <- as.integer(unlist(strsplit(unlist(strsplit(p, '='))[2], '-')));
    TARGET_REGION_START_REPCAP <- s_e[1];
    TARGET_REGION_END_REPCAP <- s_e[2];
  }
}


x.all.summary <- read_tsv(paste0(r_params$input_prefix, '.summary.tsv'), show_col_types = FALSE) %>%
  mutate(map_start = map_start0, map_end = map_end1) %>%
  mutate(SampleID = r_params$sample_id, .before = read_id)
write_tsv(x.all.summary, str_c(c(r_params$input_prefix, ".alignments.tsv"), collapse = ""))

x.all.err <- read_tsv(paste0(r_params$input_prefix, '.nonmatch_stat.tsv.gz'), show_col_types = FALSE) %>%
  mutate(SampleID = r_params$sample_id, .before = read_id)
x.all.read <- read_tsv(paste0(r_params$input_prefix, '.per_read.tsv'), show_col_types = FALSE) %>%
  mutate(SampleID = r_params$sample_id, .before = read_id)

x.all.err[x.all.err$type == 'D', "type"] <- 'deletion'
x.all.err[x.all.err$type == 'I', "type"] <- 'insertion'
x.all.err[x.all.err$type == 'X', "type"] <- 'mismatch'
x.all.err[x.all.err$type == 'N', "type"] <- 'gaps'


# ----------------------------------------------------------
# produce stats for vector only (ssAAV or scAAV)
# ----------------------------------------------------------

x.read.vector <- filter(x.all.read, assigned_type %in% c('scAAV', 'ssAAV'))
x.err.vector <- filter(x.all.err, read_id %in% x.read.vector$read_id)
x.summary.vector <- filter(x.all.summary, read_id %in% x.read.vector$read_id)


total_num_reads <- dim(x.read.vector)[1]

total_err <- dim(x.err.vector)[1]
x.err.vector$pos0_div <- (x.err.vector$pos0 %/% 10 * 10)
df.err.vector <- x.err.vector %>%
  group_by(pos0_div, type) %>%
  summarise(count = n())

x.err.vector$type_len_cat <- "1-10"
x.err.vector[x.err.vector$type_len > 10, "type_len_cat"] <- "11-100"
x.err.vector[x.err.vector$type_len > 100, "type_len_cat"] <- "100-500"
x.err.vector[x.err.vector$type_len > 500, "type_len_cat"] <- ">500"
x.err.vector$type_len_cat <- ordered(x.err.vector$type_len_cat, levels = c('1-10', '11-100', '100-500', '>500'))
write_tsv(x.err.vector, str_c(c(r_params$input_prefix, ".sequence-error.tsv"), collapse = ""))

x.read.vector$subtype <- x.read.vector$assigned_subtype
# Call "other" all complex and multi-part AAV alignments
x.read.vector[!x.read.vector$subtype %in% valid_subtypes, "subtype"] <- 'other'
# Call "snapback" any scAAV with an ITR-containing partial alignment but no full
# alignment
if (r_params$vector_type == "ssaav") {
  x.read.vector[((x.read.vector$assigned_type == "scAAV") &
    (x.read.vector$assigned_subtype == "left-partial")
  ), "subtype"] <- "left-snapback"
  x.read.vector[((x.read.vector$assigned_type == "scAAV") &
    (x.read.vector$assigned_subtype == "right-partial")
  ), "subtype"] <- "right-snapback"
}

total_read_count.vector <- sum(x.read.vector$effective_count)
df.read.vector1 <- x.read.vector %>%
  group_by(assigned_type) %>%
  summarise(e_count = sum(effective_count)) %>%
  mutate(freq = round(e_count * 100 / total_read_count.vector, 2))
df.read.vector1 <- df.read.vector1[order(-df.read.vector1$freq), ]


df.read.vector2 <- x.read.vector %>%
  group_by(assigned_type, assigned_subtype) %>%
  summarise(e_count = sum(effective_count)) %>%
  mutate(freq = round(e_count * 100 / total_read_count.vector, 2))
df.read.vector2 <- df.read.vector2[order(-df.read.vector2$freq), ]


x.all.read[is.na(x.all.read$assigned_type), "assigned_type"] <- 'unmapped'
x.all.read[grep("|", as.character(x.all.read$assigned_type), fixed = T), "assigned_type"] <- 'chimeric'
x.all.read[!(x.all.read$assigned_type %in% valid_types), "assigned_type"] <- 'other'
x.all.read[!(x.all.read$assigned_subtype %in% valid_subtypes), "assigned_subtype"] <- 'other'
write_tsv(x.all.read, str_c(c(r_params$input_prefix, ".readsummary.tsv"), collapse = ""))


# ----------------------------------------------------------
# Stats and plot for flip/flop analysis (if available)
# ----------------------------------------------------------

if (file.exists(r_params$flipflop_summary)) {
  df.read.ssaav <- dplyr::filter(df.read.vector2, assigned_type == 'ssAAV') %>%
    filter(
      assigned_subtype == 'full' |
        assigned_subtype == 'right-partial' |
        assigned_subtype == 'left-partial'
    ) %>%
    select(e_count) %>%
    as.data.frame()
  total_ssaav <- sum(df.read.ssaav$e_count)

  data.flipflop <- read.table(r_params$flipflop_summary,
    sep = '\t',
    header = T
  )
  df.flipflop <- data.flipflop %>%
    group_by(type, subtype, leftITR, rightITR) %>%
    summarise(count = n())
  scff <- filter(df.flipflop, type == 'scAAV')
  ssff <- filter(df.flipflop, type == 'ssAAV')

  # Double single-stranded counts if appropriate
  numssff <- sum(ssff$count)
  if (is.numeric(numssff) & is.numeric(total_ssaav)) {
    if (total_ssaav > 0 & numssff > 0 & numssff * 2 == total_ssaav) {
      ssff <- ssff %>% mutate(count = count * 2)
    } else {
      ssff <- ssff %>% mutate(count = count)
    }
  }
  # Write TSV of flip flop configurations
  fftbl <- bind_rows(scff, ssff)
  write_tsv(fftbl, str_c(c(r_params$input_prefix, ".flipflop.tsv"), collapse = ""))
}


## For downstream consumers

save.image(file = paste0(r_params$input_prefix, ".Rdata"))
