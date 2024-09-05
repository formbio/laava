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
message("Parameters:")
print(r_params)


# Sort order
valid_types <- c('ssAAV', 'scAAV', 'host', 'repcap', 'helper', 'lambda', 'unmapped', 'chimeric')
valid_subtypes <- c('full', 'full-gap', 'left-partial', 'right-partial', 'wtITR-partial', 'mITR-partial', 'partial', 'backbone', 'vector+backbone')



# =========================
# annotation.txt -> TARGET_REGION_{START,END}(_REPCAP)
# -------------------------
# Read the annotation file to find the vector target region
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


# =====================================
# alignments.tsv.gz
# -------------------------------------

x.all.summary <- read_tsv(paste0(r_params$input_prefix, '.alignments.tsv.gz'), show_col_types = FALSE)


# ---------------------------------------
# per_read.tsv.gz
# ---------------------------------------

x.all.read <- read_tsv(paste0(r_params$input_prefix, '.per_read.tsv.gz'), show_col_types = FALSE)

# XXX only in Rdata
total_read_count_all <- sum(x.all.read$effective_count) # dim(x.all.read)[1]
# "Assigned types by read alignment characteristics"
df.read1 <- x.all.read %>%
  group_by(reference_label, assigned_type) %>%
  summarise(e_count = sum(effective_count)) %>%
  mutate(freq = round(e_count * 100 / total_read_count_all, 2))
df.read1 <- df.read1[order(df.read1$reference_label, df.read1$freq, decreasing=TRUE), ]

# Filter to ssAAV/scAAV vector only
# NB: also used by nonmatch.tsv.gz below
# XXX this df is only in Rdata, for report display
x.read.vector <- filter(x.all.read, reference_label == "vector")

# XXX only in Rdata
total_read_count_vector <- sum(x.read.vector$effective_count)
df.read.vector1 <- x.read.vector %>%
  group_by(assigned_type) %>%
  summarise(e_count = sum(effective_count)) %>%
  mutate(freq = round(e_count * 100 / total_read_count_vector, 2))
df.read.vector1 <- df.read.vector1[order(-df.read.vector1$freq), ]

# XXX only in Rdata; also used in flipflop code below
df.read.vector2 <- x.read.vector %>%
  group_by(assigned_type, assigned_subtype) %>%
  summarise(e_count = sum(effective_count)) %>%
  mutate(freq = round(e_count * 100 / total_read_count_vector, 2))
df.read.vector2 <- df.read.vector2[order(-df.read.vector2$freq), ]

# XXX only in Rdata, for the report plots -- but source is alignments.tsv
x.summary.primary = filter(x.all.summary, is_mapped == "Y", is_supp == "N")
x.joined.vector = left_join(x.read.vector, x.summary.primary, by = "read_id", multiple = "first")


# ==================================================
# nonmatch.tsv.gz
# --------------------------------------------------

x.all.err <- read_tsv(paste0(r_params$input_prefix, '.nonmatch.tsv.gz'), show_col_types = FALSE)

# Filter for ss/scAAV vector only
x.err.vector <- filter(x.all.err, read_id %in% x.read.vector$read_id)

# Spell out nonmatch types -- from CIGAR codes to words
x.err.vector[x.err.vector$type == 'D', "type"] <- 'deletion'
x.err.vector[x.err.vector$type == 'I', "type"] <- 'insertion'
x.err.vector[x.err.vector$type == 'X', "type"] <- 'mismatch'
x.err.vector[x.err.vector$type == 'N', "type"] <- 'gaps'

# Round positions down to 10s for plotting
x.err.vector$pos0_div <- (x.err.vector$pos0 %/% 10 * 10)

# XXX only in Rdata
df.err.vector <- x.err.vector %>%
  group_by(pos0_div, type) %>%
  summarise(count = n())


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
  write_tsv(fftbl, paste0(r_params$input_prefix, ".flipflop.tsv"))
}


## For downstream consumers

save.image(file = paste0(r_params$input_prefix, ".Rdata"))
