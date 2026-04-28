#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript earlgrey_like_plots.R repeatmasker.out genome_size output_prefix")
}

rm_out <- args[1]
genome_size <- as.numeric(args[2])
prefix <- args[3]

dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)

read_repeatmasker_out <- function(path) {
  x <- readLines(path)
  x <- x[!grepl("^\\s*$", x)]
  x <- x[!grepl("^\\s*SW\\s+perc", x)]
  x <- x[!grepl("^\\s*score\\s+div", x)]
  x <- x[!grepl("^\\s*-+", x)]

  df <- read.table(
    text = x,
    fill = TRUE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )

  if (ncol(df) < 14) {
    stop("RepeatMasker .out file does not look like standard RepeatMasker output")
  }

  colnames(df)[1:14] <- c(
    "sw_score", "perc_div", "perc_del", "perc_ins",
    "query", "q_start", "q_end", "q_left",
    "strand", "repeat_name", "repeat_class",
    "r_start", "r_end", "r_left"
  )

  df %>%
    mutate(
      q_start = as.numeric(q_start),
      q_end = as.numeric(q_end),
      perc_div = as.numeric(perc_div),
      length_bp = abs(q_end - q_start) + 1,
      repeat_class = ifelse(repeat_class == "" | is.na(repeat_class), "Unknown", repeat_class),
      class_major = case_when(
        grepl("^DNA", repeat_class, ignore.case = TRUE) ~ "DNA",
        grepl("^LINE", repeat_class, ignore.case = TRUE) ~ "LINE",
        grepl("^SINE", repeat_class, ignore.case = TRUE) ~ "SINE",
        grepl("^LTR", repeat_class, ignore.case = TRUE) ~ "LTR",
        grepl("Penelope", repeat_class, ignore.case = TRUE) ~ "Penelope",
        grepl("RC|Helitron|Rolling", repeat_class, ignore.case = TRUE) ~ "Rolling-circle",
        grepl("Satellite", repeat_class, ignore.case = TRUE) ~ "Satellite",
        grepl("Simple_repeat", repeat_class, ignore.case = TRUE) ~ "Simple_repeat",
        grepl("Low_complexity", repeat_class, ignore.case = TRUE) ~ "Low_complexity",
        TRUE ~ "Other"
      )
    )
}

df <- read_repeatmasker_out(rm_out)

summary_class <- df %>%
  group_by(class_major) %>%
  summarise(
    coverage_bp = sum(length_bp, na.rm = TRUE),
    copy_number = n(),
    genome_percent = 100 * coverage_bp / genome_size,
    .groups = "drop"
  ) %>%
  arrange(desc(coverage_bp))

write_tsv(summary_class, paste0(prefix, "_class_summary.tsv"))

p1 <- ggplot(summary_class, aes(x = reorder(class_major, genome_percent), y = genome_percent)) +
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(
    x = NULL,
    y = "% genome masked",
    title = "TE genome coverage by major class"
  )

ggsave(paste0(prefix, "_coverage_by_class.pdf"), p1, width = 7, height = 5)

p2 <- ggplot(summary_class, aes(x = reorder(class_major, copy_number), y = copy_number)) +
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(
    x = NULL,
    y = "Copy number",
    title = "TE copy number by major class"
  )

ggsave(paste0(prefix, "_copy_number_by_class.pdf"), p2, width = 7, height = 5)

landscape <- df %>%
  filter(!is.na(perc_div)) %>%
  mutate(div_bin = floor(perc_div)) %>%
  group_by(class_major, div_bin) %>%
  summarise(bp = sum(length_bp, na.rm = TRUE), .groups = "drop") %>%
  mutate(genome_percent = 100 * bp / genome_size)

write_tsv(landscape, paste0(prefix, "_repeat_landscape.tsv"))

p3 <- ggplot(landscape, aes(x = div_bin, y = genome_percent, fill = class_major)) +
  geom_col(width = 1) +
  theme_bw(base_size = 12) +
  labs(
    x = "Kimura divergence / RepeatMasker divergence (%)",
    y = "% genome masked",
    fill = "Class",
    title = "Repeat landscape"
  )

ggsave(paste0(prefix, "_repeat_landscape_stacked.pdf"), p3, width = 9, height = 5)

p4 <- ggplot(landscape, aes(x = div_bin, y = genome_percent)) +
  geom_col(width = 1) +
  facet_wrap(~ class_major, scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(
    x = "Kimura divergence / RepeatMasker divergence (%)",
    y = "% genome masked",
    title = "Repeat landscape by class"
  )

ggsave(paste0(prefix, "_repeat_landscape_faceted.pdf"), p4, width = 10, height = 7)

top_families <- df %>%
  group_by(repeat_name, repeat_class, class_major) %>%
  summarise(
    coverage_bp = sum(length_bp, na.rm = TRUE),
    copy_number = n(),
    genome_percent = 100 * coverage_bp / genome_size,
    .groups = "drop"
  ) %>%
  arrange(desc(coverage_bp)) %>%
  slice_head(n = 30)

write_tsv(top_families, paste0(prefix, "_top_families.tsv"))

p5 <- ggplot(top_families, aes(x = reorder(repeat_name, genome_percent), y = genome_percent)) +
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 10) +
  labs(
    x = NULL,
    y = "% genome masked",
    title = "Top TE families by genome coverage"
  )

ggsave(paste0(prefix, "_top_families.pdf"), p5, width = 8, height = 8)
