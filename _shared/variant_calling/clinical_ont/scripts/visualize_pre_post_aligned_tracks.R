#!/usr/bin/env Rscript

# Single linear genome track (post sample on pre reference coordinates)
# with gene bands and lollipop mutation markers.
#
# Outputs:
# - figures/Figure_pre_post_aligned_tracks_lollipop.png

needed <- c("readr", "dplyr", "ggplot2", "stringr", "tidyr")
for (pkg in needed) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tidyr)
})

input_tsv <- "tables/mutations.tsv"
ref_fasta <- "assemblies/pre.polished.fasta"
gff_path <- "annotation/snpeff/data/mrsa_ref_pre/genes.gff"
out_png <- "figures/Figure_pre_post_aligned_tracks_lollipop.png"
dir.create("figures", recursive = TRUE, showWarnings = FALSE)

read_contig_len <- function(fasta_path, contig_id = "contig_2") {
  x <- readLines(fasta_path, warn = FALSE)
  hdr <- which(startsWith(x, ">"))
  if (length(hdr) == 0) stop("No FASTA headers found: ", fasta_path)
  nm <- sub("^>(\\S+).*", "\\1", x[hdr])
  idx <- which(nm == contig_id)
  if (length(idx) == 0) stop("Contig not found in FASTA: ", contig_id)
  h <- hdr[idx[1]]
  next_h <- if (idx[1] < length(hdr)) hdr[idx[1] + 1] - 1 else length(x)
  seq_lines <- x[(h + 1):next_h]
  seq_lines <- seq_lines[nzchar(seq_lines)]
  sum(nchar(seq_lines))
}

read_gene_ranges <- function(path) {
  if (!file.exists(path)) return(tibble())
  read_tsv(
    path,
    col_names = c("chrom", "source", "type", "start", "end", "score", "strand", "phase", "attrs"),
    comment = "#",
    show_col_types = FALSE,
    col_types = cols(.default = col_character())
  ) %>%
    filter(chrom == "contig_2") %>%
    transmute(
      locus_tag = str_match(attrs, "locus_tag=([^;]+)")[, 2],
      start = as.numeric(start),
      end = as.numeric(end)
    ) %>%
    filter(!is.na(locus_tag), !is.na(start), !is.na(end)) %>%
    mutate(start = pmin(start, end), end = pmax(start, end)) %>%
    distinct(locus_tag, .keep_all = TRUE)
}

effect_class <- function(effect_txt) {
  e <- str_to_lower(coalesce(effect_txt, ""))
  case_when(
    str_detect(e, "missense") ~ "missense",
    str_detect(e, "upstream") ~ "upstream",
    str_detect(e, "synonymous") ~ "synonymous",
    TRUE ~ "other"
  )
}

derive_pheno_category <- function(gene_id, effect_txt, annotation_txt) {
  g <- str_to_lower(coalesce(gene_id, ""))
  e <- str_to_lower(coalesce(effect_txt, ""))
  a <- str_to_lower(coalesce(annotation_txt, ""))
  case_when(
    str_detect(g, "debjfdcl_01664") ~ "regulator",
    str_detect(g, "debjfdcl_00876") | str_detect(g, "debjfdcl_01461") ~ "metabolism",
    str_detect(g, "debjfdcl_02530") ~ "envelope_transport",
    str_detect(g, "debjfdcl_00981") ~ "surface_access",
    str_detect(e, "rrna") | str_detect(a, "rrna") ~ "rrna_adjacent",
    TRUE ~ "other"
  )
}

mut <- read_tsv(input_tsv, show_col_types = FALSE) %>%
  filter(chrom == "contig_2") %>%
  mutate(
    pos = as.numeric(pos),
    effect_plot = effect_class(effect_detail),
    pheno_category = derive_pheno_category(gene, effect_detail, annotation),
    label = case_when(
      str_detect(gene, "DEBJFDCL_01664") ~ "agrA",
      str_detect(gene, "DEBJFDCL_00876") ~ "yclM",
      str_detect(gene, "DEBJFDCL_01461") ~ "fumC-upstream",
      str_detect(gene, "DEBJFDCL_02530") ~ "mdrP",
      str_detect(gene, "DEBJFDCL_00981") ~ "ebh",
      TRUE ~ ""
    )
  )

if (nrow(mut) == 0) stop("No contig_2 entries in mutations.tsv")

genome_len <- read_contig_len(ref_fasta, "contig_2")
gene_ranges <- read_gene_ranges(gff_path)

track_y <- 0.0
track_half_h <- 0.06

pts <- mut %>%
  filter(present_in_post %in% TRUE) %>%
  transmute(pos, effect_plot, pheno_category, label) %>%
  arrange(pos) %>%
  mutate(
    lbl_idx = row_number(),
    alt_dir = if_else(lbl_idx %% 2 == 0L, -1L, 1L),
    y_track = track_y,
    y_head = y_track + if_else(label != "", alt_dir * 0.18, 0.12),
    y_text = y_track + if_else(label != "", alt_dir * 0.27, 0.2)
  ) %>%
  mutate(
    effect_plot = factor(effect_plot, levels = c("missense", "upstream", "synonymous", "other")),
    fill_key = if_else(label != "", "relevant_black", "muted"),
    y_edge = if_else(y_head >= y_track, y_track + track_half_h, y_track - track_half_h)
  )

# Color-coded gene regions on each track (phenotype-coded, unlabeled).
# Use a display window large enough to be visible at whole-genome scale.
region_window <- 20000
gene_regions <- mut %>%
  filter(label != "") %>%
  transmute(
    locus_tag = gene,
    pheno_category,
    pos,
    gene_label = label
  ) %>%
  left_join(gene_ranges, by = "locus_tag") %>%
  mutate(
    start = if_else(is.na(start), pmax(0, pos - region_window), start),
    end = if_else(is.na(end), pmin(genome_len, pos + region_window), end)
  ) %>%
  select(pheno_category, start, end, gene_label) %>%
  distinct() %>%
  transmute(
    pheno_category,
    gene_label,
    start, end,
    y = track_y,
    band_ymin = y - track_half_h,
    band_ymax = y + track_half_h
  ) %>%
  distinct()

cat("Mutation count (post track):\n")
print(nrow(pts))
cat("\nMarker preview:\n")
print(pts %>% arrange(pos) %>% head(15))

cat_palette <- c(
  regulator = "#E69F00",
  metabolism = "#009E73",
  envelope_transport = "#0072B2",
  surface_access = "#6c757d",
  rrna_adjacent = "#CC79A7",
  other = "#9a9a9a",
  muted = "#d4d7dc",
  relevant_black = "#1f1f1f"
)
shape_map <- c(missense = 24, upstream = 22, synonymous = 21, other = 23)

track_df <- tibble(y = track_y)

p <- ggplot() +
  geom_rect(
    data = track_df,
    aes(xmin = 0, xmax = genome_len, ymin = y - track_half_h, ymax = y + track_half_h),
    fill = "#e8e8e8",
    color = NA
  ) +
  geom_rect(
    data = gene_regions,
    aes(
      xmin = start, xmax = end,
      ymin = band_ymin, ymax = band_ymax,
      fill = pheno_category
    ),
    color = NA,
    alpha = 0.95
  ) +
  geom_point(
    data = pts,
    aes(
      x = pos,
      y = y_head,
      shape = effect_plot,
      fill = fill_key,
      color = fill_key
    ),
    size = 2.1,
    stroke = 0.55,
    show.legend = TRUE
  ) +
  geom_segment(
    data = pts,
    aes(
      x = pos, xend = pos,
      y = y_edge, yend = y_head,
      color = fill_key
    ),
    linewidth = 0.8,
    alpha = 0.95
  ) +
  geom_text(
    data = filter(pts, label != "", y_text > y_track),
    aes(x = pos, y = y_text, label = label),
    size = 3.1,
    hjust = 0.5,
    vjust = 0.15,
    lineheight = 0.9,
    color = "#1a1a1a",
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  geom_text(
    data = filter(pts, label != "", y_text <= y_track),
    aes(x = pos, y = y_text, label = label),
    size = 3.1,
    hjust = 0.5,
    vjust = 1,
    lineheight = 0.9,
    color = "#1a1a1a",
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  scale_shape_manual(values = shape_map, name = "Mutation type", drop = FALSE) +
  scale_fill_manual(values = cat_palette, guide = "none") +
  scale_color_manual(values = cat_palette, guide = "none") +
  scale_x_continuous(
    limits = c(0, genome_len * 1.01),
    breaks = seq(0, ceiling(genome_len / 5e5) * 5e5, by = 5e5),
    labels = function(x) sprintf("%.1f Mb", x / 1e6),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = track_y,
    labels = "Post (pre coordinates)",
    limits = c(-0.38, 0.38),
    expand = c(0, 0)
  ) +
  labs(
    title = "Post-treatment variants along pre reference (contig_2)",
    subtitle = "Gene bands: functional category; lollipop shape: mutation class",
    x = "Genomic position",
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  ) +
  guides(shape = guide_legend(order = 1))

ggsave(out_png, p, width = 9, height = 4.6, dpi = 300)
cat("\nSaved:", out_png, "\n")
