#!/usr/bin/env Rscript

# Simple linear genome map for MRSA reference assembly.
# Output:
# - figures/Figure_linear_genome_map.png

needed <- c("dplyr", "ggplot2", "readr", "stringr")
for (pkg in needed) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
})

fasta_path <- "assemblies/pre.polished.fasta"
gff_path <- "annotation/snpeff/data/mrsa_ref_pre/genes.gff"
out_png <- "figures/Figure_linear_genome_map.png"
out_simple_png <- "figures/Figure_linear_genome_map_simple.png"

dir.create("figures", recursive = TRUE, showWarnings = FALSE)

read_fasta_lengths <- function(path) {
  if (!file.exists(path)) stop("FASTA not found: ", path)
  x <- readLines(path, warn = FALSE)
  hdr_idx <- which(startsWith(x, ">"))
  if (length(hdr_idx) == 0) stop("No FASTA headers found in: ", path)

  contig_names <- character(length(hdr_idx))
  contig_lens <- numeric(length(hdr_idx))

  for (i in seq_along(hdr_idx)) {
    h <- hdr_idx[i]
    next_h <- if (i < length(hdr_idx)) hdr_idx[i + 1] - 1 else length(x)
    contig_names[i] <- sub("^>(\\S+).*", "\\1", x[h])
    seq_lines <- x[(h + 1):next_h]
    seq_lines <- seq_lines[nzchar(seq_lines)]
    contig_lens[i] <- sum(nchar(seq_lines))
  }

  tibble(contig = contig_names, length_bp = contig_lens)
}

read_key_features <- function(path, main_contig) {
  if (!file.exists(path)) return(tibble())
  g <- read_tsv(
    path,
    col_names = c("chrom", "source", "type", "start", "end", "score", "strand", "phase", "attrs"),
    comment = "#",
    show_col_types = FALSE,
    col_types = cols(.default = col_character())
  ) %>%
    filter(chrom == main_contig) %>%
    mutate(
      start = as.numeric(start),
      end = as.numeric(end),
      gene_symbol = coalesce(
        str_match(attrs, "gene=([^;]+)")[, 2],
        str_match(attrs, "Name=([^;]+)")[, 2],
        str_match(attrs, "locus_tag=([^;]+)")[, 2]
      ),
      is_rrna = type == "rRNA" | str_detect(str_to_lower(attrs), "rrna")
    )

  key_genes <- c("agrA", "yclM", "fumC", "mdrP", "ebh")
  key <- g %>%
    filter(str_to_lower(gene_symbol) %in% str_to_lower(key_genes)) %>%
    transmute(
      feature = gene_symbol,
      start = pmin(start, end),
      end = pmax(start, end),
      feature_type = "key_gene"
    ) %>%
    distinct(feature, .keep_all = TRUE)

  rrna <- g %>%
    filter(is_rrna) %>%
    transmute(
      feature = "16S/5S rRNA loci",
      start = pmin(start, end),
      end = pmax(start, end),
      feature_type = "rRNA"
    ) %>%
    arrange(start) %>%
    mutate(
      prev_end = lag(end, default = first(start)),
      new_grp = if_else(start - prev_end > 20000, 1L, 0L),
      grp = cumsum(new_grp)
    ) %>%
    group_by(grp) %>%
    summarise(
      feature = "16S/5S rRNA loci",
      start = min(start, na.rm = TRUE),
      end = max(end, na.rm = TRUE),
      feature_type = "rRNA",
      .groups = "drop"
    )

  bind_rows(key, rrna)
}

contigs <- read_fasta_lengths(fasta_path) %>%
  arrange(desc(length_bp)) %>%
  mutate(
    contig_f = factor(contig, levels = rev(contig)),
    label = sprintf("%s (%.2f Mb)", contig, length_bp / 1e6)
  )

main_contig <- contigs$contig[1]
features <- read_key_features(gff_path, main_contig)

cat("Contig summary:\n")
print(contigs)
cat("\nMain-contig feature bands preview:\n")
print(features)

base_plot <- ggplot(contigs, aes(x = 0, xend = length_bp, y = contig_f, yend = contig_f)) +
  geom_segment(linewidth = 8, color = "#9ea09b", lineend = "round")

if (nrow(features) > 0) {
  feature_plot <- ggplot(features) +
    geom_segment(
      aes(x = start, xend = end, y = 1.15, yend = 1.15, color = feature_type),
      linewidth = 5, lineend = "round"
    ) +
    geom_text(
      aes(x = (start + end) / 2, y = 1.27, label = feature, color = feature_type),
      size = 3
    ) +
    scale_color_manual(values = c(key_gene = "#1b9e77", rRNA = "#cc79a7"), guide = "none")
} else {
  feature_plot <- NULL
}

p <- base_plot +
  labs(
    title = "Linear genome map of MRSA reference assembly",
    x = "Genomic position",
    y = NULL
  ) +
  scale_x_continuous(
    labels = function(x) sprintf("%.1f Mb", x / 1e6),
    expand = c(0.01, 0.01)
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(face = "bold")
  )

# Overlay feature bands on the top/main contig only.
if (!is.null(feature_plot)) {
  feature_text <- features %>%
    filter(feature_type == "key_gene") %>%
    mutate(xmid = (start + end) / 2)
  rrna_text <- features %>%
    filter(feature_type == "rRNA") %>%
    summarise(xmid = min(start, na.rm = TRUE), .groups = "drop") %>%
    mutate(feature = "rRNA loci", feature_type = "rRNA")

  p <- p +
    geom_segment(
      data = features,
      aes(x = start, xend = end, y = max(as.numeric(contigs$contig_f)) + 0.32, yend = max(as.numeric(contigs$contig_f)) + 0.32, color = feature_type),
      linewidth = 4, lineend = "round", inherit.aes = FALSE
    ) +
    geom_text(
      data = feature_text,
      aes(x = xmid, y = max(as.numeric(contigs$contig_f)) + 0.55, label = feature, color = feature_type),
      size = 3, inherit.aes = FALSE
    ) +
    geom_text(
      data = rrna_text,
      aes(x = xmid, y = max(as.numeric(contigs$contig_f)) + 0.55, label = feature, color = feature_type),
      size = 3, inherit.aes = FALSE
    ) +
    scale_color_manual(values = c(key_gene = "#1b9e77", rRNA = "#cc79a7"), name = "Top-track feature")
}

ggsave(out_png, p, width = 7.2, height = 4.2, dpi = 300)
cat("\nSaved:", out_png, "\n")

# Minimal map-only version (no top feature track).
p_simple <- base_plot +
  labs(
    title = "Linear genome map of MRSA reference assembly",
    x = "Genomic position",
    y = NULL
  ) +
  scale_x_continuous(
    labels = function(x) sprintf("%.1f Mb", x / 1e6),
    breaks = seq(0, ceiling(max(contigs$length_bp) / 5e5) * 5e5, by = 5e5),
    expand = c(0.01, 0.01)
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(face = "bold")
  )

ggsave(out_simple_png, p_simple, width = 7.2, height = 3.6, dpi = 300)
cat("Saved:", out_simple_png, "\n")
