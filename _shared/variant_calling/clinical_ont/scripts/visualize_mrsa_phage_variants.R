#!/usr/bin/env Rscript

# Publication-ready MRSA pre/post variant visualizations.
# Resolves locus_tag -> gene symbol via reference GFF so categories and labels
# match biology (agrA, scpB, etc.) instead of DEBJFDCL_* IDs.
#
# Outputs:
# - figures/Figure_genomic_variants.png
# - tables/variants_candidates_for_figure.tsv

needed <- c("readr", "dplyr", "ggplot2", "ggrepel", "patchwork", "stringr")
for (pkg in needed) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(stringr)
})

input_tsv <- "tables/variants_for_paper.tsv"
gff_default <- "annotation/snpeff/data/mrsa_ref_pre/genes.gff"
out_png <- "figures/Figure_genomic_variants.png"
out_jpg <- "figures/Figure_genomic_variants.jpg"
out_candidates <- "tables/variants_candidates_for_figure.tsv"
out_panelA_main_png <- "figures/Figure_genomic_variants_PanelA_main.png"
out_panelA_main_jpg <- "figures/Figure_genomic_variants_PanelA_main.jpg"
out_panelA_mini_png <- "figures/Figure_genomic_variants_PanelA_mini.png"
out_panelA_mini_jpg <- "figures/Figure_genomic_variants_PanelA_mini.jpg"
out_panelB_png <- "figures/Figure_genomic_variants_PanelB.png"
out_panelB_jpg <- "figures/Figure_genomic_variants_PanelB.jpg"

# Key loci for emphasis (match resolved gene_symbol from GFF; raw TSV `gene` is often locus_tag).
highlight_gene_symbols <- c("agrA", "yclM", "fumC", "mdrP", "ebh")

dir.create("figures", recursive = TRUE, showWarnings = FALSE)
dir.create("tables", recursive = TRUE, showWarnings = FALSE)

# --- GFF: locus_tag -> gene symbol (Prokka gene= / Name=) -----------------
read_locus_gene_map <- function(gff_path) {
  if (!file.exists(gff_path)) {
    warning("GFF not found; gene symbols will be missing: ", gff_path)
    return(tibble(locus_tag = character(), gene_symbol = character()))
  }
  gff <- read_tsv(
    gff_path,
    col_names = c(
      "seqid", "source", "type", "start", "end",
      "score", "strand", "phase", "attrs"
    ),
    comment = "#",
    show_col_types = FALSE,
    col_types = cols(.default = col_character())
  )
  gff %>%
    filter(!is.na(attrs)) %>%
    transmute(
      locus_tag = str_match(attrs, "locus_tag=([^;]+)")[, 2],
      gene_attr = str_match(attrs, "gene=([^;]+)")[, 2],
      name_attr = str_match(attrs, "Name=([^;]+)")[, 2],
      gene_symbol = coalesce(gene_attr, name_attr, locus_tag)
    ) %>%
    filter(!is.na(locus_tag)) %>%
    distinct(locus_tag, .keep_all = TRUE) %>%
    select(locus_tag, gene_symbol)
}

read_gene_features <- function(gff_path) {
  if (!file.exists(gff_path)) {
    return(tibble(
      seqid = character(),
      start = numeric(),
      end = numeric(),
      strand = character(),
      locus_tag = character(),
      gene_symbol = character(),
      product_annot = character()
    ))
  }
  read_tsv(
    gff_path,
    col_names = c(
      "seqid", "source", "type", "start", "end",
      "score", "strand", "phase", "attrs"
    ),
    comment = "#",
    show_col_types = FALSE,
    col_types = cols(.default = col_character())
  ) %>%
    filter(type == "CDS", !is.na(attrs)) %>%
    transmute(
      seqid = seqid,
      start = as.numeric(start),
      end = as.numeric(end),
      strand = strand,
      locus_tag = str_match(attrs, "locus_tag=([^;]+)")[, 2],
      gene_attr = str_match(attrs, "gene=([^;]+)")[, 2],
      name_attr = str_match(attrs, "Name=([^;]+)")[, 2],
      product_annot = str_match(attrs, "product=([^;]+)")[, 2],
      gene_symbol = coalesce(gene_attr, name_attr, locus_tag)
    ) %>%
    filter(!is.na(locus_tag))
}

shorten_text <- function(x, maxl = 24) {
  out <- x
  long <- !is.na(out) & nchar(out) > maxl
  out[long] <- paste0(str_sub(out[long], 1, maxl - 3), "...")
  out
}

# Categories chosen to match manuscript narrative; tweak keywords as needed.
derive_plot_category <- function(gene_sym, product, effect) {
  g <- str_to_lower(coalesce(gene_sym, ""))
  p <- str_to_lower(coalesce(product, ""))
  e <- str_to_lower(coalesce(effect, ""))

  case_when(
    g %in% c("agra", "cdar", "agrb", "agrc", "agrd") ~ "regulator",
    g %in% c("yclm", "fumc") ~ "metabolism",
    str_detect(e, "rrna") | str_detect(p, "rrna") ~ "rRNA_adjacent",
    g == "mdrp" |
      str_detect(p, "antiporter|transporter|efflux|symporter|channel|porin|major facilitator|abc transporter") ~ "envelope/transport",
    str_detect(p, "aspartokinase|fumarate|hydratase|hydratase|synthase|kinase|dehydrogenase|reductase") &
      !str_detect(p, "hypothetical") ~ "metabolism",
    TRUE ~ "other"
  )
}

derive_pheno_cluster <- function(gene_sym, product, effect) {
  g <- str_to_lower(coalesce(gene_sym, ""))
  p <- str_to_lower(coalesce(product, ""))
  e <- str_to_lower(coalesce(effect, ""))
  case_when(
    g %in% c("agra", "agrb", "agrc", "agrd") ~ "Quorum/Virulence regulation",
    g == "mdrp" | str_detect(p, "antiporter|transporter|efflux|symporter|channel|porin") ~ "Ion homeostasis/Membrane transport",
    g == "ebh" ~ "Surface architecture/access",
    g %in% c("yclm", "fumc") | str_detect(p, "aspartokinase|fumarate|hydratase|dehydrogenase|kinase|reductase|synthase") ~ "Cell-wall precursor/envelope metabolism",
    str_detect(e, "rrna") | str_detect(p, "rrna") ~ "rRNA-adjacent/structural context",
    TRUE ~ "Other coding background"
  )
}

mini_protein_label <- function(gene_sym, effect_txt) {
  g <- str_to_lower(coalesce(gene_sym, ""))
  e <- str_to_lower(coalesce(effect_txt, ""))
  case_when(
    g == "agra" ~ "AgrA accessory regulator",
    g == "yclm" ~ "YclM aspartokinase 3",
    g == "fumc" & str_detect(e, "upstream") ~ "fumC upstream regulatory",
    g == "fumc" ~ "FumC fumarate hydratase",
    g == "mdrp" ~ "MdrP Na Li K H antiporter",
    g == "ebh" ~ "Ebh surface protein",
    TRUE ~ gene_sym
  )
}

gff_path <- Sys.getenv("MRSA_VARIANTS_GFF", unset = gff_default)
locus_map <- read_locus_gene_map(gff_path)
gene_features <- read_gene_features(gff_path)

# Post-only rows for key gene symbols (paper TSV may omit some Tier2 loci still in mutations.tsv).
supplement_variants_from_mutations <- function(
    variants_tbl,
    mut_path,
    locus_map_df,
    key_syms
) {
  if (!file.exists(mut_path)) {
    return(variants_tbl)
  }
  mut <- read_tsv(mut_path, show_col_types = FALSE)
  if (!"unique_to_post" %in% names(mut)) {
    return(variants_tbl)
  }
  mut <- mut %>%
    filter(unique_to_post %in% TRUE) %>%
    left_join(locus_map_df, by = c("gene" = "locus_tag")) %>%
    mutate(
      gsym = coalesce(gene_symbol, gene),
      gsym_l = str_to_lower(gsym),
      eff_l = str_to_lower(coalesce(effect_detail, effect, ""))
    )
  key_l <- str_to_lower(key_syms)
  mut_f <- mut %>%
    filter(
      gsym_l %in% key_l |
        (gsym_l == "fumc" & str_detect(eff_l, "upstream"))
    )
  if (nrow(mut_f) == 0) {
    return(variants_tbl)
  }
  paper_cols <- names(variants_tbl)
  ex <- mut_f %>%
    transmute(
      chrom = chrom,
      pos = as.numeric(pos),
      ref = as.character(ref),
      alt = as.character(alt),
      gene = gene,
      locus_tag = gene,
      product = coalesce(product, ""),
      effect = coalesce(effect_detail, ""),
      impact = case_when(
        str_detect(eff_l, "missense") ~ "MODERATE",
        str_detect(eff_l, "synonymous") ~ "LOW",
        TRUE ~ "MODIFIER"
      ),
      category = "other",
      tier = if_else(as.integer(pos) == 1646692L, "Tier 1", "Tier 2"),
      focus_set = "yes",
      caller_support = "mutations_supplement",
      pass_both_callers = NA,
      medaka_qual = as.numeric(coalesce(qual, NA_real_)),
      medaka_dp = NA_real_,
      medaka_af = NA_real_,
      clair3_qual = NA_real_,
      clair3_dp = NA_real_,
      clair3_af = NA_real_,
      region_flag = "other"
    )
  for (cn in setdiff(paper_cols, names(ex))) {
    ex[[cn]] <- NA
  }
  ex <- ex[, paper_cols]
  bind_rows(variants_tbl, ex) %>%
    distinct(chrom, pos, ref, alt, .keep_all = TRUE)
}

variants <- read_tsv(input_tsv, show_col_types = FALSE)
variants <- supplement_variants_from_mutations(
  variants,
  "tables/mutations.tsv",
  locus_map,
  highlight_gene_symbols
)
cat("Columns:\n")
print(colnames(variants))
cat("\nPreview (up to 10 rows):\n")
print(head(variants, 10))

variants_chr <- variants %>%
  filter(chrom == "contig_2") %>%
  mutate(
    pos = as.numeric(pos),
    locus_tag = if ("locus_tag" %in% names(.)) coalesce(locus_tag, gene) else gene,
    tier_raw = tier,
    tier = factor(
      str_replace_all(as.character(tier), "\\s+", ""),
      levels = c("Tier1", "Tier2", "Tier3"),
      ordered = TRUE
    )
  ) %>%
  left_join(locus_map, by = "locus_tag") %>%
  mutate(
    gene_symbol = coalesce(gene_symbol, locus_tag),
    category_in = coalesce(category, ""),
    plot_category = if_else(
      !is.na(category_in) & nzchar(category_in) & category_in != "other",
      category_in,
      derive_plot_category(gene_symbol, product, effect)
    ),
    pheno_cluster = derive_pheno_cluster(gene_symbol, product, effect),
    plot_category = factor(
      plot_category,
      levels = c(
        "regulator", "envelope/transport", "metabolism",
        "rRNA_adjacent", "other"
      )
    ),
    gene_label = case_when(
      str_to_lower(coalesce(gene_symbol, "")) == "fumc" &
        str_detect(str_to_lower(coalesce(effect, "")), "upstream") ~ "fumC (upstream)",
      !is.na(gene_symbol) &
        nzchar(gene_symbol) &
        !str_detect(gene_symbol, "^DEBJFDCL_") ~ gene_symbol,
      !is.na(product) & nzchar(product) ~ shorten_text(product, 24),
      TRUE ~ coalesce(locus_tag, "")
    ),
    # Prefer symbol for highlight; `gene` in TSV is usually locus_tag but may be a symbol.
    highlight =
      (str_to_lower(gene_symbol) %in% str_to_lower(highlight_gene_symbols)) |
        (str_to_lower(coalesce(gene, "")) %in% str_to_lower(highlight_gene_symbols)) |
        (tier == "Tier1"),
    hl = if_else(highlight, "highlight", "background")
  )

if (nrow(variants_chr) == 0) {
  stop("No variants found for chrom == 'contig_2'.")
}

cat("\nvariants_chr (chrom, pos, gene, product, tier, category, gene_label):\n")
print(
  variants_chr %>%
    select(chrom, pos, gene, product, tier, category, gene_label)
)

cat("\nSame rows + resolved symbol, plot category, highlight:\n")
print(
  variants_chr %>%
    select(
      chrom, pos, gene_symbol, plot_category,
      highlight, tier, gene_label
    )
)

label_genes <- c("agrA", "cdaR", "mdrP", "yclM", "fumC", "ebh")
variants_labeled <- variants_chr %>%
  filter(highlight, nzchar(gene_label))
mini_label_genes <- c("agrA", "yclM", "fumC", "mdrP", "ebh")
variants_labeled_mini <- variants_chr %>%
  filter(str_to_lower(gene_symbol) %in% str_to_lower(mini_label_genes), nzchar(gene_label)) %>%
  mutate(mini_label = mini_protein_label(gene_symbol, effect))

genome_len <- max(variants_chr$pos, na.rm = TRUE) * 1.02
highlight_blocks <- variants_chr %>%
  filter(highlight) %>%
  transmute(
    gene_symbol = gene_symbol,
    locus_tag = locus_tag,
    gene_label = gene_label
  ) %>%
  distinct()

gene_blocks <- gene_features %>%
  filter(seqid == "contig_2") %>%
  inner_join(
    highlight_blocks,
    by = c("locus_tag", "gene_symbol")
  ) %>%
  transmute(
    xmin = pmin(start, end),
    xmax = pmax(start, end),
    gene_label = gene_label
  ) %>%
  distinct()

# Make region/gene blocks visible at whole-genome scale.
min_block_bp <- max(20000, genome_len * 0.015)
gene_blocks <- gene_blocks %>%
  mutate(
    mid = (xmin + xmax) / 2,
    xmin_plot = pmax(0, mid - min_block_bp / 2),
    xmax_plot = pmin(genome_len, mid + min_block_bp / 2)
  )

# Functional region blocks (group-level) so category regions are easy to see.
region_blocks <- variants_chr %>%
  filter(highlight) %>%
  group_by(pheno_cluster) %>%
  summarise(
    xmin = pmax(0, min(pos, na.rm = TRUE) - min_block_bp * 0.8),
    xmax = pmin(genome_len, max(pos, na.rm = TRUE) + min_block_bp * 0.8),
    .groups = "drop"
  ) %>%
  mutate(
    region_label = as.character(pheno_cluster),
    region_short = case_when(
      region_label == "Quorum/Virulence regulation" ~ "Putative WTA access regulatory axis",
      region_label == "Ion homeostasis/Membrane transport" ~ "Ion membrane homeostasis axis",
      region_label == "Cell-wall precursor/envelope metabolism" ~ "Envelope metabolism axis",
      region_label == "Surface architecture/access" ~ "Surface accessibility axis",
      region_label == "rRNA-adjacent/structural context" ~ "rRNA-adjacent",
      TRUE ~ "Other background"
    )
  ) %>%
  arrange(xmin) %>%
  mutate(
    label_row = row_number() %% 4,
    track_row = row_number() %% 4,
    track_mid = case_when(
      track_row == 0 ~ -0.145,
      track_row == 1 ~ -0.120,
      track_row == 2 ~ -0.095,
      TRUE ~ -0.070
    ),
    track_h = 0.010,
    track_ymin = track_mid - track_h,
    track_ymax = track_mid + track_h,
    label_y = case_when(
      label_row == 0 ~ -0.170,
      label_row == 1 ~ -0.153,
      label_row == 2 ~ -0.136,
      TRUE ~ -0.119
    )
  )

shape_map <- c(Tier1 = 17, Tier2 = 16, Tier3 = 15)

cat_palette <- c(
  "regulator" = "#E69F00",
  "envelope/transport" = "#0072B2",
  "metabolism" = "#009E73",
  "rRNA_adjacent" = "#CC79A7",
  "other" = "#999999"
)

cluster_palette <- c(
  "Surface accessibility axis" = "#bdbdbd",
  "Envelope metabolism axis" = "#6cc4b5",
  "Putative WTA access regulatory axis" = "#f3c86b",
  "Ion membrane homeostasis axis" = "#87b7e3",
  "rRNA-adjacent" = "#caa3d8",
  "Other background" = "#d9d9d9"
)

p_lollipop <- ggplot() +
  annotate(
    "rect",
    xmin = 0, xmax = genome_len,
    ymin = -0.12, ymax = -0.08,
    fill = "#b0b0ac", color = NA
  ) +
  geom_rect(
    data = region_blocks,
    aes(xmin = xmin, xmax = xmax, ymin = -0.22, ymax = -0.19, fill = region_short),
    alpha = 0.65, color = NA
  ) +
  geom_text(
    data = region_blocks,
    aes(x = (xmin + xmax) / 2, y = -0.205, label = region_label),
    size = 2.7, color = "#2f2f2f"
  ) +
  geom_rect(
    data = gene_blocks,
    aes(xmin = xmin_plot, xmax = xmax_plot, ymin = -0.18, ymax = -0.03),
    fill = "#6b8ec1", color = "white", linewidth = 0.3
  ) +
  geom_text(
    data = gene_blocks,
    aes(x = (xmin_plot + xmax_plot) / 2, y = -0.105, label = gene_label),
    size = 3.1, color = "white"
  ) +
  geom_segment(
    data = variants_chr,
    aes(
      x = pos, xend = pos, y = -0.03, yend = 0.25,
      color = plot_category,
      alpha = hl,
      linewidth = hl
    )
  ) +
  geom_point(
    data = variants_chr,
    aes(
      y = 0.25,
      x = pos,
      color = plot_category,
      shape = tier,
      alpha = hl,
      size = hl
    )
  ) +
  geom_text_repel(
    data = variants_labeled,
    aes(x = pos, y = 0.25, label = gene_label, color = plot_category),
    size = 3.2,
    alpha = 0.95,
    max.overlaps = 100,
    box.padding = 0.35,
    point.padding = 0.25,
    segment.alpha = 0.5,
    show.legend = FALSE
  ) +
  scale_linewidth_manual(
    values = c(background = 0.35, highlight = 1.15),
    guide = "none"
  ) +
  scale_alpha_manual(
    values = c(background = 0.3, highlight = 1),
    guide = "none"
  ) +
  scale_size_manual(
    values = c(background = 2.3, highlight = 5),
    guide = "none"
  ) +
  scale_shape_manual(values = shape_map, drop = FALSE) +
  scale_color_manual(values = cat_palette, drop = FALSE) +
  scale_fill_manual(values = cluster_palette, drop = FALSE, guide = "none") +
  scale_x_continuous(
    limits = c(0, genome_len),
    expand = c(0, 0)
  ) +
  coord_cartesian(ylim = c(-0.24, 0.34), clip = "off") +
  labs(
    title = "A. High-confidence variants along contig_2",
    x = "Genomic position (bp)",
    y = NULL,
    color = "Category",
    shape = "Tier"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

p_lollipop_mini <- ggplot() +
  annotate(
    "rect",
    xmin = 0, xmax = genome_len,
    ymin = -0.02, ymax = 0.02,
    fill = "#b7b9b4", color = NA
  ) +
  geom_rect(
    data = region_blocks,
    aes(xmin = xmin, xmax = xmax, ymin = track_ymin, ymax = track_ymax, fill = region_short),
    alpha = 0.70, color = NA
  ) +
  geom_text_repel(
    data = region_blocks,
    aes(x = (xmin + xmax) / 2, y = label_y, label = region_short),
    color = "#4a4a4a",
    size = 2.5,
    direction = "x",
    min.segment.length = 0,
    segment.alpha = 0.35,
    box.padding = 0.18,
    point.padding = 0.1,
    max.overlaps = 100
  ) +
  geom_rect(
    data = gene_blocks,
    aes(xmin = xmin_plot, xmax = xmax_plot, ymin = -0.06, ymax = 0.06),
    fill = "#9fb5d9", color = "white", linewidth = 0.2, alpha = 0.9
  ) +
  geom_segment(
    data = variants_chr,
    aes(
      x = pos, xend = pos, y = 0.02,
      yend = if_else(highlight, 0.17, 0.12),
      alpha = hl
    ),
    color = "#497d3c",
    linewidth = 0.4
  ) +
  geom_point(
    data = variants_chr,
    aes(x = pos, y = if_else(highlight, 0.17, 0.12), alpha = hl),
    color = "#497d3c",
    size = 2
  ) +
  geom_text_repel(
    data = variants_labeled_mini,
    aes(x = pos, y = 0.17, label = mini_label),
    color = "#1f3b16",
    size = 2.6,
    max.overlaps = 100,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.alpha = 0.4
  ) +
  scale_alpha_manual(
    values = c(background = 0.25, highlight = 1),
    guide = "none"
  ) +
  scale_fill_manual(values = cluster_palette, drop = FALSE, guide = "none") +
  scale_x_continuous(
    limits = c(0, genome_len),
    expand = c(0, 0)
  ) +
  coord_cartesian(ylim = c(-0.19, 0.24), clip = "off") +
  labs(x = "Genomic position (bp)", y = NULL) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_blank()
  )

cat_counts <- variants_chr %>%
  count(plot_category, name = "n") %>%
  arrange(desc(n))

p_bar <- ggplot(cat_counts, aes(x = plot_category, y = n, fill = plot_category)) +
  geom_col(width = 0.8, show.legend = FALSE) +
  geom_text(aes(label = n), vjust = -0.35, size = 3.4) +
  scale_fill_manual(values = cat_palette, drop = FALSE) +
  expand_limits(y = max(cat_counts$n) * 1.12) +
  labs(
    title = "B. Distribution of variants by functional category",
    x = "Category",
    y = "Variant count"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

candidate_genes <- label_genes

candidates_tbl <- variants_chr %>%
  mutate(gene_match = str_to_lower(gene_symbol)) %>%
  filter(tier == "Tier1" | gene_match %in% str_to_lower(candidate_genes)) %>%
  select(any_of(c(
    "chrom", "pos", "ref", "alt", "gene_symbol", "locus_tag", "product",
    "effect", "impact", "tier", "plot_category"
  ))) %>%
  arrange(chrom, pos)

write_tsv(candidates_tbl, out_candidates)

combined_plot <- wrap_plots(
  list(p_lollipop, p_lollipop_mini, p_bar),
  ncol = 1,
  heights = c(2.2, 1, 1.2)
)
ggsave(out_png, combined_plot, width = 7, height = 7, dpi = 300)
ggsave(out_jpg, combined_plot, width = 7, height = 7, dpi = 300)
ggsave(out_panelA_main_png, p_lollipop, width = 7, height = 3.8, dpi = 300)
ggsave(out_panelA_main_jpg, p_lollipop, width = 7, height = 3.8, dpi = 300)
ggsave(out_panelA_mini_png, p_lollipop_mini, width = 7, height = 3.1, dpi = 300)
ggsave(out_panelA_mini_jpg, p_lollipop_mini, width = 7, height = 3.1, dpi = 300)
ggsave(out_panelB_png, p_bar, width = 7, height = 3, dpi = 300)
ggsave(out_panelB_jpg, p_bar, width = 7, height = 3, dpi = 300)

info <- file.info(out_png)
cat("\nSaved figure:", out_png, "\n")
cat("File size (bytes):", info$size, "\n")
cat("Saved figure:", out_jpg, " size:", file.info(out_jpg)$size, "\n")
cat("Saved panel:", out_panelA_main_png, "\n")
cat("Saved panel:", out_panelA_mini_png, "\n")
cat("Saved panel:", out_panelB_png, "\n")
cat("Saved candidates table:", out_candidates, "\n")
