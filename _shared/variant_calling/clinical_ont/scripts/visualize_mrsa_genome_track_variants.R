#!/usr/bin/env Rscript

# Linear genome-track view of high-confidence MRSA variants.
# Outputs:
# - figures/Figure_genome_track_variants.png

needed <- c("readr", "dplyr", "ggplot2", "stringr", "ggrepel")
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
  library(ggrepel)
})

input_tsv <- "tables/variants_for_paper.tsv"
gff_path <- "annotation/snpeff/data/mrsa_ref_pre/genes.gff"
out_png <- "figures/Figure_genome_track_variants.png"
dir.create("figures", recursive = TRUE, showWarnings = FALSE)

key_genes <- c("agrA", "yclM", "fumC", "mdrP", "ebh")
key_gene_labels <- c("agrA", "yclM", "fumC-upstream", "mdrP", "ebh")

read_gene_features <- function(path) {
  if (!file.exists(path)) {
    return(tibble())
  }
  read_tsv(
    path,
    col_names = c(
      "seqid", "source", "type", "start", "end",
      "score", "strand", "phase", "attrs"
    ),
    comment = "#",
    show_col_types = FALSE,
    col_types = cols(.default = col_character())
  ) %>%
    transmute(
      chrom = seqid,
      type = type,
      start = as.numeric(start),
      end = as.numeric(end),
      locus_tag = str_match(attrs, "locus_tag=([^;]+)")[, 2],
      gene_symbol = coalesce(
        str_match(attrs, "gene=([^;]+)")[, 2],
        str_match(attrs, "Name=([^;]+)")[, 2],
        str_match(attrs, "locus_tag=([^;]+)")[, 2]
      ),
      product_annot = str_match(attrs, "product=([^;]+)")[, 2],
      attrs = attrs
    )
}

derive_category <- function(gene_sym, product_txt, effect_txt) {
  g <- str_to_lower(coalesce(gene_sym, ""))
  p <- str_to_lower(coalesce(product_txt, ""))
  e <- str_to_lower(coalesce(effect_txt, ""))
  case_when(
    g %in% c("agra", "cdar", "agrb", "agrc", "agrd") ~ "regulator",
    g %in% c("yclm", "fumc") ~ "metabolism",
    g == "mdrp" | str_detect(p, "antiporter|transporter|efflux|symporter|channel|porin") ~ "envelope/transport",
    g == "ebh" ~ "other",
    str_detect(e, "rrna") | str_detect(p, "rrna") ~ "rRNA_adjacent",
    TRUE ~ "other"
  )
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

variants <- read_tsv(input_tsv, show_col_types = FALSE) %>%
  filter(chrom == "contig_2") %>%
  mutate(
    pos = as.numeric(pos),
    locus_tag = coalesce(locus_tag, gene),
    tier = str_replace_all(as.character(tier), "\\s+", "")
  )

if (nrow(variants) == 0) stop("No contig_2 variants found.")

gff <- read_gene_features(gff_path)
locus_map <- gff %>%
  filter(!is.na(locus_tag)) %>%
  distinct(locus_tag, .keep_all = TRUE) %>%
  select(locus_tag, gene_symbol, product_annot)

variants_chr <- variants %>%
  left_join(locus_map, by = "locus_tag") %>%
  mutate(
    gene_plot = coalesce(gene_symbol, gene),
    product_plot = coalesce(product, product_annot, ""),
    category_plot = derive_category(gene_plot, product_plot, effect),
    effect_plot = effect_class(effect)
  )

# Supplement key loci from mutations.tsv if they are absent from variants_for_paper.tsv.
mut_path <- "tables/mutations.tsv"
if (file.exists(mut_path)) {
  key_l <- str_to_lower(key_genes)
  mut <- read_tsv(mut_path, show_col_types = FALSE) %>%
    filter(chrom == "contig_2", unique_to_post %in% TRUE) %>%
    left_join(locus_map, by = c("gene" = "locus_tag")) %>%
    mutate(gene_sym = str_to_lower(coalesce(gene_symbol, gene))) %>%
    filter(gene_sym %in% key_l) %>%
    transmute(
      chrom = chrom,
      pos = as.numeric(pos),
      ref = ref,
      alt = alt,
      gene = gene,
      locus_tag = gene,
      product = coalesce(product, ""),
      effect = coalesce(effect_detail, effect),
      impact = case_when(
        str_detect(str_to_lower(coalesce(effect_detail, effect, "")), "missense") ~ "MODERATE",
        str_detect(str_to_lower(coalesce(effect_detail, effect, "")), "synonymous") ~ "LOW",
        TRUE ~ "MODIFIER"
      ),
      category = "other",
      tier = "Tier2",
      focus_set = "yes",
      caller_support = "mutations_supplement",
      pass_both_callers = NA,
      medaka_qual = as.numeric(coalesce(qual, NA_real_)),
      medaka_dp = NA_real_,
      medaka_af = NA_real_,
      clair3_qual = NA_real_,
      clair3_dp = NA_real_,
      clair3_af = NA_real_,
      region_flag = "other",
      gene_symbol = coalesce(gene_symbol, gene),
      product_annot = NA_character_
    ) %>%
    mutate(
      gene_plot = coalesce(gene_symbol, gene),
      product_plot = coalesce(product, product_annot, ""),
      category_plot = derive_category(gene_plot, product_plot, effect),
      effect_plot = effect_class(effect)
    )
  variants_chr <- bind_rows(variants_chr, mut) %>%
    distinct(chrom, pos, ref, alt, .keep_all = TRUE)
}

genome_max <- max(variants_chr$pos, na.rm = TRUE)
window_bp <- 500

# Key gene bands from GFF when available; fallback to +/-500 around variant position.
key_bands_gff <- gff %>%
  filter(chrom == "contig_2", str_to_lower(gene_symbol) %in% str_to_lower(key_genes)) %>%
  mutate(
    gene = gene_symbol,
    start = pmin(start, end),
    end = pmax(start, end),
    category = derive_category(gene_symbol, product_annot, "")
  ) %>%
  select(gene, start, end, category) %>%
  distinct(gene, .keep_all = TRUE)

key_bands_fallback <- variants_chr %>%
  filter(str_to_lower(gene_plot) %in% str_to_lower(key_genes)) %>%
  transmute(
    gene = gene_plot,
    start = pmax(0, pos - window_bp),
    end = pmin(genome_max, pos + window_bp),
    category = category_plot
  ) %>%
  distinct(gene, .keep_all = TRUE)

gene_bands <- bind_rows(key_bands_gff, key_bands_fallback) %>%
  arrange(gene) %>%
  distinct(gene, .keep_all = TRUE)

# Add explicit rRNA loci bands if available.
rna_bands <- gff %>%
  filter(chrom == "contig_2") %>%
  filter(type == "rRNA") %>%
  transmute(
    gene = "16S/5S rRNA",
    start = pmin(start, end),
    end = pmax(start, end),
    category = "rRNA_adjacent"
  ) %>%
  distinct(start, end, .keep_all = TRUE)

if (nrow(rna_bands) > 0) {
  rna_merged <- rna_bands %>%
    arrange(start) %>%
    mutate(
      prev_end = lag(end, default = first(start)),
      new_grp = if_else(start - prev_end > 20000, 1L, 0L),
      grp = cumsum(new_grp)
    ) %>%
    group_by(grp) %>%
    summarise(
      gene = "16S/5S rRNA",
      start = min(start, na.rm = TRUE),
      end = max(end, na.rm = TRUE),
      category = "rRNA_adjacent",
      .groups = "drop"
    ) %>%
    select(gene, start, end, category)
  gene_bands <- bind_rows(gene_bands, rna_merged)
}

gene_bands <- gene_bands %>%
  mutate(
    y0 = -0.12,
    y1 = 0.12,
    mid = (start + end) / 2
  )

marker_df <- variants_chr %>%
  mutate(
    marker_y = 0.42,
    marker_label = case_when(
      str_to_lower(gene_plot) == "agra" ~ "agrA",
      str_to_lower(gene_plot) == "yclm" ~ "yclM",
      str_to_lower(gene_plot) == "fumc" & str_detect(str_to_lower(effect), "upstream") ~ "fumC-upstream",
      str_to_lower(gene_plot) == "mdrp" ~ "mdrP",
      str_to_lower(gene_plot) == "ebh" ~ "ebh",
      TRUE ~ ""
    )
  )

cat("Gene bands preview:\n")
print(head(gene_bands, 10))
cat("\nVariant markers preview:\n")
print(head(marker_df %>% select(chrom, pos, gene_plot, effect_plot, category_plot, tier, marker_label), 12))

cat_palette <- c(
  "regulator" = "#E69F00",
  "envelope/transport" = "#0072B2",
  "metabolism" = "#009E73",
  "rRNA_adjacent" = "#CC79A7",
  "other" = "#999999"
)
shape_map <- c("missense" = 17, "upstream" = 15, "synonymous" = 16, "other" = 18)

p <- ggplot() +
  annotate("rect", xmin = 0, xmax = genome_max * 1.02, ymin = -0.03, ymax = 0.03, fill = "#9fa19c", color = NA) +
  geom_rect(
    data = gene_bands,
    aes(xmin = start, xmax = end, ymin = y0, ymax = y1, fill = category),
    alpha = 0.75, color = NA
  ) +
  geom_text(
    data = gene_bands %>% filter(gene == "16S/5S rRNA"),
    aes(x = mid, y = -0.18, label = gene),
    size = 2.7,
    color = "#5c2c4f",
    alpha = 0.9
  ) +
  geom_text_repel(
    data = gene_bands %>% filter(gene != "16S/5S rRNA"),
    aes(x = mid, y = -0.18, label = gene, color = category),
    size = 3,
    direction = "x",
    min.segment.length = 0,
    segment.alpha = 0.4,
    max.overlaps = 100,
    box.padding = 0.18,
    point.padding = 0.1,
    show.legend = FALSE
  ) +
  geom_point(
    data = marker_df,
    aes(x = pos, y = marker_y, color = category_plot, shape = effect_plot),
    size = 3.2
  ) +
  geom_text(
    data = marker_df %>% filter(marker_label != ""),
    aes(x = pos, y = marker_y + 0.08, label = marker_label, color = category_plot),
    size = 3.1
  ) +
  scale_color_manual(values = cat_palette, drop = FALSE, name = "Variant category") +
  scale_fill_manual(values = cat_palette, drop = FALSE, name = "Gene band category") +
  scale_shape_manual(values = shape_map, name = "Effect") +
  scale_x_continuous(
    limits = c(0, genome_max * 1.02),
    breaks = seq(0, floor(genome_max / 5e5) * 5e5 + 5e5, by = 5e5),
    labels = function(x) sprintf("%.1f Mb", x / 1e6),
    expand = c(0, 0)
  ) +
  coord_cartesian(ylim = c(-0.27, 0.60), clip = "off") +
  labs(
    title = "Linear genome track of high-confidence variants on contig_2",
    x = "Genomic position",
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

ggsave(out_png, p, width = 7.2, height = 3.8, dpi = 300)
cat("\nSaved:", out_png, "\n")
