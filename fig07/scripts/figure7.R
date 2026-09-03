# Figure 7 — host range of OCPT-Ludmila across 83 Staphylococcus aureus strains:
# core-genome phylogeny, host group, and the candidate receptor gene grid.
#
# Run with ../run.sh, which supplies FIG_DATA, FIG_WORK and FIG_OUT.
#
# PROVENANCE
#
# The tree and its typography are taken directly from the original figure,
# Host-Range/scripts/06_structure/R/plot_figure5a_phylogeny_grid.R: midpoint
# rooting, ladderize(right = FALSE), branch lengths rescaled to a fixed depth,
# aligned tip labels, theme_tree2() with a substitutions-per-site axis, and the
# three-column layout at widths 8 : 1.15 : 3. The host-group column is the same
# lettered tile, and the gene grid is the same eight genes in the same order.
#
# WHAT DIFFERS, AND WHY
#
# Only the gene calls. The original read merged_features.tsv, which holds raw
# Panaroo cluster calls; those are wrong for several genes in this panel, most
# visibly tagO and tarF at 3/83 when both are in fact universal. The same
# splitting is what manufactured the sraP association. So the grid is drawn from
# the corrected gene families (scripts/gene_grid.py).
#
# The per-gene statistics are deliberately NOT printed on the figure. No gene in
# this panel survives correction, and a nominal p beside a gene name invites it
# to be read as a result. tarM is the case in point: nominal p = 0.038, but
# permutation family-wise p = 1.000, p = 0.167 with the VISA panel removed, and
# p = 0.608 collapsed to independent lineages. Those numbers belong in
# work/gene_grid_stats.tsv and the supplementary table, with their corrections
# attached.
#
# THE RESULT THE FIGURE REPORTS
#
# None of the genes shown is associated with susceptibility after correction.
# Genome-wide, three families do clear both Benjamini-Hochberg and a permutation
# family-wise threshold, but scripts/clonal_confounding.py shows they are one
# co-inherited locus present in 13 of 14 VISA isolates, that the VISA panel is
# itself markedly phage-resistant (3/14 vs 50/69 susceptible, Fisher p = 5e-04),
# and that those 14 isolates are only 2 lineages at 50 SNPs. Testing once per
# independent lineage removes the association (p = 1.00, 0.62, 0.58).
#
# Inputs:
#   data/coregenome_tree.newick        IQ-TREE ML tree, 83 tips
#   data/phylogeny_tip_metadata.tsv    strain_id, host_group, ludmila_susceptible
#   work/gene_grid.tsv                 corrected presence calls
#   work/gene_grid_stats.tsv           per-gene beta, p, q
#
# Optional overrides:
#   FIG_SUBTITLE=0   suppress the methods subtitle
#   FIG_DPI          PNG resolution (default 600)
#   FIG_GENES        comma-separated gene order for the grid

if (!interactive()) pdf(NULL)

suppressPackageStartupMessages({
  library(ggplot2)
  library(ape)
  library(phytools)
  library(ggtree)
  library(patchwork)
})

need <- function(v) {
  x <- Sys.getenv(v, unset = NA_character_)
  if (is.na(x)) stop(v, " must be set. Run this figure with ./run.sh")
  x
}
env_or <- function(v, d) { x <- Sys.getenv(v, unset = ""); if (nzchar(x)) x else d }

data_dir    <- need("FIG_DATA")
results_dir <- need("FIG_OUT")
work_dir    <- env_or("FIG_WORK", file.path(dirname(data_dir), "work"))
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

FIG_DPI      <- as.numeric(env_or("FIG_DPI", "600"))
FIG_SUBTITLE <- env_or("FIG_SUBTITLE", "1") != "0"

# ---- constants, verbatim from the original ---------------------------------

PHYLO_X_MAX             <- 0.035
TXT_COL                 <- "black"
SIZE_TIPLAB             <- 3.8
SIZE_AXIS_X             <- 11
SIZE_AXIS_TITLE         <- 12
SIZE_MAIN_TITLE         <- 14
SIZE_PHENO_MARK         <- 3.2
SIZE_GENE_GRID_HDR      <- 20
GENE_GRID_LABEL_ANGLE   <- 90
GENE_GRID_MARGIN_BOTTOM <- 36
GENE_PRESENT_GREY       <- "#3a3a3a"
GENE_ABSENT_GREY        <- "#d8d8d8"

GROUP_COLORS <- c(A = "#2F6F3E", B = "#8B1E24", C = "#8B1E24",
                  D = "#4A5D4A", E = "#6B4E3D")
GROUP_MARK_COL <- c(A = "#B8E8C8", B = "#F5C4C8", C = "#F5C4C8",
                    D = "#C8DCC8", E = "#EDD8C8")

# the eight genes of the original grid, in the original left-to-right order
GENE_ORDER <- strsplit(env_or("FIG_GENES", "tarL,tarM,tarS,sdrC,sraP,sasG,tagX,tagG"), ",")[[1]]

# sdrC and sdrD merged into one family during re-calling; gene_grid.py emits it
# under the name sdrD. Display it under the original label with both names.
ALIAS <- c(sdrC = "sdrD")
ALIAS_LABEL <- c(sdrC = "sdrC/sdrD")

base_font <- function() {
  bf <- "Times New Roman"
  if (requireNamespace("systemfonts", quietly = TRUE) &&
      !bf %in% unique(systemfonts::system_fonts()$family)) bf <- "Times"
  bf
}
bf <- base_font()

scale_tree_to_depth <- function(tree, target_max) {
  H <- phytools::nodeHeights(tree)
  mx <- max(H, na.rm = TRUE)
  if (!is.finite(mx) || mx <= 0) return(tree)
  tree$edge.length <- tree$edge.length * (target_max / mx)
  tree
}

# ---- inputs ----------------------------------------------------------------

tree_path <- file.path(data_dir, "coregenome_tree.newick")
meta_path <- file.path(data_dir, "phylogeny_tip_metadata.tsv")
grid_path <- file.path(work_dir, "gene_grid.tsv")
stat_path <- file.path(work_dir, "gene_grid_stats.tsv")
for (p in c(tree_path, meta_path, grid_path, stat_path)) {
  if (!file.exists(p)) stop("Missing input: ", p,
                            "\nRun scripts/gene_grid.py first (./run.sh analysis).")
}

tree <- ape::read.tree(tree_path)
tree <- phytools::midpoint.root(tree)
tree <- ape::ladderize(tree, right = FALSE)
tree <- scale_tree_to_depth(tree, PHYLO_X_MAX)

meta  <- utils::read.delim(meta_path, stringsAsFactors = FALSE, check.names = FALSE)
grid  <- utils::read.delim(grid_path, stringsAsFactors = FALSE)
gstat <- utils::read.delim(stat_path, stringsAsFactors = FALSE)

meta <- meta[!duplicated(meta$strain_id), , drop = FALSE]
if (!setequal(tree$tip.label, meta$strain_id)) {
  stop("Tree tips and metadata strain_id sets differ")
}

resolve <- function(g) if (g %in% ALIAS) g else if (g %in% names(ALIAS)) ALIAS[[g]] else g
wanted  <- vapply(GENE_ORDER, resolve, character(1))
absent  <- GENE_ORDER[!wanted %in% grid$gene]
if (length(absent)) {
  stop("Genes requested for the grid but not present in gene_grid.tsv: ",
       paste(absent, collapse = ", "),
       "\nThey were dropped by gene_grid.py as invariant or unmappable; ",
       "set FIG_GENES to a different panel or re-run gene_grid.py.")
}

# ---- tip order, shared by all three panels ---------------------------------

meta_tip <- meta; names(meta_tip)[names(meta_tip) == "strain_id"] <- "label"
p_ord <- ggtree(tree, ladderize = FALSE) %<+% meta_tip
td <- p_ord$data
tip_df <- td[as.logical(td$isTip), , drop = FALSE]
# ASCENDING, deliberately. ggtree draws the top tip at the HIGHEST y, but a
# ggplot discrete scale draws factor level 1 at the BOTTOM. Ordering by
# decreasing y (as the original script did) therefore puts the top tip at level
# 1 and renders the host-group column and the gene grid mirrored against the
# tree — silently, because every panel is individually well formed. Verified by
# decoding the rendered PNG: ascending scores 83/83 rows, decreasing scores 9/83.
strain_levels <- tip_df$label[order(tip_df$y)]

meta_ord <- meta[match(strain_levels, meta$strain_id), , drop = FALSE]
meta_ord$strain <- factor(meta_ord$strain_id, levels = strain_levels)

# ---- panel 1: tree ---------------------------------------------------------

p_tree <- ggtree(tree, ladderize = FALSE) %<+% meta_tip +
  geom_tree(linewidth = 0.35, color = "gray15") +
  geom_tiplab(align = TRUE, size = SIZE_TIPLAB, offset = 0.0012,
              family = bf, colour = TXT_COL) +
  scale_x_continuous(
    breaks = c(0, 0.007, 0.014, 0.021, 0.028, PHYLO_X_MAX),
    labels = function(x) sprintf("%.3f", x),
    expand = expansion(mult = c(0.02, 0.14))
  ) +
  labs(x = "Substitutions per site (branch length)", y = NULL) +
  theme_tree2() +
  theme(
    text         = element_text(family = bf, colour = TXT_COL),
    plot.title   = element_text(size = SIZE_MAIN_TITLE, colour = TXT_COL),
    axis.text.x  = element_text(size = SIZE_AXIS_X, colour = TXT_COL),
    axis.title.x = element_text(size = SIZE_AXIS_TITLE, colour = TXT_COL),
    axis.ticks   = element_line(colour = TXT_COL),
    axis.line.x  = element_line(colour = TXT_COL),
    legend.position = "none"
  )

# ---- panel 2: host group ---------------------------------------------------

hg_chr <- as.character(meta_ord$host_group)
col_map <- GROUP_COLORS; mark_map <- GROUP_MARK_COL
for (u in setdiff(unique(hg_chr), names(col_map))) {
  col_map[[u]] <- "#888888"; mark_map[[u]] <- "#E8E8E8"
}
pheno_df <- data.frame(
  strain = meta_ord$strain, x = "Host group",
  fill = unname(col_map[hg_chr]), mark = hg_chr,
  markcol = unname(mark_map[hg_chr]), stringsAsFactors = FALSE
)

p_pheno <- ggplot(pheno_df, aes(x = x, y = strain)) +
  geom_tile(aes(fill = fill), color = "white", linewidth = 0.12,
            width = 1, height = 1) +
  geom_text(aes(label = mark, colour = markcol), family = bf,
            size = SIZE_PHENO_MARK, fontface = "bold", vjust = 0.48,
            show.legend = FALSE) +
  scale_colour_identity() + scale_fill_identity() +
  labs(x = NULL, y = NULL) +
  theme_void() +
  theme(
    text = element_text(family = bf, colour = TXT_COL),
    axis.text.x = element_text(angle = GENE_GRID_LABEL_ANGLE, hjust = 1, vjust = 0.5,
                               size = SIZE_GENE_GRID_HDR, colour = TXT_COL, family = bf),
    plot.margin = margin(4, 2, GENE_GRID_MARGIN_BOTTOM, 2)
  )

# ---- panel 3: gene grid ----------------------------------------------------

glong <- grid[grid$gene %in% wanted, , drop = FALSE]
glong$strain <- factor(glong$strain, levels = strain_levels)
glong$fill <- ifelse(glong$present == 1L, GENE_PRESENT_GREY, GENE_ABSENT_GREY)
glong$gene <- factor(glong$gene, levels = wanted)

# Gene names only. The per-gene statistics are NOT on the figure: no gene in
# this panel survives correction, and printing a nominal p next to a gene name
# invites it to be read as a result. tarM is the case that matters -- nominal
# p = 0.038, but permutation family-wise p = 1.000, p = 0.167 with the VISA
# panel removed, and p = 0.608 collapsed to independent lineages. The numbers
# live in work/gene_grid_stats.tsv and the supplementary table, where they can
# be read with their corrections attached.
gene_label <- function(fam) {
  disp <- names(wanted)[match(fam, wanted)]
  if (is.na(disp)) disp <- fam
  if (disp %in% names(ALIAS_LABEL)) disp <- ALIAS_LABEL[[disp]]
  disp
}
lab_map <- vapply(wanted, gene_label, character(1))

p_genes <- ggplot(glong, aes(x = gene, y = strain)) +
  geom_tile(aes(fill = fill), color = "white", linewidth = 0.12,
            width = 1, height = 1) +
  scale_fill_identity() +
  scale_x_discrete(labels = function(x) lab_map[x]) +
  labs(x = NULL, y = NULL) +
  theme_void() +
  theme(
    text = element_text(family = bf, colour = TXT_COL),
    axis.text.x = element_text(angle = GENE_GRID_LABEL_ANGLE, hjust = 1, vjust = 0.5,
                               size = SIZE_GENE_GRID_HDR, colour = TXT_COL, family = bf),
    plot.margin = margin(4, 4, GENE_GRID_MARGIN_BOTTOM, 2)
  )

# ---- assemble, at the original proportions ---------------------------------

n_sig <- sum(suppressWarnings(as.numeric(gstat$q_bh)) < 0.05, na.rm = TRUE)
sub <- sprintf(
  paste0("%d isolates, %d susceptible. Gene presence re-derived from the pangenome. ",
         "No gene shown is associated with susceptibility after correction for ",
         "multiple testing and population structure (see Supplementary Table)."),
  nrow(meta), sum(meta$ludmila_susceptible == 1))

combo <- p_tree + p_pheno + p_genes +
  plot_layout(widths = c(8, 1.15, 3)) +
  plot_annotation(
    subtitle = if (FIG_SUBTITLE) sub else NULL,
    theme = theme(text = element_text(family = bf, colour = TXT_COL),
                  plot.subtitle = element_text(size = 11, hjust = 0))
  )

n <- length(strain_levels)
fig_h <- max(8, 0.24 * n + 2)
fig_w <- 18

base <- file.path(results_dir, "figure7_host_range")
ggsave(paste0(base, ".png"), combo, width = fig_w, height = fig_h,
       units = "in", dpi = FIG_DPI, bg = "white", limitsize = FALSE)
t0 <- Sys.time(); pdf_ok <- FALSE
suppressWarnings(try({
  ggsave(paste0(base, ".pdf"), combo, width = fig_w, height = fig_h,
         units = "in", device = grDevices::cairo_pdf, bg = "white", limitsize = FALSE)
  fi <- file.info(paste0(base, ".pdf"))
  pdf_ok <- !is.na(fi$size) && fi$size > 0 && fi$mtime >= t0
}, silent = TRUE))
message("Wrote: figure7_host_range.png (", fig_w, " x ", round(fig_h, 1), " in)",
        if (pdf_ok) " and .pdf" else "  (no .pdf: cairo_pdf wrote no file)")
message("Genes shown: ", paste(lab_map, collapse = ", "))
