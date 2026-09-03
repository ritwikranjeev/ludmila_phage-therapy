#!/usr/bin/env Rscript
# =============================================================================
# ViPTree rectangular tree: S. aureus phages, OCPT-Ludmila highlighted
# Reads ViPTree/ViPTreeGen Newick + genome table, plots rectangular cladogram
# with OCPT-Ludmila and its clade emphasized. Exports PDF and PNG.
#
# Usage: run from LudmilaCharacterization project root. Data is read from
# data/.../viptree/; outputs (PDF/PNG) go to analyses/.../2b_Ludmila_Phylogeny/results/.
#   cd LudmilaCharacterization && Rscript analyses/Figure2_Genomic_Characterization/2b_Ludmila_Phylogeny/scripts/viptree_OCPT_Ludmila_rectangular.R
# =============================================================================

# --- Libraries ---------------------------------------------------------------
library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
library(tidyr)

# Ignore negative edge lengths in tree (ViPTree can produce them)
options(ignore.negative.edge = TRUE)

# --- Parameters (change these as needed) --------------------------------------
# Run from LudmilaCharacterization project root so these relative paths work.
# Path to the folder containing ViPTree outputs (data only; no results written here)
viptree_dir    <- "data/Figure2_Genomic_Characterization/2b_Ludmila_Phylogeny/viptree"
# Folder for all analysis outputs (PDF/PNG written here)
results_dir    <- "analyses/Figure2_Genomic_Characterization/2b_Ludmila_Phylogeny/results"

# Primary Newick tree file (full proteomic tree)
tree_file      <- "set.bionj.newick"
# ViPTree genome table (Figure 2b–style metadata: ID, host, etc.)
genome_table_file <- "genome_table/genome.tsv"

# Host string to define S. aureus phages (only these tips are kept in the tree)
host_filter    <- "Staphylococcus aureus"
# Tip label (or substring) for the focal phage to highlight and box its clade
ocpt_label     <- "OCPT-Ludmila"

# Output filenames (written into results_dir); both rectangular and circular
out_pdf_rect   <- "viptree_OCPT_Ludmila_clade_rectangular.pdf"
out_png_rect   <- "viptree_OCPT_Ludmila_clade_rectangular.png"
out_pdf_circ   <- "viptree_OCPT_Ludmila_clade_circular.pdf"
out_png_circ   <- "viptree_OCPT_Ludmila_clade_circular.png"
fig_width      <- 14   # Stretch plot horizontally
fig_height     <- 24   # Rectangular: tall for vertical spacing (153 tips)
fig_size_circ  <- 12   # Circular: square aspect
png_dpi        <- 300

# =============================================================================
# 1) Inspect available data: list files and classify
# =============================================================================
# Resolve paths: if relative, assume run from project root (LudmilaCharacterization)
if (!dir.exists(viptree_dir) && dir.exists(file.path("LudmilaCharacterization", viptree_dir))) {
  viptree_dir <- file.path("LudmilaCharacterization", viptree_dir)
}
if (!dir.exists(results_dir) && dir.exists(file.path("LudmilaCharacterization", results_dir))) {
  results_dir <- file.path("LudmilaCharacterization", results_dir)
}
wd <- normalizePath(viptree_dir, mustWork = TRUE)
message("ViPTree directory (data): ", wd)
message("Results directory (outputs): ", normalizePath(results_dir, mustWork = FALSE))

all_files <- list.files(wd, recursive = TRUE, full.names = FALSE)
# Newick: common extensions
newick_files <- all_files[grepl("\\.newick$|\\.nwk$|\\.tree$", all_files, ignore.case = TRUE)]
# Matrices: common patterns
matrix_files <- all_files[grepl("matrix|dist\\.|sim\\.|similarity|distance", all_files, ignore.case = TRUE)]
# Metadata / label files (TSV, CSV, genome table)
meta_files    <- all_files[grepl("genome\\.tsv|genome_table|metadata|labels?\\.(tsv|csv|txt)", all_files, ignore.case = TRUE)]

# Comment in output (and in code) which are which
message("Newick tree file(s): ", paste(newick_files, collapse = ", "))
message("Matrix file(s): ", paste(matrix_files, collapse = ", "))
message("Metadata/label file(s): ", paste(meta_files, collapse = ", "))

# Resolve paths for primary inputs
tree_path <- file.path(wd, tree_file)
meta_path <- file.path(wd, genome_table_file)
if (!file.exists(tree_path)) stop("Tree file not found: ", tree_path)

# =============================================================================
# 2) Read tree and metadata
# =============================================================================
# Load the primary Newick tree into a phylo object
tree <- read.tree(tree_path)
message("Tree: ", Ntip(tree), " tips, ", Nnode(tree), " internal nodes")

# Load genome table (ViPTree Figure 2b–style: ID, Hname, name, etc.)
meta <- data.frame()
if (file.exists(meta_path)) {
  meta <- read.delim(meta_path, stringsAsFactors = FALSE, check.names = FALSE)
  # Use first column as ID if it's not named ID (ViPTree uses "ID")
  if (!"ID" %in% names(meta) && ncol(meta) > 0) names(meta)[1] <- "ID"
  message("Metadata: ", nrow(meta), " rows, columns: ", paste(head(names(meta), 12), collapse = ", "))
} else {
  # Minimal metadata from tip labels only
  meta <- data.frame(ID = tree$tip.label, stringsAsFactors = FALSE)
  meta$Hname <- NA_character_
  meta$name <- meta$ID
}

# Ensure every tip has a row: merge tree tips with metadata
tip_df <- data.frame(label = tree$tip.label, stringsAsFactors = FALSE)
# Join key: tree tip label may be ID or seq_ids in genome table
if ("ID" %in% names(meta)) {
  tip_meta <- left_join(tip_df, meta, by = c("label" = "ID"))
} else if ("seq_ids" %in% names(meta)) {
  tip_meta <- left_join(tip_df, meta, by = c("label" = "seq_ids"))
} else {
  tip_meta <- tip_df
}
# Host column for coloring (ViPTree uses Hname)
if (!"Hname" %in% names(tip_meta)) tip_meta$Hname <- NA_character_
# Flag S. aureus phages (include focal phage even if Hname is missing)
tip_meta$is_Saureus <- grepl(host_filter, tip_meta$Hname, fixed = TRUE) |
  grepl(ocpt_label, tip_meta$label, fixed = TRUE)
# Focal phage for bold/color/point
tip_meta$is_OCPT_Ludmila <- grepl(ocpt_label, tip_meta$label, fixed = TRUE)

# =============================================================================
# 3) Prune tree to S. aureus phages only
# =============================================================================
tips_to_keep <- tip_meta$label[tip_meta$is_Saureus]
n_saureus    <- length(tips_to_keep)
if (n_saureus == 0) stop("No S. aureus phages found. Check host_filter and metadata Hname.")
tree <- keep.tip(tree, tips_to_keep)
tip_meta <- tip_meta %>% filter(label %in% tree$tip.label)
message("Pruned to S. aureus phages only: ", Ntip(tree), " tips")

# =============================================================================
# 4) Locate OCPT-Ludmila and its clade node (in pruned tree)
# =============================================================================
ocpt_tip_idx <- which(tip_meta$is_OCPT_Ludmila)
if (length(ocpt_tip_idx) == 0) stop("No tip matching '", ocpt_label, "' found. Check ocpt_label and tree tip labels.")

# Smallest monophyletic clade containing OCPT-Ludmila = parent node of its branch
# In ape: edge[,2] is descendant; edge[,1] is ancestor
ocpt_tip_node <- ocpt_tip_idx  # tip index = node number for tips 1..Ntip
parent_edges  <- which(tree$edge[, 2] == ocpt_tip_node)
clade_node    <- tree$edge[parent_edges[1], 1]
message("OCPT-Ludmila tip index: ", ocpt_tip_node, "; clade (parent) node: ", clade_node)

# Phages in the photo (OCPT-Ludmila clade list): tip labels from that panel
# These match the ~24 phages in the ViPTree clade view (Rountreeviridae / Rosenblumvirus group)
tips_in_clade_photo <- c(
  "NC_047814", "NC_047813", "NC_070880", "NC_070879", "NC_055802", "NC_047919",
  "NC_023550", "NC_009875", "NC_048159", "NC_055916", "NC_070878", "NC_055814",
  "NC_055906", "NC_047855", "NC_007046", "NC_004678", "NC_004679", "NC_031046",
  "NC_016565", "OCPT-Ludmila-scaffold_1", "NC_031008", "AB626963", "NC_047788", "NC_048107"
)
# Keep only those present in the pruned tree
tips_in_clade <- intersect(tips_in_clade_photo, tree$tip.label)
message("Phages in photo / clade (to box, excluding OCPT-Ludmila): ", length(tips_in_clade))

# =============================================================================
# 5) Tip metadata: display names (strip "Staphylococcus phage "), bold vs boxed
# =============================================================================
# OCPT-Ludmila: bold + red only (no point). Box: clade (photo) except OCPT-Ludmila, plus outgroups K, Romulus, EW.
boxed_phages   <- c("K", "EW")
boxed_contains  <- "Romulus"

if (!"name" %in% names(tip_meta)) tip_meta$name <- NA_character_
strip_prefix <- function(x) gsub("^Staphylococcus phage ", "", x, ignore.case = TRUE)
tip_meta <- tip_meta %>% mutate(
  name_short = case_when(
    is_OCPT_Ludmila ~ "OCPT-Ludmila",
    !is.na(name) & name != "" & name != "user_virus1" ~ strip_prefix(name),
    TRUE ~ label
  ),
  display_label = name_short,
  in_clade_photo = label %in% tips_in_clade,
  is_boxed = !is_OCPT_Ludmila & (
    in_clade_photo |
    (display_label %in% boxed_phages) |
    grepl(boxed_contains, display_label, fixed = TRUE)
  ),
  # Names to highlight (grey box on label only): clade + outgroups
  is_highlight_name = in_clade_photo | (display_label %in% boxed_phages) | grepl(boxed_contains, display_label, fixed = TRUE),
  label_color = if_else(is_OCPT_Ludmila, "#B22222", "black"),
  label_font = if_else(is_OCPT_Ludmila, "bold", "plain"),
  label_size = if_else(is_OCPT_Ludmila, 3.5, 3)
)
rownames(tip_meta) <- tip_meta$label

# Outgroup tip labels (K, Romulus, EW)
outgroup_labels <- tip_meta$label[tip_meta$display_label %in% boxed_phages | grepl(boxed_contains, tip_meta$display_label, fixed = TRUE)]

# Manuscript-style theme: serif font, clean, no grid
theme_manuscript <- function(base_size = 11) {
  theme_minimal(base_size = base_size, base_family = "serif") +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    axis.title       = element_blank(),
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    panel.grid       = element_blank(),
    strip.background = element_blank(),
    plot.margin      = margin(6, 6, 6, 6, "mm")
  )
}

# Build tree plot for a given layout (rectangular or circular)
tiplab_offset <- 0.002
label_extent  <- 0.18   # how far right the box extends from label start (over text only)
label_pad_y   <- 0.45   # vertical padding around label

build_tree_plot <- function(layout = "rectangular") {
  tree_fort <- ggtree::fortify(tree)
  # Boxes over label text only (not branches): use tip positions + offset
  clade_tips <- tree_fort %>% filter(.data$isTip, .data$label %in% tips_in_clade)
  xmax_boxes <- max(clade_tips$x) + tiplab_offset + label_extent
  # One big box over clade labels only
  rect_clade <- data.frame(
    xmin = min(clade_tips$x) + tiplab_offset,
    xmax = xmax_boxes,
    ymin = min(clade_tips$y) - label_pad_y,
    ymax = max(clade_tips$y) + label_pad_y
  )
  p <- ggtree(tree, layout = layout, size = 0.35) %<+% tip_meta +
    theme_manuscript(11)
  # Boxes over text only (drawn first so text is on top)
  p <- p + geom_rect(
    data = rect_clade,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "grey88", alpha = 0.7, inherit.aes = FALSE
  )
  # Outgroup boxes: same width (xmax), over each outgroup label only
  outgroup_tips <- tree_fort %>% filter(.data$isTip, .data$label %in% outgroup_labels)
  if (nrow(outgroup_tips) > 0) {
    rect_out <- outgroup_tips %>%
      mutate(
        xmin = .data$x + tiplab_offset,
        xmax = xmax_boxes,
        ymin = .data$y - label_pad_y,
        ymax = .data$y + label_pad_y
      )
    p <- p + geom_rect(
      data = rect_out,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "grey88", alpha = 0.7, inherit.aes = FALSE
    )
  }
  # Box around OCPT-Ludmila label only (same width as other boxes)
  ocpt_label <- tip_meta$label[tip_meta$is_OCPT_Ludmila]
  ocpt_tip <- tree_fort %>% filter(.data$isTip, .data$label %in% ocpt_label)
  if (nrow(ocpt_tip) > 0) {
    rect_ocpt <- ocpt_tip %>%
      mutate(
        xmin = .data$x + tiplab_offset,
        xmax = xmax_boxes,
        ymin = .data$y - label_pad_y,
        ymax = .data$y + label_pad_y
      )
    p <- p + geom_rect(
      data = rect_ocpt,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "grey88", alpha = 0.7, inherit.aes = FALSE
    )
  }
  # Tip labels on top
  p <- p + geom_tiplab(
    aes(label = display_label, color = label_color, fontface = label_font, size = label_size),
    geom = "text",
    align = TRUE,
    linesize = 0.25,
    offset = tiplab_offset,
    family = "serif"
  ) +
    scale_color_identity() +
    scale_size_identity() +
    hexpand(0.25) +
    theme(plot.margin = margin(6, 36, 6, 6, "mm"))
  p
}

# =============================================================================
# 6) Rectangular tree: build and save
# =============================================================================
p_rect <- build_tree_plot("rectangular")

# =============================================================================
# 7) Circular tree: build and save (same data, circular layout)
# =============================================================================
p_circ <- build_tree_plot("circular")

# =============================================================================
# 8) Export PDF and PNG (rectangular and circular)
# =============================================================================
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
# Rectangular
ggsave(file.path(results_dir, out_pdf_rect), plot = p_rect, width = fig_width, height = fig_height, device = "pdf")
message("Saved: ", file.path(results_dir, out_pdf_rect))
ggsave(file.path(results_dir, out_png_rect), plot = p_rect, width = fig_width, height = fig_height, dpi = png_dpi, bg = "white")
message("Saved: ", file.path(results_dir, out_png_rect))
# Circular (square dimensions)
ggsave(file.path(results_dir, out_pdf_circ), plot = p_circ, width = fig_size_circ, height = fig_size_circ, device = "pdf")
message("Saved: ", file.path(results_dir, out_pdf_circ))
ggsave(file.path(results_dir, out_png_circ), plot = p_circ, width = fig_size_circ, height = fig_size_circ, dpi = png_dpi, bg = "white")
message("Saved: ", file.path(results_dir, out_png_circ))

message("Done. Rectangular and circular trees (with phage names, manuscript style) exported.")
