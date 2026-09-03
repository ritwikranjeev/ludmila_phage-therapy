# Supplementary Figure 2a — ViPTree proteomic tree of Staphylococcus phages,
# with OCPT-Ludmila and its clade highlighted.
#
# Run with ../../run.sh, which supplies FIG_DATA and FIG_OUT.
#
# The clade highlight is COMPUTED, not transcribed. The original script carried
# a hand-copied list of 24 accessions "matching the phages in the ViPTree clade
# view"; four of them are Staphylococcus epidermidis phages that the host filter
# removes, so the drawn box silently contained 20. Here the clade is defined by
# the genomic similarity ViPTree itself reports (S_G to the focal genome), and
# the script checks that the resulting set is monophyletic before drawing it.
#
# Optional overrides:
#   S2_TREE_FILE   Newick file inside FIG_DATA (default set.bionj.newick)
#   S2_META_FILE   ViPTree genome table inside FIG_DATA (default genome.tsv)
#   S2_HOST        host string that defines the tips to keep
#   S2_FOCAL       substring identifying the focal phage tip
#   S2_SG_MIN      minimum S_G to the focal genome for clade membership
#   S2_OUTGROUPS   comma-separated display names to box as outgroups
#   S2_LAYOUTS     comma-separated ggtree layouts to render

if (!interactive()) pdf(NULL)

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggplot2)
})

options(ignore.negative.edge = TRUE)   # ViPTree can emit negative edge lengths

# Typography. The set is standardised on Times New Roman; "Times" and "serif"
# resolve to a DIFFERENT file on macOS (Times.ttc, the Linotype face), so mixing
# the names silently sets panels in two faces. Fall back only if the font is
# genuinely absent, e.g. on a Linux reproducibility run.
FONT_FAMILY <- "Times New Roman"
if (requireNamespace("systemfonts", quietly = TRUE)) {
  if (!FONT_FAMILY %in% unique(systemfonts::system_fonts()$family)) FONT_FAMILY <- "Times"
} else {
  FONT_FAMILY <- "Times"
}

need <- function(v) {
  x <- Sys.getenv(v, unset = NA_character_)
  if (is.na(x)) stop(v, " must be set. Run this figure with ./run.sh")
  x
}
env_or <- function(v, default) {
  x <- Sys.getenv(v, unset = "")
  if (!nzchar(x)) default else x
}

# Raster resolution for every PNG in this figure. Panels are assembled in Canva
# and exported to PDF from there, so PNG is the delivered output. Width and
# height are NOT configurable: they match the published originals.
FIG_DPI <- as.numeric(env_or("FIG_DPI", 600))

# Best-effort vector companion to the PNG.
#
# cairo_pdf() can fail at runtime with a WARNING rather than an error: it writes
# no file, tryCatch(error=) never fires, and the caller would otherwise report a
# PDF that does not exist (see CLAUDE.md, Environment). So attempt the write and
# then confirm a fresh, non-empty file actually landed. Returns TRUE only then.
#
# There is deliberately no base pdf() fallback: that device rejects
# 'Times New Roman' and would reset the panel in a second face.
save_pdf_best_effort <- function(plot, path, ...) {
  t0 <- Sys.time()
  tryCatch(
    suppressWarnings(ggsave(path, plot = plot, device = grDevices::cairo_pdf, ...)),
    error = function(e) invisible(NULL))
  isTRUE(file.exists(path) && file.size(path) > 0 && file.mtime(path) >= t0)
}
env_csv <- function(v, default) {
  x <- Sys.getenv(v, unset = "")
  if (!nzchar(x)) return(default)
  trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
}

data_dir    <- need("FIG_DATA")
results_dir <- need("FIG_OUT")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

TREE_FILE <- env_or("S2_TREE_FILE", "set.bionj.newick")
META_FILE <- env_or("S2_META_FILE", "genome.tsv")
HOST      <- env_or("S2_HOST", "Staphylococcus aureus")
FOCAL     <- env_or("S2_FOCAL", "OCPT-Ludmila")
SG_MIN    <- as.numeric(env_or("S2_SG_MIN", "0.1"))
OUTGROUPS <- env_csv("S2_OUTGROUPS", c("K", "EW", "Romulus"))
LAYOUTS   <- env_csv("S2_LAYOUTS", c("rectangular", "circular"))

# ---- read ------------------------------------------------------------------

tree <- ape::read.tree(file.path(data_dir, TREE_FILE))
message("Tree: ", ape::Ntip(tree), " tips, ", ape::Nnode(tree), " internal nodes")

meta <- utils::read.delim(file.path(data_dir, META_FILE), check.names = FALSE,
                          stringsAsFactors = FALSE)
if (!"ID" %in% names(meta)) names(meta)[1L] <- "ID"

focal_tips <- grep(FOCAL, tree$tip.label, fixed = TRUE, value = TRUE)
if (length(focal_tips) != 1L) {
  stop("Expected exactly one tip matching '", FOCAL, "', found ", length(focal_tips))
}
focal_id <- focal_tips[1L]

# ViPTree names this column "<i>S</i><sub>G</sub> to <focal>"; match on the tail.
sg_col <- grep(paste0("to ", focal_id), names(meta), fixed = TRUE, value = TRUE)
if (length(sg_col) != 1L) {
  stop("Could not find the S_G column for ", focal_id,
       " in ", META_FILE, " (columns: ", paste(names(meta), collapse = ", "), ")")
}
meta$sg <- suppressWarnings(as.numeric(meta[[sg_col]]))

idx  <- match(tree$tip.label, meta$ID)
host <- ifelse(is.na(idx), NA_character_, meta$Hname[idx])
sg   <- ifelse(is.na(idx), NA_real_, meta$sg[idx])
names(sg) <- tree$tip.label

# ---- prune to the host of interest -----------------------------------------

is_host  <- grepl(HOST, host, fixed = TRUE)
is_focal <- tree$tip.label == focal_id
keep     <- tree$tip.label[is_host | is_focal]

# Genomes similar to the focal phage that the host filter removes. These are
# real relatives; the figure excludes them by choice of host, not by distance.
excluded_relatives <- tree$tip.label[!(is_host | is_focal) & !is.na(sg) & sg >= SG_MIN]
if (length(excluded_relatives)) {
  message("Close relatives excluded by the '", HOST, "' host filter (",
          length(excluded_relatives), "):")
  for (t in excluded_relatives) {
    r <- meta[match(t, meta$ID), , drop = FALSE]
    message("  ", t, "  S_G = ", signif(sg[[t]], 3), "  host = ", r$Hname,
            "  ", r$name)
  }
}

tree <- ape::keep.tip(tree, keep)
message("Pruned to ", ape::Ntip(tree), " tips matching host '", HOST, "' (plus the focal phage)")

# ---- clade -----------------------------------------------------------------

sg_pruned  <- sg[tree$tip.label]
clade_tips <- tree$tip.label[!is.na(sg_pruned) & sg_pruned >= SG_MIN]
if (length(clade_tips) < 2L) {
  stop("S_G >= ", SG_MIN, " selects ", length(clade_tips),
       " tip(s); nothing to highlight. Check S2_SG_MIN.")
}

mrca       <- ape::getMRCA(tree, clade_tips)
clade_desc <- ape::extract.clade(tree, mrca)$tip.label
extra      <- setdiff(clade_desc, clade_tips)

message("Clade highlighted: ", length(clade_tips), " tips with S_G >= ", SG_MIN,
        " (range ", signif(min(sg_pruned[clade_tips]), 3), " to ",
        signif(max(sg_pruned[clade_tips]), 3), ")")
if (length(extra)) {
  message("NOTE: the smallest clade containing them also contains ", length(extra),
          " tip(s) below the threshold, so the highlighted set is not monophyletic:")
  for (t in extra) message("  ", t, "  S_G = ", signif(sg_pruned[[t]], 3))
  message("Box drawn over the monophyletic clade, which includes those tips.")
  clade_tips <- clade_desc
}

# ---- tip labels ------------------------------------------------------------

midx <- match(tree$tip.label, meta$ID)
name <- ifelse(is.na(midx), NA_character_, meta$name[midx])
strip_prefix <- function(x) sub("^Staphylococcus phage ", "", x, ignore.case = TRUE)

display <- ifelse(
  tree$tip.label == focal_id, FOCAL,
  ifelse(!is.na(name) & nzchar(name) & name != "user_virus1",
         strip_prefix(name), tree$tip.label)
)

tip_meta <- data.frame(
  label         = tree$tip.label,
  display_label = display,
  is_focal      = tree$tip.label == focal_id,
  in_clade      = tree$tip.label %in% clade_tips,
  stringsAsFactors = FALSE
)
tip_meta$is_outgroup <- !tip_meta$is_focal &
  (tip_meta$display_label %in% OUTGROUPS |
     Reduce(`|`, lapply(OUTGROUPS, function(o) grepl(o, tip_meta$display_label, fixed = TRUE))))
tip_meta$label_color <- ifelse(tip_meta$is_focal, "#B22222", "black")
tip_meta$label_font  <- ifelse(tip_meta$is_focal, "bold", "plain")
tip_meta$label_size  <- ifelse(tip_meta$is_focal, 3.5, 3)

found_outgroups <- tip_meta$display_label[tip_meta$is_outgroup]
message("Outgroups boxed: ",
        if (length(found_outgroups)) paste(found_outgroups, collapse = ", ") else "none found")

# ---- plot ------------------------------------------------------------------

theme_manuscript <- function(base_size = 11) {
  theme_minimal(base_size = base_size, base_family = FONT_FAMILY) +
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

tiplab_offset <- 0.002
label_extent  <- 0.18
label_pad_y   <- 0.45

build_tree_plot <- function(layout) {
  fort <- ggtree::fortify(tree)
  tips <- fort[fort$isTip, , drop = FALSE]

  clade_rows <- tips[tips$label %in% clade_tips, , drop = FALSE]
  xmax_boxes <- max(clade_rows$x) + tiplab_offset + label_extent

  box <- function(rows) {
    data.frame(
      xmin = rows$x + tiplab_offset,
      xmax = xmax_boxes,
      ymin = rows$y - label_pad_y,
      ymax = rows$y + label_pad_y
    )
  }

  # One box over the clade labels, one per outgroup label, one over the focal.
  rect_clade <- data.frame(
    xmin = min(clade_rows$x) + tiplab_offset,
    xmax = xmax_boxes,
    ymin = min(clade_rows$y) - label_pad_y,
    ymax = max(clade_rows$y) + label_pad_y
  )
  out_rows   <- tips[tips$label %in% tip_meta$label[tip_meta$is_outgroup], , drop = FALSE]
  focal_rows <- tips[tips$label == focal_id, , drop = FALSE]

  # `size` for line width was deprecated in ggplot2 3.4.0; ggtree forwards it to
  # geom_tree, so the warning surfaced from there. `linewidth` builds an
  # identical layer (verified: same linewidth value, no warning).
  p <- ggtree(tree, layout = layout, linewidth = 0.35) %<+% tip_meta +
    theme_manuscript(11)

  for (r in list(rect_clade,
                 if (nrow(out_rows))   box(out_rows)   else NULL,
                 if (nrow(focal_rows)) box(focal_rows) else NULL)) {
    if (is.null(r)) next
    p <- p + geom_rect(data = r,
                       aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                       fill = "grey88", alpha = 0.7, inherit.aes = FALSE)
  }

  p +
    geom_tiplab(
      aes(label = display_label, colour = label_color,
          fontface = label_font, size = label_size),
      geom = "text", align = TRUE, linesize = 0.25,
      offset = tiplab_offset, family = FONT_FAMILY
    ) +
    scale_colour_identity() +
    scale_size_identity() +
    hexpand(0.25) +
    theme(plot.margin = margin(6, 36, 6, 6, "mm"))
}

# Rectangular is tall (one row per tip); circular is square.
dims <- list(rectangular = c(14, 24), circular = c(12, 12))

for (layout in LAYOUTS) {
  d <- if (!is.null(dims[[layout]])) dims[[layout]] else c(12, 12)
  p <- build_tree_plot(layout)
  out_base <- file.path(results_dir, paste0("proteomic_tree_", layout))
  ggsave(paste0(out_base, ".png"), plot = p, width = d[1], height = d[2],
         dpi = FIG_DPI, bg = "white", limitsize = FALSE)
  if (save_pdf_best_effort(p, paste0(out_base, ".pdf"), width = d[1], height = d[2],
                           limitsize = FALSE)) {
    message("Wrote: ", out_base, ".png and .pdf")
  } else {
    message("Wrote: ", out_base, ".png   (no .pdf: cairo_pdf wrote no file)")
  }
}
