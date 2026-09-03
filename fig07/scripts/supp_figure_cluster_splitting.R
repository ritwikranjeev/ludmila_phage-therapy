# Supplementary figure — why the pangenome presence/absence matrix cannot
# support gene association claims in this dataset.
#
# Run with ../run.sh, which supplies FIG_DATA and FIG_OUT.
#
# Panaroo assigns a gene to a cluster by sequence similarity. When alleles of one
# gene diverge past that threshold the gene is split across clusters, and each
# arm then looks like an accessory gene absent from the strains carrying the
# other arm. Testing the arms separately produces two significant associations in
# OPPOSITE directions; testing the gene produces nothing.
#
# rplI is the clearest case: a 50S ribosomal protein, present in 82 of 83 strains,
# split into arms of 24 and 59. metXA behaves identically. sraP is the case that
# matters for this paper, because the published association is one arm of a split.
#
# Input: cluster_splitting_stats.tsv from scripts/supp_outputs.py

if (!interactive()) pdf(NULL)

suppressPackageStartupMessages({
  library(ggplot2)
})

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

# FIG_SUBTITLE=0 suppresses the methods subtitle without removing it — the text
# is still built by the caller, it just is not drawn.
FIG_SUBTITLE <- as.integer(env_or("FIG_SUBTITLE", 1))

data_dir    <- need("FIG_DATA")
results_dir <- need("FIG_OUT")
# This panel's input is an analysis intermediate, so it comes from work/, not
# data/. FIG_DATA used to point at work/; it now points at data/, hence FIG_WORK.
work_dir    <- env_or("FIG_WORK", file.path(dirname(data_dir), "work"))
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

STATS_FILE <- env_or("F7_SPLIT_STATS", "cluster_splitting_stats.tsv")
ALPHA      <- as.numeric(env_or("F7_SPLIT_ALPHA", "0.05"))

d <- utils::read.delim(file.path(work_dir, STATS_FILE), stringsAsFactors = FALSE)
d$p <- as.numeric(d$p)
d$feature <- factor(d$feature, levels = c("merged", "arm 2", "arm 1"))
d$significant <- ifelse(d$p < ALPHA, "p < 0.05", "not significant")
d$label <- sprintf("n = %d,  p = %s", d$strains_with,
                   formatC(d$p, format = "g", digits = 2))

theme_manuscript <- function(base_size = 11) {
  theme_minimal(base_size = base_size, base_family = "Times New Roman") +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      strip.text = element_text(face = "bold", size = base_size - 1, hjust = 0),
      plot.title = element_text(face = "bold", size = base_size + 1, hjust = 0),
      plot.subtitle = element_text(size = base_size - 2, colour = "grey30", hjust = 0),
      axis.text = element_text(colour = "black"),
      axis.line.x = element_line(colour = "black", linewidth = 0.4),
      axis.ticks.x = element_line(colour = "black", linewidth = 0.4)
    )
}

rng <- max(abs(d$beta)) * 1.9

p <- ggplot(d, aes(x = beta, y = feature)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = beta, y = feature, yend = feature),
               colour = "grey70", linewidth = 0.4) +
  geom_point(aes(fill = significant, shape = feature), size = 3.1,
             colour = "grey20", stroke = 0.35) +
  geom_text(aes(label = label,
                hjust = ifelse(beta > 0, -0.12, 1.12)),
            size = 2.5, family = "Times New Roman", colour = "grey25") +
  facet_wrap(~ gene, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("p < 0.05" = "#d1495b",
                               "not significant" = "grey85"), name = NULL) +
  scale_shape_manual(values = c("arm 1" = 21, "arm 2" = 21, "merged" = 23),
                     guide = "none") +
  scale_x_continuous(limits = c(-rng, rng)) +
  labs(x = "Effect on susceptibility (mixed-model beta)", y = NULL,
       title = "Splitting one gene across clusters manufactures associations",
       subtitle = if (FIG_SUBTITLE)
         paste("Each arm of a split gene associates significantly and in",
               "opposite directions; the merged gene does not.") else NULL) +
  theme_manuscript(11) +
  theme(legend.position = "top", legend.justification = "left",
        legend.text = element_text(size = 8))

base <- file.path(results_dir, "supp_figure7_cluster_splitting")
ggsave(paste0(base, ".png"), p, width = 7.2, height = 6.4, units = "in",
       dpi = FIG_DPI, bg = "white")
if (save_pdf_best_effort(p, paste0(base, ".pdf"), width = 7.2, height = 6.4,
                         units = "in", bg = "white")) {
  message("Wrote: ", base, ".png and .pdf")
} else {
  message("Wrote: ", base, ".png   (no .pdf: cairo_pdf wrote no file)")
}
