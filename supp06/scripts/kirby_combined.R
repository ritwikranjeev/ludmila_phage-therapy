# Supplementary Figure 6 — Kirby-Bauer zone diameters for every strain in one
# 3 x 2 grid: the five in vitro ancestor/resistant comparisons plus the clinical
# pre/post pair.
#
# Run with ../run.sh, which supplies FIG_DATA, FIG_OUT and FIG_SHARED.
#
# This figure has NO data and NO statistics of its own. It reads Figure 6's two
# input files and Figure 6's helpers, so the test and the multiplicity
# correction cannot drift between the main figure and the supplement. The
# original kept a third near-duplicate script (kirby_combined.R, 350 lines) that
# carried its own copy of the na_if(diameter_mm, 0) bug.
#
# One thing IS computed differently here, deliberately: the supplement shows all
# 36 comparisons at once, so Benjamini-Hochberg is applied across the whole grid
# rather than within each panel. That is the honest correction level for a
# figure a reader takes in as one display item.
#
# Optional overrides:
#   KB_DISC_MM, KB_TEST, KB_STARS, KB_ABX_ORDER, KB_STRAIN_ORDER — as in
#   Figure 6
#   S6_NCOL          panels per row (default 2)

need <- function(v) {
  x <- Sys.getenv(v, unset = NA_character_)
  if (is.na(x)) stop(v, " must be set. Run this figure with ./run.sh")
  x
}
data_dir    <- need("FIG_DATA")      # Figure 6's data/ directory
results_dir <- need("FIG_OUT")
shared      <- need("FIG_SHARED")    # Figure 6's scripts/_shared.R
source(shared)

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

DISC_MM      <- as.numeric(env_or("KB_DISC_MM", "6"))
TEST         <- env_or("KB_TEST", "welch")
STAR_P       <- env_or("KB_STARS", "p")
ABX_ORDER    <- env_csv("KB_ABX_ORDER",    c("VA", "LEV", "MEM", "SXT", "MH", "LZD"))
STRAIN_ORDER <- env_csv("KB_STRAIN_ORDER", c("USA300", "SA412", "OCPT02SA01", "OCPT04SA01", "PedCF"))
NCOL         <- as.integer(env_or("S6_NCOL", "2"))

panel_width  <- 6
panel_height <- 4.5

if (!requireNamespace("cowplot", quietly = TRUE)) {
  stop("Package cowplot is required to assemble the grid: install.packages(\"cowplot\")")
}

# ---- read both inputs through Figure 6's reader ----------------------------

mutants <- read_kirby_long(
  file.path(data_dir, "6a", env_or("KB_FILE_6A", "kirby_bauer_in_vitro_mutants.tsv")),
  group_col = "Clone", group_levels = c("Ancestor", "Resistant"),
  abx_order = ABX_ORDER, strain_order = STRAIN_ORDER
)
clinical <- read_kirby_long(
  file.path(data_dir, "6b", env_or("KB_FILE_6B", "kirby_bauer_pre_post.tsv")),
  group_col = "Clone", group_levels = c("Pre", "Post"),
  abx_order = ABX_ORDER, strain_order = NULL
)

mutants  <- check_disc_floor(mutants,  DISC_MM)
clinical <- check_disc_floor(clinical, DISC_MM)

# ---- statistics over the whole grid ----------------------------------------

st_mut <- kirby_welch(mutants,  by_strain = TRUE,  test = TEST)
st_cln <- kirby_welch(clinical, by_strain = FALSE, test = TEST)
st_cln$Strain <- "Clinical"

common <- intersect(names(st_mut), names(st_cln))
all_stats <- rbind(st_mut[, common, drop = FALSE], st_cln[, common, drop = FALSE])

# One correction across every comparison shown in this figure.
all_stats$p_adj_figure <- stats::p.adjust(all_stats$p, method = "BH")

utils::write.table(all_stats, file.path(results_dir, "kirby_combined_stats.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

n_tested <- sum(!is.na(all_stats$p))
message("Comparisons drawn: ", nrow(all_stats), "; tested: ", n_tested,
        "; test: ", TEST)
message("Nominally p <= 0.05: ",
        sum(all_stats$p <= 0.05, na.rm = TRUE),
        " | BH across the whole figure q < 0.05: ",
        sum(all_stats$p_adj_figure < 0.05, na.rm = TRUE))
sig <- all_stats[!is.na(all_stats$p_adj_figure) & all_stats$p_adj_figure < 0.05, , drop = FALSE]
if (nrow(sig)) {
  for (i in seq_len(nrow(sig))) {
    message("  ", sig$Strain[i], " ", sig$Antibiotic[i],
            "  p = ", signif(sig$p[i], 3), "  q = ", signif(sig$p_adj_figure[i], 3))
  }
}

# ---- one panel per strain, plus the clinical pair ---------------------------

# Stars in the supplement follow the figure-wide correction unless KB_STARS
# names a different column.
star_col <- if (STAR_P == "p") "p" else if (STAR_P %in% names(all_stats)) STAR_P else "p_adj_figure"

build_panel <- function(long, stats_df, strain, title, by_strain, show_y) {
  sub_long <- if (by_strain) {
    d <- long[as.character(long$Strain) == strain, , drop = FALSE]
    d$Strain <- factor(as.character(d$Strain), levels = strain)
    d
  } else long
  sub_stats <- stats_df[as.character(stats_df$Strain) == strain, , drop = FALSE]

  sum_df <- summarize_kirby(sub_long, by_strain = FALSE)
  ann    <- kirby_brackets(sub_stats, sub_long, sum_df, p_col = star_col,
                           by_strain = FALSE)

  fill_values <- stats::setNames(c("white", "gray20"), levels(sub_long$Group))
  p <- kirby_barplot(sub_long, sum_df, ann,
                     fill_values = fill_values, facet_by_strain = FALSE,
                     fig_width = panel_width, fig_height = panel_height) +
    ggtitle(title) +
    theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5,
                                    margin = margin(b = 6, unit = "pt")))
  if (!show_y) p <- p + labs(y = NULL)
  p
}

panels <- list()
for (s in STRAIN_ORDER) {
  panels[[length(panels) + 1L]] <- build_panel(
    mutants, all_stats, s, s, by_strain = TRUE,
    show_y = (length(panels) %% NCOL == 0L))
}
panels[[length(panels) + 1L]] <- build_panel(
  clinical, all_stats, "Clinical", "Clinical pre vs post", by_strain = FALSE,
  show_y = (length(panels) %% NCOL == 0L))

nrow_grid <- ceiling(length(panels) / NCOL)
combined  <- cowplot::plot_grid(plotlist = panels, ncol = NCOL, align = "hv",
                                labels = NULL)

fig_width  <- panel_width * NCOL
fig_height <- panel_height * nrow_grid

out_base <- file.path(results_dir, "kirby_combined")
ggsave(paste0(out_base, ".png"), combined, width = fig_width, height = fig_height,
       units = "in", dpi = FIG_DPI, bg = "white", limitsize = FALSE)
if (save_pdf_best_effort(combined, paste0(out_base, ".pdf"), width = fig_width,
                         height = fig_height, units = "in", bg = "white",
                         limitsize = FALSE)) {
  message("Wrote: ", out_base, ".png and .pdf")
} else {
  message("Wrote: ", out_base, ".png   (no .pdf: cairo_pdf wrote no file)")
}

save_group_legend(stats::setNames(c("white", "gray20"), c("Ancestor", "Resistant")),
                  file.path(results_dir, "kirby_combined_legend_mutants.png"))
save_group_legend(stats::setNames(c("white", "gray20"), c("Pre", "Post")),
                  file.path(results_dir, "kirby_combined_legend_clinical.png"))
