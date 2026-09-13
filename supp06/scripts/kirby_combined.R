# Supplementary Figure 6 — Kirby-Bauer zone diameters for the four in vitro
# backgrounds not shown in the main figure, in a 2 x 2 grid.
#
# WHICH STRAINS AND WHY
#
# OCPT04SA01 and the clinical pre/post pair are deliberately absent: they are
# Figure 6c and Figure 6b respectively, each drawn full size in the main figure.
# The supplement carries the remaining four backgrounds — USA300, SA412,
# OCPT02SA01 and PedCF — so nothing is duplicated between the two figures.
#
# Run with ../run.sh, which supplies FIG_DATA, FIG_OUT and FIG_SHARED.
#
# This figure has NO data and NO statistics of its own. It reads Figure 6's
# input file and Figure 6's helpers, so the test and the multiplicity
# correction cannot drift between the main figure and the supplement.
#
# THE CORRECTION FAMILY IS THE FULL IN VITRO SET, NOT THE FOUR DRAWN
#
# The five-strain sheet is read and tested in full, and only four strains are
# drawn. Testing the subset instead would make these four their own family of
# 24 comparisons and the adjusted p-values would stop matching Supplementary
# Table 3 and the Methods, which state that the in vitro comparisons are
# corrected together.
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
# p_adj_all is the per-panel Benjamini-Hochberg value carried over from
# kirby_welch: across the 29 testable in vitro comparisons for the mutant
# panels, and across the 6 clinical comparisons for the clinical panel. That is
# the two-family correction the manuscript reports (clinical meropenem q =
# 0.005). The alternative, p_adj_figure, pools all 36 into one family and gives
# q = 0.028 for the same comparison. Both are defensible; they must not be
# mixed between Figure 6 and this supplement.
STAR_P       <- env_or("KB_STARS", "p_adj_all")
ABX_ORDER    <- env_csv("KB_ABX_ORDER",    c("VA", "LEV", "MEM", "SXT", "MH", "LZD"))
# Strains DRAWN. OCPT04SA01 is excluded because it is Figure 6c; the clinical
# pair is excluded because it is Figure 6b. The statistics below are still
# computed over the full five-strain sheet — see the header.
STRAIN_ORDER <- env_csv("KB_STRAIN_ORDER", c("USA300", "SA412", "OCPT02SA01", "PedCF"))
NCOL         <- as.integer(env_or("S6_NCOL", "2"))

panel_width  <- 6
panel_height <- 4.5

if (!requireNamespace("cowplot", quietly = TRUE)) {
  stop("Package cowplot is required to assemble the grid: install.packages(\"cowplot\")")
}

# ---- read both inputs through Figure 6's reader ----------------------------

# strain_order is NULL deliberately: every strain in the sheet is read, so the
# test below runs over the full in vitro set. STRAIN_ORDER selects what is drawn.
mutants <- read_kirby_long(
  file.path(data_dir, "6a", env_or("KB_FILE_6A", "kirby_bauer_in_vitro_mutants.tsv")),
  group_col = "Clone", group_levels = c("Ancestor", "Resistant"),
  abx_order = ABX_ORDER, strain_order = NULL
)
mutants <- check_disc_floor(mutants, DISC_MM)

missing <- setdiff(STRAIN_ORDER, levels(mutants$Strain))
if (length(missing)) {
  stop("KB_STRAIN_ORDER names strains not in the sheet: ",
       paste(missing, collapse = ", "),
       ". Available: ", paste(levels(mutants$Strain), collapse = ", "))
}

# ---- statistics over the whole grid ----------------------------------------

# Tested over every strain in the sheet, not just the four drawn, so p_adj_all
# is the same in vitro correction family the manuscript reports.
all_stats <- kirby_welch(mutants, by_strain = TRUE, test = TEST)

# Alternative correction, computed over the in vitro comparisons only. The
# clinical pair is not in this figure, so it is not in this family either.
all_stats$p_adj_figure <- stats::p.adjust(all_stats$p, method = "BH")

utils::write.table(all_stats, file.path(results_dir, "kirby_combined_stats.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

n_tested <- sum(!is.na(all_stats$p))
message("Comparisons tested: ", n_tested, " across ",
        length(levels(mutants$Strain)), " strains; drawn: ",
        length(STRAIN_ORDER), " strains (", paste(STRAIN_ORDER, collapse = ", "), ")")
message("Nominally p <= 0.05: ",
        sum(all_stats$p <= 0.05, na.rm = TRUE),
        " | BH across the in vitro set q < 0.05: ",
        sum(all_stats$p_adj_figure < 0.05, na.rm = TRUE))
sig <- all_stats[!is.na(all_stats$p_adj_figure) & all_stats$p_adj_figure < 0.05, , drop = FALSE]
if (nrow(sig)) {
  for (i in seq_len(nrow(sig))) {
    message("  ", sig$Strain[i], " ", sig$Antibiotic[i],
            "  p = ", signif(sig$p[i], 3), "  q = ", signif(sig$p_adj_figure[i], 3))
  }
}

# ---- one panel per drawn strain ---------------------------------------------

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
    theme(plot.title = element_text(face = "bold", size = 12 * FIG_TEXT_SCALE, hjust = 0.5,
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

# Only one legend now: every panel in this figure is an ancestor/resistant
# comparison. The Pre/Post legend belonged to the clinical panel, which has
# moved to Figure 6b and carries its own.
save_group_legend(stats::setNames(c("white", "gray20"), c("Ancestor", "Resistant")),
                  file.path(results_dir, "kirby_combined_legend.png"))
