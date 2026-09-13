# Figure 6c — Kirby-Bauer zone diameters for OCPT04SA01, ancestor vs its
# phage-resistant mutants, drawn in the same single-panel format as 6b.
#
# WHY THIS PANEL EXISTS
#
# 6a shows all five in vitro backgrounds faceted together, which makes any one
# of them small. OCPT04SA01 is the patient's own pre-therapy isolate, so its
# ancestor-vs-resistant comparison is the direct in vitro counterpart to the
# clinical pre/post comparison in 6b. Drawing it at the same size and in the
# same format lets the two be read side by side.
#
# THE STATISTICS ARE NOT RECOMPUTED ON THE SUBSET
#
# This is the point that matters. Filtering to one strain and then running the
# test would make OCPT04SA01 its own correction family of six comparisons, and
# the adjusted p-values would no longer match Supplementary Table 3 or the
# Methods, which state that the in vitro comparisons are corrected together.
# So the full five-strain table is read and tested exactly as 6a does, and only
# the OCPT04SA01 rows are drawn. The q-values shown are therefore the same
# numbers the manuscript reports.
#
# Run with ../../run.sh, which supplies FIG_DATA, FIG_OUT and FIG_SHARED.
#
# Optional overrides (display or policy, never data):
#   KB_FILE, KB_DISC_MM, KB_STARS, KB_ABX_ORDER, KB_GROUP_ORDER — see 6a.
#   KB_PANEL_STRAIN  which strain to draw (default OCPT04SA01)

need <- function(v) {
  x <- Sys.getenv(v, unset = NA_character_)
  if (is.na(x)) stop(v, " must be set. Run this figure with ./run.sh")
  x
}
data_dir    <- need("FIG_DATA")
results_dir <- need("FIG_OUT")
shared      <- need("FIG_SHARED")
source(shared)

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

INFILE      <- env_or("KB_FILE", "kirby_bauer_in_vitro_mutants.tsv")
DISC_MM     <- as.numeric(env_or("KB_DISC_MM", "6"))
STAR_P      <- env_or("KB_STARS", "p_adj_all")
TEST        <- env_or("KB_TEST", "welch")
STRAIN      <- env_or("KB_PANEL_STRAIN", "OCPT04SA01")
ABX_ORDER   <- env_csv("KB_ABX_ORDER",   c("VA", "LEV", "MEM", "SXT", "MH", "LZD"))
GROUP_ORDER <- env_csv("KB_GROUP_ORDER", c("Ancestor", "Resistant"))

# Matches 6b so the two panels sit together at the same scale.
fig_width  <- 10
fig_height <- 5.5

fill_values <- stats::setNames(c("white", "gray20"), GROUP_ORDER)

# ---- read and test the FULL in vitro set -----------------------------------

long <- read_kirby_long(
  file.path(data_dir, INFILE),
  group_col    = "Clone",
  group_levels = GROUP_ORDER,
  abx_order    = ABX_ORDER,
  strain_order = NULL
)
long <- check_disc_floor(long, DISC_MM)

if (!STRAIN %in% levels(long$Strain)) {
  stop("KB_PANEL_STRAIN '", STRAIN, "' is not in ", INFILE,
       ". Available: ", paste(levels(long$Strain), collapse = ", "))
}

stats_df <- kirby_welch(long, by_strain = TRUE, test = TEST)

# ---- draw one strain --------------------------------------------------------

sub_long <- long[as.character(long$Strain) == STRAIN, , drop = FALSE]
sub_long$Strain <- factor(as.character(sub_long$Strain), levels = STRAIN)
sub_stats <- stats_df[as.character(stats_df$Strain) == STRAIN, , drop = FALSE]

utils::write.table(
  sub_stats,
  file.path(results_dir, paste0("kirby_", tolower(STRAIN), "_welch.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)

skipped <- sub_stats[!is.na(sub_stats$skipped), , drop = FALSE]
if (nrow(skipped)) {
  message("Comparisons not tested (", nrow(skipped), "):")
  for (i in seq_len(nrow(skipped))) {
    message("  ", skipped$Antibiotic[i], " — ", skipped$skipped[i])
  }
}

n_drawn <- sum(!is.na(sub_stats[[STAR_P]]))
n_total <- sum(!is.na(stats_df[[STAR_P]]))
message(STRAIN, ": ", n_drawn, " comparisons drawn; ", STAR_P,
        " computed across all ", n_total, " testable in vitro comparisons, ",
        "not recomputed on this strain alone")
message("Replicates per group: ",
        paste(sprintf("%s n1=%d n2=%d", sub_stats$Antibiotic,
                      sub_stats$n1, sub_stats$n2), collapse = "; "))

sum_df <- summarize_kirby(sub_long, by_strain = FALSE)
ann    <- kirby_brackets(sub_stats, sub_long, sum_df, p_col = STAR_P,
                         by_strain = FALSE)

p <- kirby_barplot(
  sub_long, sum_df, ann,
  fill_values     = fill_values,
  facet_by_strain = FALSE,
  fig_width       = fig_width,
  fig_height      = fig_height,
  subtitle = paste0(
    TEST, " test ", GROUP_ORDER[1], " vs ", GROUP_ORDER[2], " in ", STRAIN,
    "; brackets: ", STAR_P, " <= 0.05 across all ", n_total,
    " testable in vitro comparisons; n = 5 independent resistant clones",
    " (all p-values in kirby_", tolower(STRAIN), "_welch.tsv)"
  )
)

save_panel(p, file.path(results_dir, paste0("kirby_", tolower(STRAIN))),
           fig_width, fig_height)
save_group_legend(fill_values,
                  file.path(results_dir,
                            paste0("kirby_", tolower(STRAIN), "_legend.png")))
