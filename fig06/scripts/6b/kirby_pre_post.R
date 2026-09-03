# Figure 6b — Kirby-Bauer zone diameters, pre-treatment (OCPT04SA01) vs
# post-treatment (OCPT04SA03) clinical isolate, six antibiotics.
#
# Same plot as 6a without the strain facet: here the two isolates are the
# comparison, so Strain and group are the same contrast.
#
# Run with ../../run.sh, which supplies FIG_DATA, FIG_OUT and FIG_SHARED.
#
# Optional overrides (display or policy, never data):
#   KB_FILE, KB_DISC_MM, KB_STARS, KB_ABX_ORDER, KB_GROUP_ORDER — see 6a.

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

INFILE      <- env_or("KB_FILE", "kirby_bauer_pre_post.tsv")
DISC_MM     <- as.numeric(env_or("KB_DISC_MM", "6"))
STAR_P      <- env_or("KB_STARS", "p")
TEST        <- env_or("KB_TEST", "welch")
ABX_ORDER   <- env_csv("KB_ABX_ORDER",   c("VA", "LEV", "MEM", "SXT", "MH", "LZD"))
GROUP_ORDER <- env_csv("KB_GROUP_ORDER", c("Pre", "Post"))

fig_width  <- 10
fig_height <- 5.5

fill_values <- stats::setNames(c("white", "gray20"), GROUP_ORDER)

long <- read_kirby_long(
  file.path(data_dir, INFILE),
  group_col    = "Clone",
  group_levels = GROUP_ORDER,
  abx_order    = ABX_ORDER,
  strain_order = NULL          # strain identity is carried by the group column
)
long <- check_disc_floor(long, DISC_MM)

stats_df <- kirby_welch(long, by_strain = FALSE, test = TEST)
utils::write.table(
  stats_df,
  file.path(results_dir, "kirby_pre_post_welch.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)

skipped <- stats_df[!is.na(stats_df$skipped), , drop = FALSE]
if (nrow(skipped)) {
  message("Comparisons not tested (", nrow(skipped), "):")
  for (i in seq_len(nrow(skipped))) {
    message("  ", skipped$Antibiotic[i], " — ", skipped$skipped[i])
  }
}

# n = 3 per group: Welch on three replicates has very little power, and the
# replicates are technical. Recorded here so the limitation travels with the
# panel rather than living only in the Methods.
message("Replicates per group: ",
        paste(sprintf("%s n1=%d n2=%d", stats_df$Antibiotic, stats_df$n1, stats_df$n2),
              collapse = "; "))

sum_df <- summarize_kirby(long, by_strain = FALSE)
ann    <- kirby_brackets(stats_df, long, sum_df, p_col = STAR_P, by_strain = FALSE)

p <- kirby_barplot(
  long, sum_df, ann,
  fill_values     = fill_values,
  facet_by_strain = FALSE,
  fig_width       = fig_width,
  fig_height      = fig_height,
  subtitle = paste0(
    TEST, " test ", GROUP_ORDER[1], " vs ", GROUP_ORDER[2],
    "; brackets: ", STAR_P, " <= 0.05; n = 3 independent colonies per isolate",
    " (all p-values in kirby_pre_post_welch.tsv)"
  )
)

save_panel(p, file.path(results_dir, "kirby_pre_post"), fig_width, fig_height)
save_group_legend(fill_values, file.path(results_dir, "kirby_pre_post_legend.png"))
