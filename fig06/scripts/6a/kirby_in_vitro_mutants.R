# Figure 6a — Kirby-Bauer zone diameters, phage-resistant mutants vs their
# ancestors, five strains x six antibiotics.
#
# Run with ../../run.sh, which supplies FIG_DATA, FIG_OUT and FIG_SHARED.
#
# Optional overrides (display or policy, never data):
#   KB_FILE          input file name inside FIG_DATA
#   KB_DISC_MM       disc diameter; readings may not fall below it (default 6)
#   KB_STARS         p (default) | p_adj_within | p_adj_all
#   KB_ABX_ORDER     comma-separated antibiotic order
#   KB_STRAIN_ORDER  comma-separated facet order
#   KB_GROUP_ORDER   comma-separated group order; first level is the open bar

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

INFILE       <- env_or("KB_FILE", "kirby_bauer_in_vitro_mutants.tsv")
DISC_MM      <- as.numeric(env_or("KB_DISC_MM", "6"))
STAR_P       <- env_or("KB_STARS", "p")
TEST         <- env_or("KB_TEST", "welch")
ABX_ORDER    <- env_csv("KB_ABX_ORDER",    c("VA", "LEV", "MEM", "SXT", "MH", "LZD"))
STRAIN_ORDER <- env_csv("KB_STRAIN_ORDER", c("USA300", "SA412", "OCPT02SA01", "OCPT04SA01", "PedCF"))
GROUP_ORDER  <- env_csv("KB_GROUP_ORDER",  c("Ancestor", "Resistant"))

fig_width  <- 12
fig_height <- 13.5

fill_values <- stats::setNames(c("white", "gray20"), GROUP_ORDER)

long <- read_kirby_long(
  file.path(data_dir, INFILE),
  group_col    = "Clone",
  group_levels = GROUP_ORDER,
  abx_order    = ABX_ORDER,
  strain_order = STRAIN_ORDER
)
long <- check_disc_floor(long, DISC_MM)

# Rep 1-5 are treated as independent observations. Whether they are independent
# mutant clones or repeat measurements of one clone changes what the t-test
# means; the Methods need to say which. Logged here so the assumption is not
# silent.
message("Replicate structure assumed independent: ",
        paste(sort(unique(long$Rep)), collapse = ", "))

stats_df <- kirby_welch(long, by_strain = TRUE, test = TEST)
utils::write.table(
  stats_df,
  file.path(results_dir, "kirby_in_vitro_mutants_welch.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)

skipped <- stats_df[!is.na(stats_df$skipped), , drop = FALSE]
if (nrow(skipped)) {
  message("Comparisons not tested (", nrow(skipped), "):")
  for (i in seq_len(nrow(skipped))) {
    message("  ", skipped$Strain[i], " ", skipped$Antibiotic[i], " — ", skipped$skipped[i])
  }
}

sum_df <- summarize_kirby(long, by_strain = TRUE)
ann    <- kirby_brackets(stats_df, long, sum_df, p_col = STAR_P, by_strain = TRUE)

p <- kirby_barplot(
  long, sum_df, ann,
  fill_values     = fill_values,
  facet_by_strain = TRUE,
  facet_ncol      = 2,
  fig_width       = fig_width,
  fig_height      = fig_height,
  subtitle = paste0(
    TEST, " test ", GROUP_ORDER[1], " vs ", GROUP_ORDER[2],
    "; brackets: ", STAR_P, " <= 0.05; n = 5 independent resistant clones per background",
    " (all p-values in kirby_in_vitro_mutants_welch.tsv)"
  )
)

save_panel(p, file.path(results_dir, "kirby_in_vitro_mutants"), fig_width, fig_height)
save_group_legend(fill_values, file.path(results_dir, "kirby_in_vitro_mutants_legend.png"))
