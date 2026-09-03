# Supplementary Table 1 — rate of adsorption, pre- vs post-treatment isolate,
# Welch t-test at each timepoint with a Benjamini-Hochberg column.
#
# Run with ../run.sh, which supplies FIG_DATA, FIG_OUT and FIG_SHARED.
#
# The table holds no statistics of its own: it calls welch_by_timepoint() from
# Figure 4's shared helpers, the same function that places the stars on panel
# 4a. Table and panel therefore cannot disagree.
#
# Read the table with its design in mind. The six timepoints are repeated
# samples of the same three cultures per group, so the rows are not independent
# tests and BH across them is a floor rather than a correction that makes them
# so. The summary that respects the repeated-measures design is the adsorption
# slope comparison in Supplementary Figure 4a.
#
# Optional overrides:
#   T1_FILE   input file name inside FIG_DATA
#   T1_PRE    regex identifying the pre-treatment replicate columns
#   T1_POST   regex identifying the post-treatment replicate columns

need <- function(v) {
  x <- Sys.getenv(v, unset = NA_character_)
  if (is.na(x)) stop(v, " must be set. Run these tables with ./run.sh")
  x
}
env_or <- function(v, default) {
  x <- Sys.getenv(v, unset = "")
  if (!nzchar(x)) default else x
}

data_dir    <- need("FIG_DATA")
results_dir <- need("FIG_OUT")
shared      <- need("FIG_SHARED")
source(shared)

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

INFILE <- env_or("T1_FILE", "roa_pre_post.tsv")
PRE    <- env_or("T1_PRE",  "^OCPT04SA01")
POST   <- env_or("T1_POST", "^OCPT04SA03")

df <- utils::read.delim(file.path(data_dir, INFILE), check.names = FALSE,
                        stringsAsFactors = FALSE)
w  <- welch_by_timepoint(df, pre_pattern = PRE, post_pattern = POST)

out <- data.frame(
  time_min = w$Time_min,
  n1 = w$n_pre, n2 = w$n_post,
  mean1 = w$Pre_mean, sd1 = w$Pre_sd,
  mean2 = w$Post_mean, sd2 = w$Post_sd,
  t = w$t_statistic, df = w$df_Welch, p = w$p_value,
  fdr_bh = w$p_fdr_bh,
  sig_raw = w$significant_p_lt_0.05,
  sig_bh  = w$significant_fdr_lt_0.05,
  stringsAsFactors = FALSE
)

utils::write.table(out, file.path(results_dir, "supp_table01_roa_welch.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

message("Supplementary Table 1: ", nrow(out), " timepoints, ",
        sum(out$sig_raw, na.rm = TRUE), " with p < 0.05, ",
        sum(out$sig_bh, na.rm = TRUE), " surviving BH")
message("Timepoints are repeated measures of the same cultures; see ",
        "Supplementary Figure 4a for the slope comparison.")
