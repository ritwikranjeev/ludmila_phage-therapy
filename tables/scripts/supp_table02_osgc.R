# Supplementary Table 2 — one-step growth, pre- vs post-treatment isolate,
# Welch t-test at each timepoint with a Benjamini-Hochberg column.
#
# Run with ../run.sh, which supplies FIG_DATA, FIG_OUT and FIG_SHARED.
#
# Same function as Table 1 and as the stars on panel 4b, so the three cannot
# disagree. The same caution applies: the nine timepoints are repeated samples
# of the same cultures, not nine independent experiments.
#
# Optional overrides:
#   T2_FILE   input file name inside FIG_DATA
#   T2_PRE    regex identifying the pre-treatment replicate columns
#   T2_POST   regex identifying the post-treatment replicate columns

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

INFILE <- env_or("T2_FILE", "osgc_post-pre-raw.tsv")
PRE    <- env_or("T2_PRE",  "^OCPT04SA01")
POST   <- env_or("T2_POST", "^OCPT04SA03")

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

utils::write.table(out, file.path(results_dir, "supp_table02_osgc_welch.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

message("Supplementary Table 2: ", nrow(out), " timepoints, ",
        sum(out$sig_raw, na.rm = TRUE), " with p < 0.05, ",
        sum(out$sig_bh, na.rm = TRUE), " surviving BH")

# The burst-size claim rests on the plateau, which this curve never reaches —
# see Supplementary Figure 4b. Flagged here because a reader lands on this
# table looking for the numbers behind that claim.
last <- nrow(out)
if (last >= 3L) {
  tail_rows <- seq(max(1L, last - 2L), last)
  rising <- all(diff(out$mean1[tail_rows]) > 0) && all(diff(out$mean2[tail_rows]) > 0)
  if (rising) {
    message("Both curves are still rising at the last timepoint (",
            out$time_min[last], " min): burst sizes derived from these data are ",
            "lower bounds.")
  }
}
