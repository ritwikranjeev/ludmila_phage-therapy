# Figure 4 — shared helpers for panels 4a and 4b.
# Sourced by the panel scripts; sets no paths and produces no output itself.

# Ludmila one-step growth curves (OSGC) and rate-of-attachment (ROA)
# Aesthetics match Figure 4a style: theme_ludmila, mean line + SD error bars + faint replicates.
# One-step growth biological setup (documentation only): MOI 0.001, CFU 10^8/mL — does not change ROA math below.
# Run from repo root: Rscript scripts/osgc.R
# Or: PROJECT_DIR=/path/to/OSGC Rscript scripts/osgc.R

if (!interactive()) pdf(NULL)  # prevent Rplots.pdf when run via Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

env_or <- function(v, default) {
  x <- Sys.getenv(v, unset = "")
  if (!nzchar(x)) default else x
}

# Raster resolution for every PNG in this figure. Panels are assembled in Canva
# and exported to PDF from there, so PNG is the delivered output. Width and
# height are NOT configurable: they match the published originals.
FIG_DPI <- as.numeric(env_or("FIG_DPI", 600))

# ---- figure constants -------------------------------------------------------
# These were referenced by theme_ludmila() (aspect.ratio) and by
# star_annotation_df_attached() (INPUT_PHAGE_ROA) but defined nowhere, so panels
# 4a and 4b both stopped with "object not found". Panel 4c carries its own copy
# of the dimension block, which is why only it still ran.
#
# Input phage for ROA (% free phage = 100 * PFU / INPUT_PHAGE); rate-of-attachment
# assay reference (2e7 PFU/mL). Sets the "% Phage Attached" axis in panel 4a.
INPUT_PHAGE_ROA <- 2e7

# Panel size shared with 4c and fig02/2d; matches the published 8.5x6 in panels.
fig_width   <- 8.5
fig_height  <- 6
plot_margin <- ggplot2::margin(t = 12, r = 14, b = 10, l = 36, unit = "pt")

theme_ludmila <- function() {
  theme_bw(base_size = 13, base_family = "Times New Roman") +
    theme(
      axis.title       = element_text(face = "bold", colour = "black"),
      axis.title.x     = element_text(face = "bold", size = 14, colour = "black", margin = margin(t = 12, unit = "pt")),
      axis.title.y     = element_text(face = "bold", size = 14, colour = "black", margin = margin(r = 12, unit = "pt")),
      axis.text        = element_text(face = "bold", colour = "black"),
      axis.text.x      = element_text(face = "bold", colour = "black", margin = margin(t = 6, unit = "pt")),
      axis.text.y      = element_text(face = "bold", colour = "black", margin = margin(r = 6, unit = "pt")),
      panel.grid.minor = element_blank(),
      aspect.ratio     = fig_height / fig_width,
      panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      axis.ticks       = element_line(colour = "black", linewidth = 0.5),
      axis.line        = element_line(colour = "black", linewidth = 0.5)
    )
}

save_pre_post_legend <- function(p, outfile) {
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    message("Package cowplot not installed; skipped legend file: ", outfile)
    message("Install with: install.packages(\"cowplot\")")
    return(invisible(NULL))
  }
  leg <- cowplot::get_legend(
    p + theme(
      legend.position = "top",
      legend.justification = "center",
      legend.box.margin = margin(10, 14, 10, 14, unit = "pt")
    )
  )
  if (is.null(leg)) return(invisible(NULL))
  pq <- cowplot::plot_grid(leg)
  ggplot2::ggsave(outfile, pq, width = 4.25, height = 1.35, dpi = FIG_DPI, bg = "white")
  message("Wrote legend: ", outfile)
  invisible(outfile)
}

time_col_name <- function(df) names(df)[1L]

grep_cols <- function(df, pattern) {
  grep(pattern, names(df), value = TRUE)
}

summarize_replicates <- function(mat) {
  PFU <- rowMeans(mat, na.rm = TRUE)
  SD  <- apply(mat, 1L, sd, na.rm = TRUE)
  SD[is.na(SD)] <- 0
  list(PFU = PFU, SD = SD)
}

long_two_groups <- function(df, time_nm, cols_a, cols_b, lab_a, lab_b) {
  tvec <- df[[time_nm]]
  n <- nrow(df)
  out <- list()
  for (nm in cols_a) {
    out[[length(out) + 1L]] <- data.frame(
      x = tvec, value = df[[nm]], group = lab_a, stringsAsFactors = FALSE
    )
  }
  for (nm in cols_b) {
    out[[length(out) + 1L]] <- data.frame(
      x = tvec, value = df[[nm]], group = lab_b, stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

# Welch t-test of the two groups at each timepoint.
#
# This is the single implementation behind Figure 4a/4b and Supplementary
# Tables 1 and 2 — the table scripts call this function rather than repeating
# the test, so the panel and the table cannot disagree.
#
# Caution when reading the result: the timepoints come from the same cultures
# sampled repeatedly, so these tests are not independent of one another, and a
# per-timepoint table overstates how much separate evidence it contains. The
# summary that respects the design is the slope comparison in Supplementary
# Figure 4a. BH across timepoints is reported here as a floor, not as a fix.
welch_by_timepoint <- function(df, time_nm = time_col_name(df),
                               pre_pattern = "^OCPT04SA01",
                               post_pattern = "^OCPT04SA03") {
  pre_cols  <- grep_cols(df, pre_pattern)
  post_cols <- grep_cols(df, post_pattern)
  if (!length(pre_cols) || !length(post_cols)) {
    stop("welch_by_timepoint: no columns matching '", pre_pattern, "' and '",
         post_pattern, "' in ", paste(names(df), collapse = ", "))
  }
  rows <- vector("list", nrow(df))
  for (i in seq_len(nrow(df))) {
    x <- as.numeric(df[i, pre_cols, drop = TRUE])
    y <- as.numeric(df[i, post_cols, drop = TRUE])
    x <- x[!is.na(x)]; y <- y[!is.na(y)]
    tt <- if (length(x) > 1L && length(y) > 1L) {
      tryCatch(stats::t.test(x, y, var.equal = FALSE, paired = FALSE),
               error = function(e) NULL)
    } else NULL
    rows[[i]] <- data.frame(
      Time_min  = df[[time_nm]][i],
      n_pre     = length(x),
      n_post    = length(y),
      Pre_mean  = if (length(x)) mean(x) else NA_real_,
      Pre_sd    = if (length(x) > 1L) stats::sd(x) else NA_real_,
      Post_mean = if (length(y)) mean(y) else NA_real_,
      Post_sd   = if (length(y) > 1L) stats::sd(y) else NA_real_,
      t_statistic = if (is.null(tt)) NA_real_ else unname(tt$statistic),
      df_Welch    = if (is.null(tt)) NA_real_ else unname(tt$parameter),
      p_value     = if (is.null(tt)) NA_real_ else tt$p.value,
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  out$p_fdr_bh <- stats::p.adjust(out$p_value, method = "BH")
  out$significant_p_lt_0.05 <- out$p_value < 0.05
  out$significant_fdr_lt_0.05 <- out$p_fdr_bh < 0.05
  out
}

p_to_stars <- function(p) {
  ifelse(
    is.na(p) | p >= 0.05, "",
    ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", "*"))
  )
}

star_annotation_df <- function(welch_df, time_col = "Time_min") {
  stars <- p_to_stars(welch_df$p_value)
  idx <- nzchar(stars)
  if (!any(idx)) return(NULL)
  w <- welch_df[idx, , drop = FALSE]
  stars <- stars[idx]
  pre_sd <- ifelse(is.na(w$Pre_sd), 0, w$Pre_sd)
  post_sd <- ifelse(is.na(w$Post_sd), 0, w$Post_sd)
  y_top <- pmax(w$Pre_mean + pre_sd, w$Post_mean + post_sd)
  data.frame(
    Time = w[[time_col]],
    y = y_top,
    label = stars,
    stringsAsFactors = FALSE
  )
}

bump_star_y <- function(stars_df, summary_df, frac = 0.045) {
  if (is.null(stars_df)) return(NULL)
  y_max <- max(summary_df$PFU + summary_df$SD, na.rm = TRUE)
  stars_df$y <- stars_df$y + y_max * frac
  stars_df
}

star_annotation_df_attached <- function(welch_df, time_col = "Time_min") {
  stars <- p_to_stars(welch_df$p_value)
  idx <- nzchar(stars)
  if (!any(idx)) return(NULL)
  w <- welch_df[idx, , drop = FALSE]
  stars <- stars[idx]
  pre_sd <- ifelse(is.na(w$Pre_sd), 0, w$Pre_sd)
  post_sd <- ifelse(is.na(w$Post_sd), 0, w$Post_sd)
  pre_up <- (100 - 100 * w$Pre_mean / INPUT_PHAGE_ROA) + 100 * pre_sd / INPUT_PHAGE_ROA
  post_up <- (100 - 100 * w$Post_mean / INPUT_PHAGE_ROA) + 100 * post_sd / INPUT_PHAGE_ROA
  y_top <- pmax(pre_up, post_up)
  data.frame(
    Time = w[[time_col]],
    y = y_top,
    label = stars,
    stringsAsFactors = FALSE
  )
}
