# Figure 4a — rate of attachment: percent of phage adsorbed over time,
# pre-treatment vs post-treatment isolate.
#
# NOTE: the same data can be rendered as "% free phage" (the falling curve).
# That is NOT the published panel. This script produces "% phage attached".

# ---- paths -------------------------------------------------------------------
# Supplied by run.sh. This script knows nothing about the directory layout.
need <- function(v) {
  x <- Sys.getenv(v, unset = NA_character_)
  if (is.na(x)) stop(v, " must be set. Run this figure with ./run.sh")
  x
}
data_dir    <- need("FIG_DATA")
results_dir <- need("FIG_OUT")
shared      <- need("FIG_SHARED")
if (!dir.exists(data_dir)) stop("FIG_DATA does not exist: ", data_dir)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
source(shared)

plot_roa_two_attached <- function(df, outfile, title = NULL) {
  tc <- time_col_name(df)
  pre_cols  <- grep_cols(df, "^OCPT04SA01")
  post_cols <- grep_cols(df, "^OCPT04SA03")
  pre_s  <- summarize_replicates(df[, pre_cols, drop = FALSE])
  post_s <- summarize_replicates(df[, post_cols, drop = FALSE])
  pct_pre  <- 100 - 100 * pre_s$PFU / INPUT_PHAGE_ROA
  pct_post <- 100 - 100 * post_s$PFU / INPUT_PHAGE_ROA
  sd_pre  <- 100 * pre_s$SD / INPUT_PHAGE_ROA
  sd_post <- 100 * post_s$SD / INPUT_PHAGE_ROA

  summary_df <- rbind(
    data.frame(Time = df[[tc]], PFU = pct_pre, SD = sd_pre, group = "Pre"),
    data.frame(Time = df[[tc]], PFU = pct_post, SD = sd_post, group = "Post")
  )
  summary_df$group <- factor(summary_df$group, levels = c("Pre", "Post"))

  long_df <- long_two_groups(df, tc, pre_cols, post_cols, "Pre", "Post")
  long_df$value <- 100 - 100 * long_df$value / INPUT_PHAGE_ROA
  long_df$group <- factor(long_df$group, levels = c("Pre", "Post"))

  tuniq <- sort(unique(df[[tc]]))
  dt_min <- if (length(tuniq) > 1L) min(diff(tuniq)) else 1
  pd <- ggplot2::position_dodge(width = dt_min * 0.12)

  welch_tab <- welch_by_timepoint(df, tc)
  stars_df <- bump_star_y(star_annotation_df_attached(welch_tab, "Time_min"), summary_df)

  p <- ggplot(summary_df, aes(x = Time, y = PFU, color = group, shape = group, group = group)) +
    geom_line(linewidth = 0.8, position = pd) +
    geom_errorbar(
      aes(ymin = PFU - SD, ymax = PFU + SD),
      width = dt_min * 0.25,
      linewidth = 0.5,
      position = pd
    ) +
    geom_point(
      data = long_df,
      aes(x = x, y = value, color = group, shape = group),
      alpha = 0.5,
      size = 1.8,
      inherit.aes = FALSE,
      position = pd
    ) +
    geom_point(aes(color = group, shape = group), size = 2.5, position = pd) +
    {
      if (!is.null(stars_df)) {
        geom_text(
          data = stars_df,
          aes(x = Time, y = y, label = label),
          inherit.aes = FALSE,
          show.legend = FALSE,
          family = "Times New Roman",
          fontface = "bold",
          size = 5.2 * FIG_TEXT_SCALE,
          color = "black"
        )
      } else {
        NULL
      }
    } +
    scale_color_manual(values = c(Pre = "black", Post = "gray40"), name = NULL) +
    scale_shape_manual(values = c(Pre = 16, Post = 15), name = NULL) +
    scale_x_continuous(breaks = df[[tc]]) +
    scale_y_continuous(
      limits = c(0, NA),
      expand = expansion(mult = c(0, if (!is.null(stars_df)) 0.12 else 0.05)),
      name = "% Phage Attached"
    ) +
    labs(x = "Time (min)", title = title) +
    theme_ludmila() +
    theme(
      plot.margin = plot_margin,
      legend.position = "top",
      legend.justification = "left"
    ) +
    coord_cartesian(clip = "off")

  p <- p + theme(
    axis.title.y      = element_text(face = "bold", colour = "black"),
    axis.title.y.left = element_text(face = "bold", colour = "black"),
    axis.text.y       = element_text(face = "bold", colour = "black"),
    axis.text.y.left  = element_text(face = "bold", colour = "black"),
    axis.text.x       = element_text(face = "bold", colour = "black")
  )
  p_main <- p + theme(legend.position = "none")
  ggsave(outfile, plot = p_main, width = fig_width, height = fig_height, units = "in", dpi = FIG_DPI, bg = "white")
  invisible(p)
}

roa <- read.delim(file.path(data_dir, "roa_pre_post.tsv"), check.names = FALSE)

plot_roa_two_attached(roa, file.path(results_dir, "ROA_pre_post_attached.png"))

welch <- welch_by_timepoint(roa)
write.csv(welch, file.path(results_dir, "ROA_pre_post_welch_ttest.csv"), row.names = FALSE)

message("Figure 4a -> ", file.path(results_dir, "ROA_pre_post_attached.png"))
print(welch, digits = 4)
