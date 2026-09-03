# Supplementary Figure 2c — phage titre after one hour at 4, 25, 37 and 55 C.
#
# Temperature is drawn as a categorical axis, as in the published panel: the
# four temperatures are not evenly spaced, and a line between them would imply
# titres at temperatures that were never tested.
#
# Run with ../../run.sh, which supplies FIG_DATA, FIG_OUT and FIG_SHARED.
#
# Optional overrides:
#   S2_TEMP_FILE   input file name inside FIG_DATA

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

INFILE <- env_or("S2_TEMP_FILE", "thermal_stability.tsv")

fig_width  <- 8.5
fig_height <- 6

d      <- read_titer_wide(file.path(data_dir, INFILE))
sum_df <- titer_summary(d)

message("Conditions: ", paste(d$x, collapse = ", "),
        " | replicates: ", length(d$titer_cols))

write_stats(anova_tukey(d), results_dir, "thermal_stability")

lvl <- as.character(d$x)
sum_df$grp   <- factor(as.character(sum_df$x), levels = lvl)
d$long$grp   <- factor(as.character(d$long$x), levels = lvl)

p <- ggplot(sum_df, aes(x = grp, y = mean)) +
  geom_col(width = 0.6, fill = "#48494B", colour = NA) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                width = 0.2, colour = "black", linewidth = 0.5) +
  geom_point(data = d$long, aes(x = grp, y = titer),
             colour = "gray30", alpha = 0.5, size = 1.8,
             position = position_jitter(width = 0.12, seed = 1L),
             inherit.aes = FALSE) +
  titer_y_scale() +
  labs(x = "Temperature (°C)",
       subtitle = if (FIG_SUBTITLE)
         sprintf("One-way ANOVA across temperature; error bars are SD of %d technical replicates",
                 length(d$titer_cols)) else NULL) +
  theme_ludmila(fig_width, fig_height) +
  theme(plot.subtitle = element_text(size = 9, hjust = 0,
                                     margin = margin(b = 8, unit = "pt"))) +
  theme(
    legend.position = "none",
    plot.margin     = margin(t = 12, r = 14, b = 10, l = 36, unit = "pt"),
    axis.line.x     = element_blank()
  ) +
  coord_cartesian(clip = "off")

save_panel(p, file.path(results_dir, "thermal_stability"), fig_width, fig_height)
