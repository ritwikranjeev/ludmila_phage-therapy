# Supplementary Figure 2b — phage titre after one hour at pH 2 to 12.
#
# Run with ../../run.sh, which supplies FIG_DATA, FIG_OUT and FIG_SHARED.
#
# Optional overrides:
#   S2_PH_FILE   input file name inside FIG_DATA

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

INFILE <- env_or("S2_PH_FILE", "ph_stability.tsv")

fig_width  <- 8.5
fig_height <- 6

d   <- read_titer_wide(file.path(data_dir, INFILE))
sum_df <- titer_summary(d)

message("Conditions: ", paste(d$x, collapse = ", "),
        " | replicates: ", length(d$titer_cols))

write_stats(anova_tukey(d), results_dir, "ph_stability")

p <- ggplot(sum_df, aes(x = x, y = mean)) +
  geom_line(colour = "black", linewidth = 0.8) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                width = 0.35, colour = "gray40", linewidth = 0.5) +
  geom_point(data = d$long, aes(x = x, y = titer),
             colour = "gray30", alpha = 0.5, size = 1.8, inherit.aes = FALSE) +
  geom_point(colour = "black", size = 2.5) +
  scale_x_continuous(breaks = d$x) +
  titer_y_scale() +
  labs(x = "pH",
       subtitle = if (FIG_SUBTITLE)
         sprintf("One-way ANOVA across pH; error bars are SD of %d technical replicates",
                 length(d$titer_cols)) else NULL) +
  theme_ludmila(fig_width, fig_height) +
  theme(plot.subtitle = element_text(size = 9, hjust = 0,
                                     margin = margin(b = 8, unit = "pt"))) +
  theme(plot.margin = margin(t = 12, r = 14, b = 10, l = 36, unit = "pt")) +
  coord_cartesian(clip = "off")

save_panel(p, file.path(results_dir, "ph_stability"), fig_width, fig_height)
