# Figure 3b: Ludmila temperature stability — plot titer vs temperature with mean and error bars
# Reads project data, computes mean/SD from replicates, optional ANOVA/Tukey, saves to results/

if (!interactive()) pdf(NULL)  # prevent Rplots.pdf when run via Rscript

library(ggplot2)
library(dplyr)

# Shared figure aesthetics (consistent across Figure 3 & 4 scripts)
# Same dimensions and margins so the plot box (panel) is consistent across figures.
fig_width   <- 8.5
fig_height  <- 6
plot_margin <- margin(t = 12, r = 14, b = 10, l = 36, unit = "pt")
theme_ludmila <- function() {
  theme_bw(base_size = 13, base_family = "Times New Roman") +
  theme(
    axis.title       = element_text(face = "bold"),
    axis.title.x     = element_text(face = "bold", size = 14, margin = margin(t = 12, unit = "pt")),
    axis.title.y     = element_text(face = "bold", size = 14, margin = margin(r = 12, unit = "pt")),
    axis.text        = element_text(face = "bold"),
    axis.text.x      = element_text(face = "bold", margin = margin(t = 6, unit = "pt")),
    axis.text.y      = element_text(face = "bold", margin = margin(r = 6, unit = "pt")),
    panel.grid.minor = element_blank(),
    aspect.ratio     = fig_height / fig_width  # same panel proportions in all figures
  )
}

# Project root: from PROJECT_DIR env, or from script path when run via Rscript
project_dir <- Sys.getenv("PROJECT_DIR", "")
if (!nzchar(project_dir)) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[grepl("^--file=", args)]
  if (length(file_arg)) {
    script_path <- sub("^--file=", "", file_arg[1L])
    project_dir <- normalizePath(file.path(dirname(script_path), "..", "..", "..", ".."), mustWork = FALSE)
  } else {
    project_dir <- normalizePath(file.path(getwd(), "..", "..", "..", ".."), mustWork = FALSE)
  }
}
project_dir <- normalizePath(project_dir, mustWork = FALSE)

data_dir    <- file.path(project_dir, "data", "Figure3_Stability_Testing", "3b_Temp_Stability")
results_dir <- file.path(project_dir, "analyses", "Figure3_Stability_Testing", "3b_Temp_Stability", "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(data_dir, "3b_Ludmila_Temp_Data.tsv")
temp_wide <- read.delim(input_file, check.names = FALSE)

# Titer columns and x (temperature) column
titer_cols <- grep("^Titer_R[0-9]", names(temp_wide), value = TRUE)
if (length(titer_cols) < 2L) stop("Need at least 2 titer replicate columns in ", input_file)
x_col <- names(temp_wide)[1L]

# Long form for ANOVA (one row per replicate)
raw <- data.frame(
  Temp = rep(temp_wide[[x_col]], each = length(titer_cols)),
  PFU  = c(t(as.matrix(temp_wide[, titer_cols, drop = FALSE])))
)
raw$Temp <- factor(as.character(raw$Temp), levels = unique(as.character(temp_wide[[x_col]])))

# Summary stats: mean and SD per temperature
summary_stats <- raw %>%
  group_by(Temp) %>%
  summarise(
    mean_PFU = mean(PFU, na.rm = TRUE),
    sd_PFU   = sd(PFU, na.rm = TRUE),
    .groups = "drop"
  )
summary_stats$sd_PFU[is.na(summary_stats$sd_PFU)] <- 0

# Optional: ANOVA and Tukey HSD
anova_result <- aov(PFU ~ Temp, data = raw)
summary(anova_result)
TukeyHSD(anova_result)

# Plot — black bars, no border; error bars and replicate points overlay on top
p <- ggplot(summary_stats, aes(x = Temp, y = mean_PFU)) +
  geom_bar(stat = "identity", width = 0.6, fill = "#48494B", color = NA) +
  geom_errorbar(aes(ymin = mean_PFU - sd_PFU, ymax = mean_PFU + sd_PFU), width = 0.2, color = "black", linewidth = 0.5) +
  geom_point(data = raw, aes(x = Temp, y = PFU), color = "gray30", alpha = 0.5, size = 1.8, position = position_jitter(width = 0.12, seed = 1L), inherit.aes = FALSE) +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05)),
    labels = function(x) x / 1e8,
    name = expression(bold("Phage Titer (")*bold(10^8)*bold(" PFU/mL)"))
  ) +
  labs(x = "Temperature (°C)", tag = "B") +
  theme_ludmila() +
  theme(
    legend.position = "none",
    plot.margin = plot_margin,
    axis.line.x = element_blank(),
    plot.tag = element_text(size = 18, face = "bold", margin = margin(r = 12, unit = "pt")),
    plot.tag.position = "topleft"
  ) +
  coord_cartesian(clip = "off")

# Final theme layer so y-axis bold applies (Times New Roman; include .left for ggplot2 4.x)
p <- p + theme(
  axis.title.y      = element_text(face = "bold"),
  axis.title.y.left = element_text(face = "bold"),
  axis.text.y       = element_text(face = "bold"),
  axis.text.y.left  = element_text(face = "bold")
)
if (interactive()) print(p)
ggsave(file.path(results_dir, "Ludmila_temp.png"), plot = p, width = fig_width, height = fig_height, units = "in", dpi = 300, bg = "white")
