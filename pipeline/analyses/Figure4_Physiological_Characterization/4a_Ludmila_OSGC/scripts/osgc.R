# Figure 4a: Ludmila one-step growth curve — plot titer vs time with mean and error bars
# Reads project data, computes mean/SD from replicates, saves to results/

if (!interactive()) pdf(NULL)  # prevent Rplots.pdf when run via Rscript

library(ggplot2)

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

data_dir    <- file.path(project_dir, "data", "Figure4_Physiological_Characterization", "4a_Ludmila_OSGC")
results_dir <- file.path(project_dir, "analyses", "Figure4_Physiological_Characterization", "4a_Ludmila_OSGC", "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(data_dir, "4a_Ludmila_OSGC_Input.tsv")
osgc <- read.delim(input_file, check.names = FALSE)

# Titer columns and time column
titer_cols <- grep("^Titer_R[0-9]", names(osgc), value = TRUE)
if (length(titer_cols) < 2L) stop("Need at least 2 titer replicate columns in ", input_file)
time_col <- names(osgc)[1L]

osgc$PFU <- rowMeans(osgc[, titer_cols, drop = FALSE], na.rm = TRUE)
osgc$SD  <- apply(osgc[, titer_cols, drop = FALSE], 1L, sd, na.rm = TRUE)
osgc$SD[is.na(osgc$SD)] <- 0

# Long format for replicate points (faint)
osgc_long <- data.frame(
  x = rep(osgc[[time_col]], length(titer_cols)),
  value = c(as.matrix(osgc[, titer_cols, drop = FALSE]))
)

p <- ggplot(osgc, aes(x = .data[[time_col]], y = PFU)) +
  geom_line(color = "black", linewidth = 0.8) +
  geom_errorbar(aes(ymin = PFU - SD, ymax = PFU + SD), width = 2.5, color = "gray40", linewidth = 0.5) +
  geom_point(data = osgc_long, aes(x = x, y = value), color = "gray30", alpha = 0.5, size = 1.8, inherit.aes = FALSE) +
  geom_point(color = "black", size = 2.5) +
  scale_x_continuous(breaks = osgc[[time_col]]) +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05)),
    labels = function(x) x / 1e7,
    name = expression(bold("Phage Titer (")*bold(10^7)*bold(" PFU/mL)"))
  ) +
  labs(x = "Time (min)", tag = "A") +
  theme_ludmila() +
  theme(
    plot.margin = plot_margin,
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
ggsave(file.path(results_dir, "Ludmila_OSGC.png"), plot = p, width = fig_width, height = fig_height, units = "in", dpi = 300, bg = "white")
