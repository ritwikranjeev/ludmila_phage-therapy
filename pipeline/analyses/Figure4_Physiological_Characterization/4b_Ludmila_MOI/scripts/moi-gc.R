# Figure 4b: Ludmila MOI curves — OD600 over time by MOI, mean and SD ribbons
# Reads project data, averages wells per MOI/condition, saves plots to results/

if (!interactive()) pdf(NULL)  # prevent Rplots.pdf when run via Rscript

library(tidyverse)

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

data_dir    <- file.path(project_dir, "data", "Figure4_Physiological_Characterization", "4b_Ludmila_MOI")
results_dir <- file.path(project_dir, "analyses", "Figure4_Physiological_Characterization", "4b_Ludmila_MOI", "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(data_dir, "4b_Ludmila_MOI_Data.txt")
df <- read.delim(input_file, check.names = FALSE)

# Fix time format (e.g. "0:00:00" -> hours)
df <- df %>%
  mutate(Time = lubridate::hms(Time),
         Hours = as.numeric(Time) / 3600)

# Function to get average and SD of OD600 for each MOI and condition
get_avg_data <- function(df, cols, condition_name) {
  mois <- LETTERS[1:8]
  map_dfr(1:8, function(i) {
    row <- mois[i]
    if (row == "H") {
      if (identical(cols, 1:3)) wells <- c("H1", "H2")
      else if (identical(cols, 4:6)) wells <- c("H4", "H5")
    } else {
      wells <- paste0(row, cols)
    }
    df %>%
      select(Hours, all_of(wells)) %>%
      rowwise() %>%
      mutate(
        avg = mean(c_across(all_of(wells)), na.rm = TRUE),
        sd = sd(c_across(all_of(wells)), na.rm = TRUE),
        MOI = row,
        Condition = condition_name
      ) %>%
      select(Hours, avg, sd, MOI, Condition)
  })
}

ocpt_data   <- get_avg_data(df, 1:3, "OCPT04SA01 vs Ludmila")
sa412_data  <- get_avg_data(df, 4:6, "SA412 vs Ludmila")
blank_data  <- df %>%
  select(Hours, H3, H6) %>%
  rowwise() %>%
  mutate(
    avg = mean(c_across(c(H3, H6)), na.rm = TRUE),
    sd = sd(c_across(c(H3, H6)), na.rm = TRUE),
    MOI = "Blank",
    Condition = "Blank"
  ) %>%
  select(Hours, avg, sd, MOI, Condition)

plot_data <- bind_rows(ocpt_data, sa412_data, blank_data)
moi_values <- c("1", "0.1", "0.01", "0.001", "0.0001", "0.00001", "0.000001", "Bacteria Only", "Blank")
plot_data$MOI <- factor(plot_data$MOI, levels = c(LETTERS[1:8], "Blank"), labels = moi_values)
plot_data <- filter(plot_data, MOI %in% c("1", "0.1", "0.01", "0.001", "Bacteria Only", "Blank") & Hours <= 16)

# Normalize to blank: subtract blank mean at each time point (blank becomes 0)
blank_by_time <- plot_data %>% filter(MOI == "Blank") %>% select(Hours, blank_avg = avg)
plot_data <- plot_data %>%
  left_join(blank_by_time, by = "Hours") %>%
  mutate(avg = avg - blank_avg) %>%
  select(-blank_avg)

palette <- c(
  "1" = "#D55E00", "0.1" = "#0072B2", "0.01" = "#009E73", "0.001" = "#CC79A7",
  "Bacteria Only" = "#A6761D", "Blank" = "#666666"
)

d1 <- filter(plot_data, Condition == "OCPT04SA01 vs Ludmila")
d2 <- filter(plot_data, Condition == "SA412 vs Ludmila")
# Plot 1: OCPT04SA01 vs Ludmila
p1 <- ggplot(d1, aes(x = Hours, y = avg, color = MOI)) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  geom_line(linewidth = 0.8) +
  geom_ribbon(aes(ymin = avg - sd, ymax = avg + sd, fill = MOI), alpha = 0.2, color = NA) +
  labs(x = "Time (hours)", y = "Absorbance (OD600)", tag = "B") +
  scale_x_continuous(limits = c(0, 16), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_ludmila() +
  theme(
    legend.position = "none",
    plot.margin = plot_margin,
    plot.tag = element_text(size = 18, face = "bold", margin = margin(r = 12, unit = "pt")),
    plot.tag.position = "topleft"
  ) +
  coord_cartesian(clip = "off")
# Final theme layer so y-axis bold applies (Times New Roman; include .left for ggplot2 4.x)
p1 <- p1 + theme(
  axis.title.y      = element_text(face = "bold"),
  axis.title.y.left = element_text(face = "bold"),
  axis.text.y       = element_text(face = "bold"),
  axis.text.y.left  = element_text(face = "bold")
)
ggsave(file.path(results_dir, "Ludmila_MOI_OCPT04SA01.png"), plot = p1, width = fig_width, height = fig_height, units = "in", dpi = 300, bg = "white")

# Plot 2: SA412 vs Ludmila
p2 <- ggplot(d2, aes(x = Hours, y = avg, color = MOI)) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  geom_line(linewidth = 0.8) +
  geom_ribbon(aes(ymin = avg - sd, ymax = avg + sd, fill = MOI), alpha = 0.2, color = NA) +
  labs(x = "Time (hours)", y = "Absorbance (OD600)", tag = "B") +
  scale_x_continuous(limits = c(0, 16), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_ludmila() +
  theme(
    legend.position = "none",
    plot.margin = plot_margin,
    plot.tag = element_text(size = 18, face = "bold", margin = margin(r = 12, unit = "pt")),
    plot.tag.position = "topleft"
  ) +
  coord_cartesian(clip = "off")
# Final theme layer so y-axis bold applies (Times New Roman; include .left for ggplot2 4.x)
p2 <- p2 + theme(
  axis.title.y      = element_text(face = "bold"),
  axis.title.y.left = element_text(face = "bold"),
  axis.text.y       = element_text(face = "bold"),
  axis.text.y.left  = element_text(face = "bold")
)
ggsave(file.path(results_dir, "Ludmila_MOI_SA412.png"), plot = p2, width = fig_width, height = fig_height, units = "in", dpi = 300, bg = "white")

if (interactive()) { print(p1); print(p2) }
