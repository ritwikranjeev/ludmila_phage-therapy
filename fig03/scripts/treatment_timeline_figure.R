# Publication-style treatment timeline (ggplot2 + patchwork)
# Data: events.tsv, labs.tsv, courses.tsv (same folder as this script, or under data_dir)

library(tidyverse)
library(lubridate)
library(patchwork)
library(scales)
library(cowplot)

# --- paths -------------------------------------------------------------------
# Supplied by run.sh. Nothing here knows where it sits in the tree.
data_dir    <- Sys.getenv("FIG_DATA", unset = NA_character_)
results_dir <- Sys.getenv("FIG_OUT",  unset = NA_character_)
if (is.na(data_dir) || is.na(results_dir))
  stop("FIG_DATA and FIG_OUT must be set. Run this figure with ./run.sh, or:\n",
       "  FIG_DATA=data FIG_OUT=output Rscript scripts/treatment_timeline_figure.R")
if (!dir.exists(data_dir)) stop("FIG_DATA does not exist: ", data_dir)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

env_or <- function(v, default) {
  x <- Sys.getenv(v, unset = "")
  if (!nzchar(x)) default else x
}

# Raster resolution for every PNG in this figure. Panels are assembled in Canva
# and exported to PDF from there, so PNG is the delivered output. Width and
# height are NOT configurable: they match the published originals.
FIG_DPI <- as.numeric(env_or("FIG_DPI", 600))

# Best-effort vector companion to the PNG.
#
# cairo_pdf() can fail at runtime with a WARNING rather than an error: it writes
# no file, tryCatch(error=) never fires, and the caller would otherwise report a
# PDF that does not exist (see CLAUDE.md, Environment). So attempt the write and
# then confirm a fresh, non-empty file actually landed. Returns TRUE only then.
#
# There is deliberately no base pdf() fallback: that device rejects
# 'Times New Roman' and would reset the panel in a second face.
save_pdf_best_effort <- function(plot, path, ...) {
  t0 <- Sys.time()
  tryCatch(
    suppressWarnings(ggsave(path, plot = plot, device = grDevices::cairo_pdf, ...)),
    error = function(e) invisible(NULL))
  isTRUE(file.exists(path) && file.size(path) > 0 && file.mtime(path) >= t0)
}

# input filenames, overridable without editing this script
F_EVENTS  <- Sys.getenv("TL_EVENTS",  unset = "events_relative_days.tsv")
F_LABS    <- Sys.getenv("TL_LABS",    unset = "labs_relative_days.tsv")
F_COURSES <- Sys.getenv("TL_COURSES", unset = "courses_relative_days.tsv")
F_RANGES  <- Sys.getenv("TL_RANGES",  unset = "normal_ranges.tsv")
for (f in c(F_EVENTS, F_LABS, F_COURSES, F_RANGES))
  if (!file.exists(file.path(data_dir, f)))
    stop("Input not found: ", file.path(data_dir, f),
         "\n  Present in that directory: ", paste(list.files(data_dir), collapse = ", "))

path_tsv <- function(name) file.path(data_dir, name)

# -----------------------------------------------------------------------------
# Read data
# -----------------------------------------------------------------------------
events <- read_tsv(path_tsv(F_EVENTS), show_col_types = FALSE)
labs <- read_tsv(path_tsv(F_LABS), show_col_types = FALSE)
courses <- read_tsv(path_tsv(F_COURSES), show_col_types = FALSE)
norm_ranges_raw <- read_tsv(path_tsv(F_RANGES), show_col_types = FALSE)

# -----------------------------------------------------------------------------
# Study days
#
# Inputs are de-identified: calendar dates have been replaced by study day, with
# day 0 = first intravenous vancomycin dose. No dates or identifiers are present.
# See ../../docs/DEIDENTIFICATION.md.
# -----------------------------------------------------------------------------
events  <- events  %>% mutate(day = as.numeric(day))
labs    <- labs    %>% mutate(day = as.numeric(day)) %>% mutate(Date = day)
courses <- courses %>% mutate(day_start = as.numeric(day_start),
                              day_end   = as.numeric(day_end))
day0 <- 0

# -----------------------------------------------------------------------------
# Requested labs only; parse_numeric strips "<2.5"-style text in other columns
# -----------------------------------------------------------------------------
lab_cols <- c(
  "Monocyte Count",
  "Neutrophil Count",
  "WBC",
  "Creatinine",
  "BUN",
  "Platelet ct",
  "Hemoglobin ct"
)

labs_long_raw <- labs %>%
  select(Date, day, all_of(lab_cols)) %>%
  pivot_longer(
    cols = all_of(lab_cols),
    names_to = "lab",
    values_to = "value_raw"
  ) %>%
  mutate(value = readr::parse_number(as.character(value_raw))) %>%
  filter(!is.na(value))

# Raw long table for supplemental / alternate figures
labs_long_values <- labs_long_raw

# Facet / y-axis labels: analyte spelled out; second line = SI-style unit abbreviations
lab_facet_map <- c(
  "Monocyte Count"   = "Monocyte count\n(cells/\u00b5L)",
  "Neutrophil Count" = "Neutrophil count\n(\u00d7 10\u2079/L)",
  "WBC"              = "White blood cell count\n(\u00d7 10\u2079/L)",
  "Creatinine"       = "Serum creatinine\n(mg/dL)",
  "BUN"              = "Blood urea nitrogen\n(mg/dL)",
  "Platelet ct"      = "Platelet count\n(\u00d7 10\u00b3/\u00b5L)",
  "Hemoglobin ct"    = "Hemoglobin\n(g/dL)"
)

labs_long_plot <- labs_long_raw %>%
  mutate(
    lab_facet = lab_facet_map[lab],
    lab_facet = factor(lab_facet, levels = unname(lab_facet_map[lab_cols]))
  )

# Normal-range bands (per facet): light green, behind traces
lim_col <- names(norm_ranges_raw)[1]
lab_norm_bands <- norm_ranges_raw %>%
  rename(limit_type = all_of(lim_col)) %>%
  filter(limit_type %in% c("Lower_Limit", "Upper_Limit")) %>%
  pivot_longer(-limit_type, names_to = "lab", values_to = "bound") %>%
  mutate(bound = readr::parse_number(as.character(bound))) %>%
  pivot_wider(names_from = limit_type, values_from = bound) %>%
  filter(lab %in% lab_cols) %>%
  mutate(
    ymin = Lower_Limit,
    ymax = Upper_Limit,
    # Plotted neutrophils are ×10⁹/L; reference sheet uses 1500–7700 → scale to SI
    ymin = if_else(lab == "Neutrophil Count" & ymax > 100, ymin / 1000, ymin),
    ymax = if_else(lab == "Neutrophil Count" & ymax > 100, ymax / 1000, ymax)
  ) %>%
  filter(!is.na(ymin), !is.na(ymax), ymax > ymin) %>%
  mutate(
    lab_facet = lab_facet_map[lab],
    lab_facet = factor(lab_facet, levels = levels(labs_long_plot$lab_facet))
  ) %>%
  select(lab_facet, ymin, ymax)

# Optional: z-scored long data for a supplemental normalized figure
labs_long_z <- labs_long_raw %>%
  group_by(lab) %>%
  mutate(
    mu = mean(value, na.rm = TRUE),
    sigma = stats::sd(value, na.rm = TRUE),
    z_value = if_else(sigma > 0, (value - mu) / sigma, 0)
  ) %>%
  ungroup() %>%
  select(-mu, -sigma)

# -----------------------------------------------------------------------------
# Treatment panel: two rows — Antibiotics vs Phage; fill = specific agent
# -----------------------------------------------------------------------------
treatment_colors <- c(
  "Vancomycin IV" = "#B55D5C",
  "Linezolid PO" = "#5F7C99",
  "Ludmila" = "#5D8A6A"
)

# Display names (data column `item` stays "Ludmila" from TSV)
treatment_fill_labels <- c(
  "Vancomycin IV",
  "Linezolid PO",
  "Ludmila IN"
)

courses <- courses %>%
  mutate(
    therapy_group = if_else(stack == "Phage", "Phage", "Antibiotics"),
    therapy_group = factor(
      therapy_group,
      levels = c("Phage", "Antibiotics")
    ),
    item = factor(item, levels = names(treatment_colors))
  )

# -----------------------------------------------------------------------------
# X-axis: day 0–50; 7-day ticks
# -----------------------------------------------------------------------------
xmin <- 0
xmax_plot <- 50

blood_draw_days <- events %>%
  filter(type == "Blood draw") %>%
  pull(day)

blood_draw_days <- blood_draw_days[blood_draw_days >= xmin & blood_draw_days <= xmax_plot]

ticks_every_7 <- function(from, to) {
  seq(7 * ceiling(from / 7), 7 * floor(to / 7), by = 7)
}

x_breaks_all <- unique(c(ticks_every_7(0, xmax_plot - 1), xmax_plot))
x_limits <- c(xmin, xmax_plot)

# Treatment bars: half-height, equal gap between rows and same inset above/below pair
bar_half <- 0.27
inter_bar_gap <- 1 - 2 * bar_half
edge_pad <- inter_bar_gap
y_top_treat <- 2 + bar_half + edge_pad
y_phage_bar_bottom <- 1 - bar_half
y_treat_margin_bottom <- y_phage_bar_bottom - edge_pad

# Red markers sit below that margin (lower than previous offset)
y_blood_draw <- y_treat_margin_bottom - 0.11

blood_draw_markers <- tibble(
  day = blood_draw_days,
  y = y_blood_draw,
  key = "Blood draw"
)

# -----------------------------------------------------------------------------
# Theme: white panel, light y grid, Times New Roman
# -----------------------------------------------------------------------------
font_serif <- "Times New Roman"

# Panels are assembled several to a page in Canva and shrunk again, so default
# type reads far too small in the finished layout. Every size below is a
# multiple of FIG_TEXT_SCALE; figure dimensions are unaffected.
FIG_TEXT_SCALE <- as.numeric(Sys.getenv("FIG_TEXT_SCALE", unset = "2"))
.s <- FIG_TEXT_SCALE

# theme_classic rather than theme_bw: axes only, no gridlines, no panel box.
theme_clinical <- theme_classic(base_size = 11 * .s, base_family = font_serif) +
  theme(
    text = element_text(family = font_serif),
    axis.line = element_line(colour = "black", linewidth = 0.6),
    axis.ticks = element_line(colour = "black"),
    axis.title = element_text(size = 11 * .s, color = "black", family = font_serif),
    axis.text = element_text(size = 10 * .s, color = "black", family = font_serif),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", family = font_serif),
    legend.title = element_text(size = 10 * .s, family = font_serif),
    legend.text = element_text(size = 9 * .s, family = font_serif),
    plot.title = element_text(face = "bold", size = 12 * .s, family = font_serif),
    plot.margin = margin(8, 12, 8, 8)
  )

# -----------------------------------------------------------------------------
# Top: treatment courses (muted fills by agent; black outlines)
# -----------------------------------------------------------------------------
make_p_treat <- function(legend_position = "none") {
  ggplot() +
    geom_rect(
      data = courses,
      aes(
        xmin = day_start,
        xmax = day_end,
        ymin = as.numeric(therapy_group) - bar_half,
        ymax = as.numeric(therapy_group) + bar_half,
        fill = item
      ),
      color = "black",
      linewidth = 0.4,
      alpha = 0.92
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "black",
      linewidth = 0.5
    ) +
    geom_point(
      data = blood_draw_markers,
      aes(x = day, y = y, color = key),
      inherit.aes = FALSE,
      shape = 19,
      size = 3,
      stroke = 0.25
    ) +
    scale_x_continuous(
      limits = x_limits,
      breaks = x_breaks_all,
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    scale_y_continuous(
      breaks = seq_along(levels(courses$therapy_group)),
      labels = levels(courses$therapy_group),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_fill_manual(
      values = treatment_colors,
      breaks = names(treatment_colors),
      labels = treatment_fill_labels,
      name = "Treatment",
      drop = FALSE,
      guide = guide_legend(order = 1)
    ) +
    scale_color_manual(
      name = "Monitoring",
      values = c("Blood draw" = "#e31a1c"),
      breaks = "Blood draw",
      labels = "Blood draw",
      guide = guide_legend(order = 2, override.aes = list(size = 3, shape = 19))
    ) +
    coord_cartesian(
      ylim = c(y_blood_draw - 0.12, y_top_treat),
      clip = "off"
    ) +
    labs(
      title = NULL,
      y = NULL,
      x = NULL
    ) +
    theme_clinical +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = legend_position,
      legend.box = "vertical",
      legend.spacing.y = unit(4, "pt")
    )
}

p_treat <- make_p_treat(legend_position = "none")

# -----------------------------------------------------------------------------
# Stacked panels: one lab per row, raw y scale, shared x from Day 0
# -----------------------------------------------------------------------------
p_labs <- ggplot(labs_long_plot, aes(x = day, y = value, group = lab)) +
  geom_rect(
    data = lab_norm_bands,
    aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "#7ccd7c",
    alpha = 0.28,
    color = NA
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.5
  ) +
  geom_line(linewidth = 0.65, color = "black") +
  geom_point(size = 2, color = "black", shape = 19) +
  facet_wrap(
    vars(lab_facet),
    ncol = 1,
    scales = "free_y",
    strip.position = "left"
  ) +
  scale_x_continuous(
    limits = x_limits,
    breaks = x_breaks_all,
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  coord_cartesian(clip = "on") +
  labs(
    x = "Days from start of intravenous vancomycin (day 0)",
    y = NULL
  ) +
  theme_clinical +
  theme(
    legend.position = "none",
    strip.placement = "outside",
    strip.clip = "off",
    # Sideways strip text; two lines + lineheight limits within-line crowding
    strip.text.y.left = element_text(
      face = "bold",
      size = 10.75 * .s,
      angle = 90,
      hjust = 0.5,
      vjust = 0.5,
      lineheight = 1.12,
      family = font_serif,
      margin = margin(0, 11, 0, 11)
    ),
    panel.spacing.y = unit(18, "pt"),
    plot.margin = margin(8, 10, 8, 6),
    axis.text.y = element_text(family = font_serif, margin = margin(0, 5, 0, 0))
  )

# -----------------------------------------------------------------------------
# Combine
# -----------------------------------------------------------------------------
# Treatment row height = one lab facet row (7 lab panels → weights 1 + 7)
n_lab_panels <- length(lab_cols)

final_plot <- p_treat / p_labs +
  plot_layout(heights = c(1, n_lab_panels)) &
  theme(legend.position = "none")

final_plot <- final_plot +
  plot_annotation(
    title = "Treatment timeline and laboratory trends",
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 13 * .s,
        hjust = 0,
        family = font_serif
      )
    )
  )

# -----------------------------------------------------------------------------
# Save. PNG is the delivered output; the PDF is best-effort.
#
# capabilities("cairo") used to select the device here, but it is a compile-time
# flag and returns TRUE even when cairo cannot load (see CLAUDE.md, Environment),
# so it selected a device that silently wrote nothing. The "pdf" arm of that
# choice is gone too: base pdf() rejects 'Times New Roman'.
# -----------------------------------------------------------------------------
fig_w <- 12.75
fig_h <- 13.5

ggsave(
  filename = file.path(results_dir, "treatment_timeline_figure.png"),
  plot = final_plot,
  width = fig_w,
  height = fig_h,
  units = "in",
  dpi = FIG_DPI
)

if (!save_pdf_best_effort(final_plot,
                          file.path(results_dir, "treatment_timeline_figure.pdf"),
                          width = fig_w, height = fig_h, units = "in")) {
  message("No .pdf for treatment_timeline_figure: cairo_pdf wrote no file")
}

# -----------------------------------------------------------------------------
# Treatment legend only (for figure assembly in layout software)
# -----------------------------------------------------------------------------
legend_grob <- get_legend(
  make_p_treat(legend_position = "right") +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 11 * .s, family = font_serif),
      legend.text = element_text(size = 10 * .s, family = font_serif),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      legend.box.margin = margin(4, 4, 4, 4)
    )
)

p_legend_only <- ggdraw(legend_grob)

legend_w <- 3.4
legend_h <- 2.35

ggsave(
  filename = file.path(results_dir, "treatment_timeline_legend.png"),
  plot = p_legend_only,
  width = legend_w,
  height = legend_h,
  units = "in",
  dpi = FIG_DPI
)

if (!save_pdf_best_effort(p_legend_only,
                          file.path(results_dir, "treatment_timeline_legend.pdf"),
                          width = legend_w, height = legend_h, units = "in")) {
  message("No .pdf for treatment_timeline_legend: cairo_pdf wrote no file")
}

# Avoid Rplots.pdf when running non-interactively (e.g. Rscript)
if (interactive()) {
  print(final_plot)
}
