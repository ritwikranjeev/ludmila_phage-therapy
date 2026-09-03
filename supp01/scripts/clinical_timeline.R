# Supplementary Figure 1 — clinical timeline of MRSA sinusitis and adjunctive
# phage therapy, from first presentation through post-treatment follow-up.
#
# Run with ../run.sh, which supplies FIG_DATA and FIG_OUT.
#
# Time is expressed in days relative to the first intravenous vancomycin dose
# (day 0), the same anchor Figure 3 uses. No calendar dates exist anywhere in
# this figure or its input; see the note on the input file below.
#
# Optional overrides (display or policy, never data):
#   S1_FILE          input file name inside FIG_DATA
#   S1_SPAN_DAYS     width drawn for an entry with no recorded end (default 14)
#   S1_HIDE_LABELS   comma-separated labels to draw as bars/points but not text
#   S1_BREAK_MONTHS  x-axis tick spacing in months (default 6)
#   S1_XLIM_DAYS     "min,max" in days; default is the data range plus padding
#   S1_TITLE         plot title; set empty to omit

if (!interactive()) pdf(NULL)

suppressPackageStartupMessages({
  library(ggplot2)
})

# Typography. The set is standardised on Times New Roman; "Times" and "serif"
# resolve to a DIFFERENT file on macOS (Times.ttc, the Linotype face), so mixing
# the names silently sets panels in two faces. Fall back only if the font is
# genuinely absent, e.g. on a Linux reproducibility run.
FONT_FAMILY <- "Times New Roman"
if (requireNamespace("systemfonts", quietly = TRUE)) {
  if (!FONT_FAMILY %in% unique(systemfonts::system_fonts()$family)) FONT_FAMILY <- "Times"
} else {
  FONT_FAMILY <- "Times"
}

need <- function(v) {
  x <- Sys.getenv(v, unset = NA_character_)
  if (is.na(x)) stop(v, " must be set. Run this figure with ./run.sh")
  x
}
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
env_csv <- function(v, default = NULL) {
  x <- Sys.getenv(v, unset = "")
  if (!nzchar(x)) return(default)
  trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
}

data_dir    <- need("FIG_DATA")
results_dir <- need("FIG_OUT")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

INFILE       <- env_or("S1_FILE", "clinical_timeline_relative_days.tsv")
SPAN_DAYS    <- as.numeric(env_or("S1_SPAN_DAYS", "14"))
HIDE_LABELS  <- env_csv("S1_HIDE_LABELS", c("Fevers and Sinus Pressure", "Endoscopic Debriment"))
BREAK_MONTHS <- as.numeric(env_or("S1_BREAK_MONTHS", "6"))
XLIM_DAYS    <- as.numeric(env_csv("S1_XLIM_DAYS", character(0)))
TITLE        <- env_or("S1_TITLE",
                       "Clinical timeline of MRSA sinusitis and adjunctive phage therapy")

DAYS_PER_MONTH <- 30.4375   # mean Gregorian month, for tick placement only

fig_width  <- 8.5
fig_height <- 5.2

# ---- read ------------------------------------------------------------------

df <- utils::read.delim(file.path(data_dir, INFILE), check.names = FALSE,
                        stringsAsFactors = FALSE)
for (nm in c("track", "type", "day_start", "day_end", "label", "color_group")) {
  if (!nm %in% names(df)) stop("Column '", nm, "' missing from ", INFILE)
}
df$day_start <- as.numeric(df$day_start)
df$day_end   <- as.numeric(df$day_end)

blank_to_na <- function(x) { x <- trimws(x); x[x == "" | x == "NA"] <- NA; x }
df$label       <- blank_to_na(df$label)
df$group       <- if ("group" %in% names(df)) blank_to_na(df$group) else NA_character_
df$color_group <- blank_to_na(df$color_group)

# ---- normalise grouping ----------------------------------------------------

canon <- function(x) {
  x <- tolower(x)
  x[x == "debriment"] <- "surgery"
  x
}
df$group       <- canon(df$group)
df$color_group <- canon(df$color_group)
df$group       <- ifelse(is.na(df$group), df$type, df$group)
df$color_group <- ifelse(is.na(df$color_group), df$group, df$color_group)

# Labels that are drawn as bars or points but carry no text in the figure.
df$label[df$label %in% HIDE_LABELS] <- NA_character_

# ---- spans -----------------------------------------------------------------

# An entry with no recorded end is drawn SPAN_DAYS wide. That width is an
# assumption, not a measurement: the source records only a start for these.
# Anything relying on bar width for these entries is reading a placeholder.
spanned <- df$type %in% c("drug", "phage", "symptom_band")
assumed <- spanned & is.na(df$day_end)
df$day_end[assumed] <- df$day_start[assumed] + SPAN_DAYS
df$day_end[is.na(df$day_end)] <- df$day_start[is.na(df$day_end)]
df$span_source <- ifelse(assumed, "assumed", "recorded")

if (any(assumed)) {
  message("Entries drawn with an assumed ", SPAN_DAYS, "-day span (", sum(assumed), "):")
  for (i in which(assumed)) {
    message("  day ", df$day_start[i], "  ", df$type[i], "  ",
            if (is.na(df$label[i])) df$color_group[i] else df$label[i])
  }
}

# ---- layout ----------------------------------------------------------------

df$plot_track <- ifelse(df$track == "Disease", "Disease", "Combined")
df$plot_track <- factor(df$plot_track, levels = c("Disease", "Combined"))

# Vertical position within each track. Layout only.
y_disease  <- c(culture = 0.12, symptom_band = 0.32)
y_combined <- c(imaging = 0.82, surgery = 0.82, pathology = 0.82,
                phage = 0.64, linezolid = 0.42, vancomycin = 0.28,
                mupirocin = 0.14, doxycycline = 0.00)

df$y_plot <- ifelse(
  df$plot_track == "Disease",
  unname(y_disease[df$type]),
  unname(y_combined[df$color_group])
)
df$y_plot[is.na(df$y_plot)] <- 0

pal <- c(
  vancomycin  = "#1f77b4",
  linezolid   = "#9467bd",
  doxycycline = "#17becf",
  mupirocin   = "#ff7f0e",
  phage       = "#2ca02c",
  culture     = "#d62728",
  symptom     = "#ffcc80",
  imaging     = "#6c757d",
  surgery     = "#8c564b",
  pathology   = "#6c757d"
)
unknown <- setdiff(unique(df$color_group), names(pal))
if (length(unknown)) {
  message("No palette entry for: ", paste(unknown, collapse = ", "), " — drawn grey")
}

shape_map <- c(culture = 16, imaging = 18, pathology = 18, surgery = 8)

# ---- geometry --------------------------------------------------------------

df$y_base <- as.numeric(df$plot_track) + df$y_plot

rectangles <- df[spanned & !is.na(df$day_start) & !is.na(df$day_end), , drop = FALSE]
rectangles$xmin <- rectangles$day_start
rectangles$xmax <- rectangles$day_end
rectangles$ymin <- rectangles$y_base - 0.035
rectangles$ymax <- rectangles$y_base + 0.035

points <- df[df$type %in% c("culture", "imaging", "pathology", "surgery") &
               !is.na(df$day_start), , drop = FALSE]

rect_labels <- rectangles[!is.na(rectangles$label), , drop = FALSE]
left_labelled <- rect_labels$color_group %in% c("vancomycin", "mupirocin", "doxycycline")
rect_labels$label_x   <- ifelse(left_labelled, rect_labels$xmin - 8, rect_labels$xmax + 8)
rect_labels$hjust_val <- ifelse(left_labelled, 1, 0)
rect_labels$label_y   <- (rect_labels$ymin + rect_labels$ymax) / 2

point_labels <- points[!is.na(points$label), , drop = FALSE]

# ---- axis ------------------------------------------------------------------

if (length(XLIM_DAYS) == 2L) {
  x_limits <- XLIM_DAYS
} else {
  rng <- range(c(df$day_start, df$day_end), na.rm = TRUE)
  pad <- diff(rng) * 0.06
  x_limits <- c(rng[1] - pad, rng[2] + pad)
}

step <- BREAK_MONTHS * DAYS_PER_MONTH
month_breaks <- seq(
  from = ceiling(x_limits[1] / step) * step,
  to   = floor(x_limits[2] / step) * step,
  by   = step
)
month_labels <- sprintf("%+d", round(month_breaks / DAYS_PER_MONTH))
month_labels[round(month_breaks) == 0] <- "0"

# ---- plot ------------------------------------------------------------------

p <- ggplot() +
  geom_rect(
    data = rectangles,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = color_group),
    colour = "black", linewidth = 0.25
  ) +
  geom_point(
    data = points,
    aes(x = day_start, y = y_base, colour = color_group, shape = type),
    size = 2.8, stroke = 0.9
  ) +
  geom_text(
    data = rect_labels,
    aes(x = label_x, y = label_y, label = label, hjust = hjust_val),
    family = FONT_FAMILY, size = 2.9, colour = "black"
  )

if (nrow(point_labels)) {
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      data = point_labels,
      aes(x = day_start, y = y_base, label = label),
      family = FONT_FAMILY, size = 2.8, max.overlaps = 50,
      box.padding = 0.2, point.padding = 0.15,
      segment.color = "grey50", min.segment.length = 0
    )
  } else {
    message("Package ggrepel not installed; point labels drawn without repulsion")
    p <- p + geom_text(
      data = point_labels,
      aes(x = day_start, y = y_base, label = label),
      family = FONT_FAMILY, size = 2.8, vjust = -1
    )
  }
}

p <- p +
  scale_fill_manual(values = pal, na.value = "grey80") +
  scale_colour_manual(values = pal, na.value = "black") +
  scale_shape_manual(values = shape_map) +
  scale_x_continuous(
    limits = x_limits,
    breaks = month_breaks,
    labels = month_labels,
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  scale_y_continuous(
    breaks = c(1, 2), labels = c("", ""),
    expand = expansion(mult = c(0.12, 0.30))
  ) +
  labs(
    title = if (nzchar(TITLE)) TITLE else NULL,
    x = "Months relative to first intravenous vancomycin dose (day 0)",
    y = NULL,
    fill = "Treatment / event",
    colour = "Point event",
    shape = "Event type"
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 10, base_family = FONT_FAMILY) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    axis.text.x        = element_text(size = 8.5, colour = "black", family = FONT_FAMILY),
    axis.title.x       = element_text(size = 9, colour = "black", family = FONT_FAMILY,
                                      margin = margin(t = 8, unit = "pt")),
    legend.position    = "bottom",
    legend.box         = "vertical",
    legend.title       = element_text(size = 9, face = "bold", family = FONT_FAMILY),
    legend.text        = element_text(size = 8, family = FONT_FAMILY),
    plot.title         = element_text(face = "bold", size = 11, hjust = 0.5, family = FONT_FAMILY),
    plot.margin        = margin(10, 40, 10, 40)
  )

out_base <- file.path(results_dir, "clinical_timeline")
ggsave(paste0(out_base, ".png"), plot = p, width = fig_width, height = fig_height,
       units = "in", dpi = FIG_DPI, bg = "white")
if (save_pdf_best_effort(p, paste0(out_base, ".pdf"), width = fig_width,
                         height = fig_height, units = "in", bg = "white")) {
  message("Wrote: ", out_base, ".png and .pdf")
} else {
  message("Wrote: ", out_base, ".png   (no .pdf: cairo_pdf wrote no file)")
}
