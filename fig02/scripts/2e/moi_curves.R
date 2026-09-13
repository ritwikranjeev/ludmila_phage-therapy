# Figure 2e — OD600 over time by multiplicity of infection.
#
# Well assignments come from plate_layout.tsv. Nothing about the plate is
# hardcoded here: to change which wells, strains or MOIs are used, edit the sheet.

if (!interactive()) pdf(NULL)   # no stray Rplots.pdf under Rscript

library(tidyverse)

# --- paths -------------------------------------------------------------------
data_dir    <- Sys.getenv("FIG_DATA", unset = NA_character_)
results_dir <- Sys.getenv("FIG_OUT",  unset = NA_character_)
if (is.na(data_dir) || is.na(results_dir))
  stop("FIG_DATA and FIG_OUT must be set. Run this figure with ./run.sh, or:\n",
       "  FIG_DATA=data/2e FIG_OUT=output/2e Rscript scripts/2e/moi_curves.R")
if (!dir.exists(data_dir)) stop("FIG_DATA does not exist: ", data_dir)
for (f in c(Sys.getenv("MOI_DATA", "Ludmila_MOI_GC_Data.txt"),
            Sys.getenv("MOI_LAYOUT", "plate_layout.tsv")))
  if (!file.exists(file.path(data_dir, f)))
    stop("Input not found: ", file.path(data_dir, f),
         "\n  Present in that directory: ", paste(list.files(data_dir), collapse = ", "),
         "\n  If the file was renamed, set MOI_DATA or MOI_LAYOUT.")
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

# --- analysis parameters -----------------------------------------------------
DATA_FILE   <- Sys.getenv("MOI_DATA",   unset = "Ludmila_MOI_GC_Data.txt")
LAYOUT_FILE <- Sys.getenv("MOI_LAYOUT", unset = "plate_layout.tsv")
STRAIN    <- Sys.getenv("MOI_STRAIN",   unset = "USA300")
MAX_HOURS <- as.numeric(Sys.getenv("MOI_MAX_HOURS", unset = "16"))
SHOW_MOI  <- c("1", "0.1", "0.01", "0.001", "Bacteria Only", "Blank")

# --- figure aesthetics -------------------------------------------------------
fig_width <- 8.5; fig_height <- 6
plot_margin <- margin(t = 12, r = 14, b = 10, l = 36, unit = "pt")
# Panels are assembled several to a page in Canva and shrunk again, so default
# type reads far too small in the finished layout. Every size is a multiple of
# FIG_TEXT_SCALE; figure dimensions are unaffected.
FIG_TEXT_SCALE <- as.numeric(env_or("FIG_TEXT_SCALE", 2))

# theme_classic rather than theme_bw: axes only, no gridlines, no panel box.
theme_ludmila <- function() {
  s <- FIG_TEXT_SCALE
  theme_classic(base_size = 13 * s, base_family = "Times New Roman") +
    theme(
      axis.line        = element_line(colour = "black", linewidth = 0.6),
      axis.ticks       = element_line(colour = "black"),
      axis.title       = element_text(face = "bold"),
      axis.title.x     = element_text(face = "bold", size = 14 * s, margin = margin(t = 12, unit = "pt")),
      axis.title.y     = element_text(face = "bold", size = 14 * s, margin = margin(r = 12, unit = "pt")),
      axis.text        = element_text(face = "bold"),
      axis.text.x      = element_text(face = "bold", margin = margin(t = 6, unit = "pt")),
      axis.text.y      = element_text(face = "bold", margin = margin(r = 6, unit = "pt")),
      aspect.ratio     = fig_height / fig_width
    )
}

# --- inputs ------------------------------------------------------------------
raw    <- read.delim(file.path(data_dir, DATA_FILE),   check.names = FALSE)
layout <- read.delim(file.path(data_dir, LAYOUT_FILE), check.names = FALSE)

missing <- setdiff(layout$well, names(raw))
if (length(missing))
  stop("plate_layout.tsv names wells absent from the data file: ", paste(missing, collapse = ", "))

# --- long format, one row per well per timepoint ------------------------------
long <- raw %>%
  mutate(Hours = as.numeric(lubridate::hms(Time)) / 3600) %>%
  select(Hours, all_of(layout$well)) %>%
  pivot_longer(-Hours, names_to = "well", values_to = "od") %>%
  left_join(layout, by = "well")

# --- label each series --------------------------------------------------------
long <- long %>%
  mutate(series = case_when(
    role == "blank"           ~ "Blank",
    moi  == "bacteria_only"   ~ "Bacteria Only",
    TRUE                      ~ moi
  ))

# --- mean and SD across replicate wells ---------------------------------------
summarised <- long %>%
  group_by(Hours, strain, series) %>%
  summarise(avg = mean(od, na.rm = TRUE), sd = sd(od, na.rm = TRUE), n = sum(!is.na(od)),
            .groups = "drop")

# --- blank subtraction, per timepoint ------------------------------------------
blanks <- long %>%
  filter(role == "blank") %>%
  group_by(Hours) %>%
  summarise(blank_avg = mean(od, na.rm = TRUE), .groups = "drop")

plot_data <- summarised %>%
  left_join(blanks, by = "Hours") %>%
  mutate(avg = avg - blank_avg) %>%
  select(-blank_avg) %>%
  filter(Hours <= MAX_HOURS, series %in% SHOW_MOI,
         strain == STRAIN | series == "Blank") %>%
  mutate(series = factor(series, levels = SHOW_MOI))

palette <- c("1" = "#D55E00", "0.1" = "#0072B2", "0.01" = "#009E73",
             "0.001" = "#CC79A7", "Bacteria Only" = "#A6761D", "Blank" = "#666666")

# --- plot ---------------------------------------------------------------------
p <- ggplot(plot_data, aes(x = Hours, y = avg, color = series)) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  geom_line(linewidth = 0.8) +
  geom_ribbon(aes(ymin = avg - sd, ymax = avg + sd, fill = series), alpha = 0.2, color = NA) +
  labs(x = "Time (hours)", y = "Absorbance (OD600)", tag = "B") +
  theme_ludmila() +
  theme(legend.position = "none", plot.margin = plot_margin,
        plot.tag = element_text(size = 18 * FIG_TEXT_SCALE, face = "bold",
                                margin = margin(r = 12, unit = "pt")),
        plot.tag.position = "topleft")

out <- file.path(results_dir, paste0("Ludmila_MOI_", STRAIN))
ggsave(paste0(out, ".png"), plot = p, width = fig_width, height = fig_height,
       units = "in", dpi = FIG_DPI, bg = "white")
pdf_ok <- save_pdf_best_effort(p, paste0(out, ".pdf"), width = fig_width,
                               height = fig_height, units = "in", bg = "white")

readr::write_tsv(plot_data, file.path(results_dir, paste0("Ludmila_MOI_", STRAIN, "_source_data.tsv")))
if (pdf_ok) {
  message("wrote ", out, ".png / .pdf and the source-data table")
} else {
  message("wrote ", out, ".png and the source-data table",
          "   (no .pdf: cairo_pdf wrote no file)")
}

# --- legend as its own file, for panel assembly -------------------------------
#
# The panel draws no legend of its own (legend.position = "none") because the
# figure is assembled in Canva, where one legend serves the whole figure. The
# series are taken from plot_data rather than from `palette`, so the legend can
# never list a curve the panel does not draw: SHOW_MOI, MOI_STRAIN and the
# blank filter all apply before this point.
if (requireNamespace("cowplot", quietly = TRUE)) {
  lv <- levels(droplevels(plot_data$series))
  leg_label <- function(x) ifelse(x %in% c("Bacteria Only", "Blank"), x, paste("MOI", x))
  d_leg <- data.frame(x = seq_along(lv), y = 1,
                      series = factor(lv, levels = lv))
  p_leg <- ggplot(d_leg, aes(x = x, y = y, colour = series)) +
    geom_line(aes(group = series), linewidth = 0.8) +
    geom_point(size = 3) +
    scale_color_manual(values = palette, labels = leg_label, name = NULL) +
    theme_ludmila() +
    theme(
      legend.position = "bottom",
      legend.text     = element_text(face = "bold"),
      plot.background = element_rect(fill = "white", colour = NA)
    )
  leg <- cowplot::get_legend(p_leg)
  if (!is.null(leg)) {
    ggsave(paste0(out, "_legend.png"), cowplot::ggdraw(leg),
           width = 9, height = 1.25, units = "in", dpi = FIG_DPI, bg = "white")
    message("Wrote legend: ", basename(out), "_legend.png  (",
            paste(leg_label(lv), collapse = ", "), ")")
  }
} else {
  message("Package cowplot not installed; skipped the legend file.")
  message("Install with: install.packages(\"cowplot\")")
}
