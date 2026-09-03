#!/usr/bin/env Rscript
# EOP bar chart from pre-post-eop.tsv.
# Figure typography and dimensions match OSGC/scripts/osgc.R (Times New Roman, 8.5×6 in, margins, black frame).
# Bar fill colors unchanged from original logic (user may restyle).

suppressPackageStartupMessages({

# ---- paths -------------------------------------------------------------------
data_dir    <- Sys.getenv("FIG_DATA", unset = NA_character_)
results_dir <- Sys.getenv("FIG_OUT",  unset = NA_character_)
if (is.na(data_dir) || is.na(results_dir))
  stop("FIG_DATA and FIG_OUT must be set. Run this figure with ./run.sh")
if (!dir.exists(data_dir)) stop("FIG_DATA does not exist: ", data_dir)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
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

# Match OSGC figure constants
fig_width   <- 8.5
fig_height  <- 6
plot_margin <- ggplot2::margin(t = 12, r = 14, b = 10, l = 36, unit = "pt")

theme_ludmila_bar <- function() {
  theme_bw(base_size = 13, base_family = "Times New Roman") +
    theme(
      axis.title       = element_text(face = "bold", colour = "black"),
      axis.title.x     = element_text(face = "bold", size = 14, colour = "black", margin = margin(t = 12, unit = "pt")),
      axis.title.y     = element_text(face = "bold", size = 14, colour = "black", margin = margin(r = 12, unit = "pt")),
      axis.text        = element_text(face = "bold", colour = "black"),
      axis.text.x      = element_text(face = "bold", colour = "black", margin = margin(t = 6, unit = "pt")),
      axis.text.y      = element_text(face = "bold", colour = "black", margin = margin(r = 6, unit = "pt")),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.border       = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      axis.ticks         = element_line(colour = "black", linewidth = 0.5),
      axis.line          = element_line(colour = "black", linewidth = 0.5),
      aspect.ratio       = NULL,
      plot.margin        = plot_margin
    )
}

# Inputs and outputs come from run.sh via FIG_DATA/FIG_OUT, like every other
# panel. This script used to resolve both relative to its own location on disk
# (scripts/4c/pre-post-eop.tsv, scripts/4c/figures/), which does not exist.
data_path <- file.path(data_dir, "pre-post-eop.tsv")
if (!file.exists(data_path))
  stop("Input not found: ", data_path,
       "\n  Present in that directory: ", paste(list.files(data_dir), collapse = ", "))
out_dir <- results_dir

df <- read.delim(data_path, stringsAsFactors = FALSE)
df$Group <- factor(df$Isolate, levels = c("Pre", "Post"))

titer_cols <- grep("^Titer_[0-9]+$", names(df), value = TRUE)
if (length(titer_cols) < 2) {
  stop("Need at least 2 titer replicate columns matching ^Titer_[0-9]+$.")
}
tit <- as.matrix(df[, titer_cols, drop = FALSE])
ref <- mean(tit[df$Isolate == "Pre", ])
eop_repl <- tit / ref
df$EOP_mean <- rowMeans(eop_repl)
df$EOP_sd <- apply(eop_repl, 1, stats::sd)
df$EOP_ymin <- df$EOP_mean - df$EOP_sd
df$EOP_ymax <- df$EOP_mean + df$EOP_sd

COL_BAR_GREY <- "#3a3a3a"

lv <- levels(df$Group)
fill_vals <- setNames(rep(COL_BAR_GREY, length(lv)), lv)
fill_vals["Pre"] <- "white"

points_df <- data.frame(
  Group = factor(rep(as.character(df$Group), each = ncol(tit)), levels = lv),
  EOP = as.vector(t(tit / ref))
)

pre_vals <- as.numeric(eop_repl[df$Isolate == "Pre", ])
post_vals <- as.numeric(eop_repl[df$Isolate == "Post", ])
welch <- t.test(pre_vals, post_vals, var.equal = FALSE)
sig_label <- if (welch$p.value < 0.001) {
  "***"
} else if (welch$p.value < 0.01) {
  "**"
} else if (welch$p.value < 0.05) {
  "*"
} else {
  "ns"
}
stats_out <- file.path(out_dir, "eop_welch_ttest.tsv")
stats_df <- data.frame(
  group_1 = "Pre",
  group_2 = "Post",
  mean_pre = mean(pre_vals),
  mean_post = mean(post_vals),
  n_pre = length(pre_vals),
  n_post = length(post_vals),
  t_statistic = unname(welch$statistic),
  df = unname(welch$parameter),
  p_value = welch$p.value,
  conf_low = welch$conf.int[1],
  conf_high = welch$conf.int[2],
  stringsAsFactors = FALSE
)
write.table(stats_df, stats_out, sep = "\t", quote = FALSE, row.names = FALSE)

max_y <- max(c(df$EOP_ymax, points_df$EOP), na.rm = TRUE)
ann_base <- max_y * 1.05
ann_tip <- max_y * 1.02
ann_text <- max_y * 1.08

p <- ggplot(df, aes(x = Group, y = EOP_mean, fill = Group)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#718096", linewidth = 0.35) +
  geom_col(width = 0.65, color = "gray20", linewidth = 0.6) +
  geom_errorbar(
    aes(ymin = EOP_ymin, ymax = EOP_ymax),
    width = 0.22,
    linewidth = 0.45,
    color = "gray20"
  ) +
  geom_point(
    data = points_df,
    aes(x = Group, y = EOP),
    inherit.aes = FALSE,
    position = position_jitter(width = 0.08, height = 0),
    size = 1.8,
    alpha = 0.75,
    color = "gray30"
  ) +
  annotate("segment", x = 1, xend = 1, y = ann_tip, yend = ann_base, linewidth = 0.5, color = "gray20") +
  annotate("segment", x = 2, xend = 2, y = ann_tip, yend = ann_base, linewidth = 0.5, color = "gray20") +
  annotate("segment", x = 1, xend = 2, y = ann_base, yend = ann_base, linewidth = 0.5, color = "gray20") +
  annotate(
    "text", x = 1.5, y = ann_text, label = sig_label,
    size = 5.2, family = "Times New Roman", fontface = "bold"
  ) +
  scale_fill_manual(values = fill_vals, guide = "none") +
  scale_x_discrete(
    labels = c(Pre = "Pre-Treatment", Post = "Post-Treatment"),
    expand = expansion(add = 0.5)
  ) +
  scale_y_continuous(limits = c(0, ann_text * 1.05), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Isolate", y = "EOP relative to Pre-Treatment") +
  theme_ludmila_bar() +
  coord_cartesian(clip = "off")

p <- p + theme(
  axis.title.y      = element_text(face = "bold", colour = "black"),
  axis.title.y.left = element_text(face = "bold", colour = "black"),
  axis.text.y       = element_text(face = "bold", colour = "black"),
  axis.text.y.left  = element_text(face = "bold", colour = "black"),
  axis.text.x       = element_text(face = "bold", colour = "black")
)

ggsave(
  file.path(out_dir, "eop_bar.png"),
  p,
  width = fig_width,
  height = fig_height,
  units = "in",
  dpi = FIG_DPI,
  bg = "white"
)
pdf_ok <- save_pdf_best_effort(
  p,
  file.path(out_dir, "eop_bar.pdf"),
  width = fig_width,
  height = fig_height,
  units = "in",
  bg = "white"
)

if (pdf_ok) {
  message(
    "Wrote ", file.path(out_dir, "eop_bar.png"),
    ", ", file.path(out_dir, "eop_bar.pdf"),
    ", and ", stats_out
  )
} else {
  message(
    "Wrote ", file.path(out_dir, "eop_bar.png"),
    " and ", stats_out,
    "   (no .pdf: cairo_pdf wrote no file)"
  )
}
