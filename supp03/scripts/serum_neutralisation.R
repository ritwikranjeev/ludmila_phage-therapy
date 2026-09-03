# Supplementary Figure 3 — efficiency of plating of OCPT-Ludmila against the
# patient's serum, week by week, relative to pre-treatment serum (week 0).
#
# Run with ../run.sh, which supplies FIG_DATA and FIG_OUT.
#
# Two tests are reported side by side, because they disagree and the difference
# is the whole story of this panel:
#
#   Welch t-test on the three plate counts. This is what the published analysis
#   did. It treats the triplicate counts as independent estimates of the titre,
#   so its denominator is the spread between three plates of the same dilution
#   — around 2 to 6 plaques. Almost any week separates from week 0 on that
#   scale.
#
#   Conditional binomial (Poisson) test on the summed counts. Plaques are
#   counts, so the sampling noise on a count of ~130 is about sqrt(130) = 11
#   plaques, several times the observed spread. This test asks whether the week
#   differs from week 0 by more than counting noise alone.
#
# The script also reports the log10 reduction, which is the scale neutralisation
# is conventionally judged on, and a dispersion check that shows why the two
# tests disagree.
#
# Optional overrides:
#   S3_FILE      input file name inside FIG_DATA
#   S3_REF_WEEK  week used as the reference (default W0)

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

data_dir    <- need("FIG_DATA")
results_dir <- need("FIG_OUT")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

INFILE   <- env_or("S3_FILE", "serum_eop.tsv")
REF_WEEK <- env_or("S3_REF_WEEK", "W0")

fig_width  <- 9
fig_height <- 5.5

# ---- read ------------------------------------------------------------------

df <- utils::read.delim(file.path(data_dir, INFILE), check.names = FALSE,
                        stringsAsFactors = FALSE)
if (!"Week" %in% names(df)) stop("Column 'Week' missing from ", INFILE)

count_cols <- grep("^Plaques_[0-9]", names(df), value = TRUE)
if (length(count_cols) < 2L) {
  stop("Need at least two Plaques_* replicate columns in ", INFILE,
       " (found: ", paste(names(df), collapse = ", "), ")")
}
counts <- as.matrix(df[, count_cols, drop = FALSE])
storage.mode(counts) <- "double"
rownames(counts) <- df$Week

if (!REF_WEEK %in% df$Week) {
  stop("Reference week '", REF_WEEK, "' not present. Weeks: ",
       paste(df$Week, collapse = ", "))
}
ref_row <- match(REF_WEEK, df$Week)
ref_counts <- counts[ref_row, ]
ref_mean   <- mean(ref_counts, na.rm = TRUE)

# Titres in this file are the plate counts multiplied by a constant dilution
# factor, so they carry no information the counts do not. Working from counts
# keeps the Poisson test meaningful.
titer_cols <- grep("^Titer_[0-9]", names(df), value = TRUE)
if (length(titer_cols) == length(count_cols)) {
  ratio <- as.matrix(df[, titer_cols, drop = FALSE]) / counts
  if (diff(range(ratio, na.rm = TRUE)) > 1e-6) {
    message("NOTE: titre / count ratio is not constant across wells (",
            signif(min(ratio, na.rm = TRUE), 4), " to ",
            signif(max(ratio, na.rm = TRUE), 4),
            "); EOP is computed from counts regardless")
  }
}

# ---- EOP -------------------------------------------------------------------

eop_repl <- counts / ref_mean
df$EOP_mean <- rowMeans(eop_repl, na.rm = TRUE)
df$EOP_sd   <- apply(eop_repl, 1L, stats::sd, na.rm = TRUE)
df$Week     <- factor(df$Week, levels = df$Week)

# ---- statistics ------------------------------------------------------------

is_ref <- as.character(df$Week) == REF_WEEK
n_week <- nrow(df)

welch_p <- rep(NA_real_, n_week)
pois_p  <- rep(NA_real_, n_week)
for (i in seq_len(n_week)) {
  if (is_ref[i]) next
  x <- counts[i, ][!is.na(counts[i, ])]
  welch_p[i] <- tryCatch(
    stats::t.test(ref_counts, x, var.equal = FALSE)$p.value,
    error = function(e) NA_real_
  )
  # Conditional on the total, the split between week and reference is binomial
  # with p = 0.5 under equal plating efficiency.
  pois_p[i] <- stats::binom.test(round(sum(x)),
                                 round(sum(x) + sum(ref_counts, na.rm = TRUE)),
                                 p = 0.5)$p.value
}

stats_df <- data.frame(
  week            = as.character(df$Week),
  n               = apply(counts, 1L, function(v) sum(!is.na(v))),
  mean_count      = rowMeans(counts, na.rm = TRUE),
  sd_count        = apply(counts, 1L, stats::sd, na.rm = TRUE),
  eop             = df$EOP_mean,
  log10_reduction = -log10(df$EOP_mean),
  welch_p         = welch_p,
  welch_q_bh      = stats::p.adjust(welch_p, method = "BH"),
  poisson_p       = pois_p,
  poisson_q_bh    = stats::p.adjust(pois_p, method = "BH"),
  stringsAsFactors = FALSE, row.names = NULL
)
# Dispersion: observed spread between plates against the sqrt(count) a Poisson
# process would give. Values well below 1 mean the replicates are plates of one
# dilution rather than independent platings, which is what makes the Welch
# denominator small.
stats_df$dispersion_ratio <- stats_df$sd_count / sqrt(stats_df$mean_count)

utils::write.table(stats_df, file.path(results_dir, "serum_neutralisation_stats.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

max_red <- max(stats_df$log10_reduction, na.rm = TRUE)
message("Largest reduction in any week: ", signif(max_red, 3),
        " log10 (EOP ", signif(min(stats_df$eop, na.rm = TRUE), 3), ")")
message("Weeks differing from ", REF_WEEK, " at q < 0.05 — Welch: ",
        paste(stats_df$week[!is.na(stats_df$welch_q_bh) & stats_df$welch_q_bh < 0.05],
              collapse = ", "))
message("                                      Poisson: ",
        paste(stats_df$week[!is.na(stats_df$poisson_q_bh) & stats_df$poisson_q_bh < 0.05],
              collapse = ", "))
message("Median dispersion ratio (1 = Poisson): ",
        signif(stats::median(stats_df$dispersion_ratio, na.rm = TRUE), 3))
if (max_red < 1) {
  message("No week reaches a 1 log10 reduction, the conventional threshold for ",
          "neutralising activity.")
}

# ---- plot ------------------------------------------------------------------

lv <- levels(df$Week)
fill_vals <- stats::setNames(rep("#3a3a3a", length(lv)), lv)
fill_vals[REF_WEEK] <- "white"

points_df <- data.frame(
  Week = factor(rep(as.character(df$Week), times = ncol(eop_repl)), levels = lv),
  EOP  = as.vector(eop_repl),
  stringsAsFactors = FALSE
)
points_df <- points_df[!is.na(points_df$EOP), , drop = FALSE]

p <- ggplot(df, aes(x = Week, y = EOP_mean, fill = Week)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "#718096", linewidth = 0.35) +
  geom_col(width = 0.65, colour = "gray20", linewidth = 0.6) +
  geom_errorbar(aes(ymin = EOP_mean - EOP_sd, ymax = EOP_mean + EOP_sd),
                width = 0.22, linewidth = 0.45, colour = "gray20") +
  geom_point(data = points_df, aes(x = Week, y = EOP), inherit.aes = FALSE,
             position = position_jitter(width = 0.08, height = 0, seed = 1L),
             size = 1.8, alpha = 0.75, colour = "gray30") +
  scale_fill_manual(values = fill_vals, guide = "none") +
  scale_x_discrete(labels = function(x) sub("^W", "", x), expand = expansion(add = 0.5)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Week", y = paste0("EOP relative to week ", sub("^W", "", REF_WEEK))) +
  theme_bw(base_size = 13, base_family = FONT_FAMILY) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.margin        = margin(t = 12, r = 14, b = 10, l = 36, unit = "pt"),
    axis.title.x       = element_text(face = "bold", size = 14, family = FONT_FAMILY,
                                      margin = margin(t = 12, unit = "pt")),
    axis.title.y       = element_text(face = "bold", size = 14, family = FONT_FAMILY,
                                      margin = margin(r = 12, unit = "pt")),
    axis.text          = element_text(face = "bold", family = FONT_FAMILY),
    axis.text.x        = element_text(hjust = 0.5, family = FONT_FAMILY,
                                      margin = margin(t = 6, unit = "pt")),
    axis.text.y        = element_text(face = "bold", family = FONT_FAMILY,
                                      margin = margin(r = 6, unit = "pt"))
  )

out_base <- file.path(results_dir, "serum_neutralisation")
ggsave(paste0(out_base, ".png"), p, width = fig_width, height = fig_height,
       dpi = FIG_DPI, bg = "white")
if (save_pdf_best_effort(p, paste0(out_base, ".pdf"), width = fig_width,
                         height = fig_height, bg = "white")) {
  message("Wrote: ", out_base, ".png, .pdf and serum_neutralisation_stats.tsv")
} else {
  message("Wrote: ", out_base, ".png and serum_neutralisation_stats.tsv",
          "   (no .pdf: cairo_pdf wrote no file)")
}
