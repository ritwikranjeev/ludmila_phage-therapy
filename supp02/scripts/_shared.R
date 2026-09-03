# Supplementary Figure 2 — shared helpers for the stability panels (2b, 2c).
# Sourced by the panel scripts; sets no paths and writes no output itself.
#
# Panels 2a (proteomic tree) and 2d (VIRIDIC heatmap) are unrelated plots with
# their own dependencies and do not use this file.
#
# Base R + ggplot2 only.

if (!interactive()) pdf(NULL)

suppressPackageStartupMessages({
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

# FIG_SUBTITLE=0 suppresses the methods subtitle without removing it — the text
# is still built by the caller, it just is not drawn.
FIG_SUBTITLE <- as.integer(env_or("FIG_SUBTITLE", 1))

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

theme_ludmila <- function(fig_width, fig_height) {
  theme_bw(base_size = 13, base_family = "Times New Roman") +
    theme(
      axis.title       = element_text(face = "bold", colour = "black"),
      axis.title.x     = element_text(face = "bold", size = 14, colour = "black",
                                      margin = margin(t = 12, unit = "pt")),
      axis.title.y     = element_text(face = "bold", size = 14, colour = "black",
                                      margin = margin(r = 12, unit = "pt")),
      axis.text        = element_text(face = "bold", colour = "black"),
      axis.text.x      = element_text(face = "bold", colour = "black",
                                      margin = margin(t = 6, unit = "pt")),
      axis.text.y      = element_text(face = "bold", colour = "black",
                                      margin = margin(r = 6, unit = "pt")),
      panel.grid.minor = element_blank(),
      aspect.ratio     = fig_height / fig_width
    )
}

# ---- input -----------------------------------------------------------------

# Both stability files have the same shape: column 1 is the condition (pH or
# temperature), and every column named Titer_R<n> is one replicate titre.
# Nothing about how many replicates or which conditions exist is hard-coded.
read_titer_wide <- function(path) {
  if (!file.exists(path)) stop("Input file not found: ", path)
  wide <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  titer_cols <- grep("^Titer_R[0-9]", names(wide), value = TRUE)
  if (length(titer_cols) < 2L) {
    stop("Need at least two Titer_R* replicate columns in ", basename(path),
         " (found: ", paste(names(wide), collapse = ", "), ")")
  }
  x_col <- names(wide)[1L]
  x <- suppressWarnings(as.numeric(wide[[x_col]]))
  if (anyNA(x)) stop("Non-numeric value in the condition column '", x_col, "' of ", basename(path))

  mat <- as.matrix(wide[, titer_cols, drop = FALSE])
  storage.mode(mat) <- "double"

  long <- data.frame(
    x         = rep(x, times = length(titer_cols)),
    replicate = rep(titer_cols, each = nrow(wide)),
    titer     = as.vector(mat),
    stringsAsFactors = FALSE
  )

  list(x_col = x_col, x = x, titer_cols = titer_cols, mat = mat,
       long = long[!is.na(long$titer), , drop = FALSE])
}

titer_summary <- function(d) {
  n    <- apply(d$mat, 1L, function(v) sum(!is.na(v)))
  mean_ <- rowMeans(d$mat, na.rm = TRUE)
  sd_   <- apply(d$mat, 1L, function(v) if (sum(!is.na(v)) > 1L) stats::sd(v, na.rm = TRUE) else NA_real_)
  data.frame(x = d$x, n = n, mean = mean_, sd = sd_,
             se = sd_ / sqrt(n), stringsAsFactors = FALSE)
}

# ---- statistics ------------------------------------------------------------

# One-way ANOVA across conditions, followed by Tukey HSD on all pairs. The
# original scripts computed these and printed them to the console, where they
# were lost; here they are returned so the panel can write them next to the
# figure. Descriptive only — no stars are drawn on the panels.
anova_tukey <- function(d) {
  df <- data.frame(titer = d$long$titer,
                   grp   = factor(as.character(d$long$x),
                                  levels = as.character(d$x)))
  if (nlevels(df$grp) < 2L) return(NULL)
  fit <- stats::aov(titer ~ grp, data = df)
  s   <- summary(fit)[[1L]]

  aov_tab <- data.frame(
    term      = trimws(rownames(s)),
    df        = s[["Df"]],
    sum_sq    = s[["Sum Sq"]],
    mean_sq   = s[["Mean Sq"]],
    f_value   = s[["F value"]],
    p_value   = s[["Pr(>F)"]],
    stringsAsFactors = FALSE
  )

  tk <- stats::TukeyHSD(fit)$grp
  tukey_tab <- data.frame(
    comparison = rownames(tk),
    diff       = tk[, "diff"],
    lwr        = tk[, "lwr"],
    upr        = tk[, "upr"],
    p_adj      = tk[, "p adj"],
    stringsAsFactors = FALSE, row.names = NULL
  )
  tukey_tab <- tukey_tab[order(tukey_tab$p_adj), , drop = FALSE]

  list(anova = aov_tab, tukey = tukey_tab)
}

write_stats <- function(res, results_dir, prefix) {
  if (is.null(res)) return(invisible(NULL))
  utils::write.table(res$anova, file.path(results_dir, paste0(prefix, "_anova.tsv")),
                     sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  utils::write.table(res$tukey, file.path(results_dir, paste0(prefix, "_tukey.tsv")),
                     sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  message("ANOVA p = ", signif(res$anova$p_value[1L], 3),
          "; Tukey pairs with p adj < 0.05: ", sum(res$tukey$p_adj < 0.05),
          " of ", nrow(res$tukey))
  invisible(NULL)
}

# ---- shared plot furniture -------------------------------------------------

# Titres are plotted in units of 10^8 PFU/mL, as in the published panels.
titer_y_scale <- function(divisor = 1e8, exponent = 8) {
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05)),
    labels = function(v) v / divisor,
    name   = bquote(bold("Phage Titer (") * bold(10^.(exponent)) * bold(" PFU/mL)"))
  )
}

save_panel <- function(p, out_base, fig_width, fig_height) {
  ggsave(paste0(out_base, ".png"), plot = p, width = fig_width, height = fig_height,
         units = "in", dpi = FIG_DPI, bg = "white")
  if (save_pdf_best_effort(p, paste0(out_base, ".pdf"), width = fig_width,
                           height = fig_height, units = "in", bg = "white")) {
    message("Wrote: ", out_base, ".png and .pdf")
  } else {
    message("Wrote: ", out_base, ".png   (no .pdf: cairo_pdf wrote no file)")
  }
  invisible(out_base)
}
