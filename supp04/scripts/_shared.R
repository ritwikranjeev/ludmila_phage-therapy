# Supplementary Figure 4 — shared helpers for panels 4a and 4b.
# Sourced by the panel scripts; sets no paths and writes no output itself.
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

# FIG_SUBTITLE=0 suppresses the methods subtitle without removing it — the text
# is still built by the caller, it just is not drawn.
FIG_SUBTITLE <- as.integer(env_or("FIG_SUBTITLE", 1))
env_csv <- function(v, default) {
  x <- Sys.getenv(v, unset = "")
  if (!nzchar(x)) return(default)
  trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
}
env_num <- function(v, default) {
  x <- Sys.getenv(v, unset = "")
  if (!nzchar(x)) return(default)
  as.numeric(trimws(strsplit(x, ",", fixed = TRUE)[[1L]]))
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
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      aspect.ratio     = fig_height / fig_width
    )
}

# ---- input -----------------------------------------------------------------

# Time course in wide form: column 1 is time, every other column is one
# replicate of one group, identified by a prefix in its name.
read_timecourse <- function(path, prefixes, labels = prefixes) {
  if (!file.exists(path)) stop("Input file not found: ", path)
  d <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  time <- suppressWarnings(as.numeric(d[[1L]]))
  if (anyNA(time)) stop("Non-numeric value in the time column of ", basename(path))

  groups <- list()
  for (i in seq_along(prefixes)) {
    cols <- grep(prefixes[i], names(d), fixed = TRUE)
    if (!length(cols)) {
      stop("No columns matching '", prefixes[i], "' in ", basename(path),
           " (found: ", paste(names(d), collapse = ", "), ")")
    }
    m <- as.matrix(d[, cols, drop = FALSE])
    storage.mode(m) <- "double"
    groups[[labels[i]]] <- m
  }
  list(time = time, groups = groups, file = basename(path))
}

# ---- statistics ------------------------------------------------------------

welch <- function(x, y, label) {
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  if (length(x) < 2L || length(y) < 2L) {
    return(data.frame(label = label, n1 = length(x), n2 = length(y),
                      mean1 = mean(x), mean2 = mean(y),
                      sd1 = NA_real_, sd2 = NA_real_, t = NA_real_,
                      df = NA_real_, p = NA_real_, ratio = mean(x) / mean(y),
                      stringsAsFactors = FALSE))
  }
  tt <- stats::t.test(x, y, var.equal = FALSE)
  data.frame(
    label = label, n1 = length(x), n2 = length(y),
    mean1 = mean(x), mean2 = mean(y),
    sd1 = stats::sd(x), sd2 = stats::sd(y),
    t = unname(tt$statistic), df = unname(tt$parameter), p = tt$p.value,
    ratio = mean(x) / mean(y),
    stringsAsFactors = FALSE
  )
}

# ---- plot ------------------------------------------------------------------

# Two-bar comparison with replicate points, used by both panels.
two_group_bar <- function(values, labels, y_lab, fig_width, fig_height,
                          fill_values = c("white", "gray20"), subtitle = NULL) {
  sum_df <- data.frame(
    group = factor(labels, levels = labels),
    mean  = vapply(values, mean, numeric(1)),
    sd    = vapply(values, function(v) if (length(v) > 1L) stats::sd(v) else NA_real_, numeric(1)),
    stringsAsFactors = FALSE
  )
  pts <- do.call(rbind, lapply(seq_along(values), function(i) {
    data.frame(group = factor(labels[i], levels = labels), value = values[[i]],
               stringsAsFactors = FALSE)
  }))

  ggplot(sum_df, aes(x = group, y = mean, fill = group)) +
    geom_col(width = 0.6, colour = "gray20", linewidth = 0.6) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                  width = 0.2, colour = "gray20", linewidth = 0.5) +
    geom_point(data = pts, aes(x = group, y = value), inherit.aes = FALSE,
               position = position_jitter(width = 0.08, height = 0, seed = 1L),
               size = 1.9, alpha = 0.8, colour = "gray30") +
    scale_fill_manual(values = stats::setNames(fill_values, labels), guide = "none") +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
    labs(x = NULL, y = y_lab, subtitle = if (FIG_SUBTITLE) subtitle else NULL) +
    theme_ludmila(fig_width, fig_height) +
    theme(plot.subtitle = element_text(size = 9, hjust = 0,
                                       margin = margin(b = 8, unit = "pt")))
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

write_tsv <- function(df, path) {
  utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  message("Wrote: ", basename(path))
}
