# Figure 6 — shared helpers for panels 6a and 6b.
# Sourced by the panel scripts; sets no paths and writes no output itself.
#
# Both panels are the same plot: mean zone of inhibition +/- SE, two groups per
# antibiotic, replicate points overlaid, Welch t-test per antibiotic with a
# bracket and stars. 6a facets by strain, 6b does not.
#
# Base R + ggplot2 only. cowplot is optional and used only for the legend file.

if (!interactive()) pdf(NULL)  # prevent Rplots.pdf when run via Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

# ---- small env helpers ----------------------------------------------------

env_or <- function(v, default) {
  x <- Sys.getenv(v, unset = "")
  if (!nzchar(x)) default else x
}

env_csv <- function(v, default = NULL) {
  x <- Sys.getenv(v, unset = "")
  if (!nzchar(x)) return(default)
  trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
}

# Raster resolution for every PNG in this figure. Panels are assembled in Canva
# and exported to PDF from there, so PNG is the delivered output. Width and
# height are NOT configurable: they match the published originals.
FIG_DPI <- as.numeric(env_or("FIG_DPI", 600))

# FIG_SUBTITLE=0 suppresses the methods subtitle without removing it — the text
# is still built by the caller, it just is not drawn. Use it when the caption
# will carry that information instead.
FIG_SUBTITLE <- as.integer(env_or("FIG_SUBTITLE", 1))

# ---- theme ----------------------------------------------------------------

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
      strip.background = element_rect(fill = "white", colour = "black"),
      strip.text       = element_text(face = "bold"),
      aspect.ratio     = fig_height / fig_width
    )
}

# ---- reading --------------------------------------------------------------

# Wide -> long. Every column whose name ends in _Diam is an antibiotic.
# Nothing about which strains, groups or antibiotics exist is hard-coded: the
# levels come from the file unless an explicit order is supplied.
read_kirby_long <- function(path,
                            group_col    = "Clone",
                            group_levels = NULL,
                            abx_order    = NULL,
                            strain_order = NULL) {
  if (!file.exists(path)) stop("Input file not found: ", path)
  wide <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)

  for (nm in c("Strain", group_col, "Rep")) {
    if (!nm %in% names(wide)) {
      stop("Column '", nm, "' missing from ", basename(path),
           " (found: ", paste(names(wide), collapse = ", "), ")")
    }
  }

  diam_cols <- grep("_Diam$", names(wide), value = TRUE)
  if (!length(diam_cols)) {
    stop("No *_Diam columns in ", basename(path))
  }

  long <- do.call(rbind, lapply(diam_cols, function(cn) {
    data.frame(
      Strain      = as.character(wide$Strain),
      Group       = as.character(wide[[group_col]]),
      Rep         = as.character(wide$Rep),
      Antibiotic  = sub("_Diam$", "", cn),
      diameter_mm = suppressWarnings(as.numeric(as.character(wide[[cn]]))),
      stringsAsFactors = FALSE
    )
  }))

  abx <- if (is.null(abx_order)) sub("_Diam$", "", diam_cols) else abx_order
  missing_abx <- setdiff(unique(long$Antibiotic), abx)
  if (length(missing_abx)) {
    stop("Antibiotic(s) present in the data but absent from the requested order: ",
         paste(missing_abx, collapse = ", "))
  }
  long$Antibiotic <- factor(long$Antibiotic, levels = abx)

  strains <- if (is.null(strain_order)) unique(long$Strain) else strain_order
  missing_strains <- setdiff(unique(long$Strain), strains)
  if (length(missing_strains)) {
    stop("Strain(s) present in the data but absent from the requested order: ",
         paste(missing_strains, collapse = ", "))
  }
  long$Strain <- factor(long$Strain, levels = strains)

  grps <- if (is.null(group_levels)) unique(long$Group) else group_levels
  missing_grps <- setdiff(unique(long$Group), grps)
  if (length(missing_grps)) {
    stop("Group(s) present in the data but absent from the requested order: ",
         paste(missing_grps, collapse = ", "))
  }
  long$Group <- factor(long$Group, levels = grps)

  long
}

# Zone diameters are recorded CLSI-style: the full diameter including the 6 mm
# disc, so a plate with no visible clearing reads 6 mm, not 0. That removes the
# zero-versus-missing ambiguity the earlier version of this data carried, and
# with it the zero policy the analysis used to need.
#
# This function now only checks the convention holds and records the disc floor,
# so that a sheet reverting to 0 mm entries is caught rather than silently
# analysed on a different scale.
check_disc_floor <- function(long, disc_mm = 6) {
  vals <- long$diameter_mm[!is.na(long$diameter_mm)]
  below <- vals[vals < disc_mm]
  if (length(below)) {
    stop(length(below), " zone diameter(s) below the ", disc_mm, " mm disc: ",
         paste(sort(unique(below)), collapse = ", "),
         ". Diameters must include the disc; a plate with no zone reads ", disc_mm, " mm.")
  }
  at_floor <- sum(vals == disc_mm)
  message("Zone diameters include the ", disc_mm, " mm disc; ", at_floor,
          " reading(s) at the floor (no visible zone), ",
          sum(is.na(long$diameter_mm)), " not measured")
  attr(long, "disc_mm") <- disc_mm
  long
}

# ---- statistics -----------------------------------------------------------

# Two-sample comparison of the group levels within each antibiotic (and within
# each strain when by_strain = TRUE).
#
#   test = "welch"     Welch t-test, what the published analysis used
#   test = "wilcoxon"  Wilcoxon rank-sum, which the original kept in a separate
#                      near-duplicate script (kirby_wilcoxon.R). With n = 5 per
#                      group the smallest attainable p is 0.0079, so it cannot
#                      reach significance after correction across 30 tests —
#                      worth knowing before reading it as a null result.
#
# Returns one row per comparison attempted. A comparison is skipped, with a
# reason, when either group has fewer than 2 usable values or both groups are
# constant (the test is undefined there — e.g. an antibiotic with a zero zone in
# every replicate of both groups).
kirby_welch <- function(long, by_strain = TRUE, test = c("welch", "wilcoxon")) {
  test <- match.arg(test)
  grp <- levels(long$Group)
  if (length(grp) < 2L) stop("kirby_welch: need two group levels, found ", length(grp))
  if (length(grp) > 2L) {
    stop("kirby_welch: more than two group levels (", paste(grp, collapse = ", "), ")")
  }

  keys <- if (by_strain) {
    expand.grid(Strain = levels(long$Strain), Antibiotic = levels(long$Antibiotic),
                stringsAsFactors = FALSE)
  } else {
    data.frame(Strain = NA_character_, Antibiotic = levels(long$Antibiotic),
               stringsAsFactors = FALSE)
  }

  rows <- vector("list", nrow(keys))
  for (i in seq_len(nrow(keys))) {
    sel <- long$Antibiotic == keys$Antibiotic[i] &
      (if (by_strain) as.character(long$Strain) == keys$Strain[i] else TRUE)
    x <- long$diameter_mm[sel & long$Group == grp[1L]]
    y <- long$diameter_mm[sel & long$Group == grp[2L]]
    x <- x[!is.na(x)]
    y <- y[!is.na(y)]

    skip <- NA_character_
    if (length(x) < 2L || length(y) < 2L) {
      skip <- "fewer than 2 usable replicates in one group"
    } else if (stats::sd(x) == 0 && stats::sd(y) == 0) {
      skip <- paste0("both groups constant (no zone in either); ", test, " test undefined")
    }

    tt <- NULL
    if (is.na(skip)) {
      tt <- tryCatch(
        if (test == "welch") {
          stats::t.test(x, y, var.equal = FALSE, paired = FALSE)
        } else {
          suppressWarnings(stats::wilcox.test(x, y, exact = FALSE, correct = TRUE))
        },
        error = function(e) { skip <<- conditionMessage(e); NULL })
    }

    rows[[i]] <- data.frame(
      Strain      = keys$Strain[i],
      Antibiotic  = keys$Antibiotic[i],
      group1      = grp[1L],
      group2      = grp[2L],
      n1          = length(x),
      n2          = length(y),
      mean1       = if (length(x)) mean(x) else NA_real_,
      mean2       = if (length(y)) mean(y) else NA_real_,
      sd1         = if (length(x) > 1L) stats::sd(x) else NA_real_,
      sd2         = if (length(y) > 1L) stats::sd(y) else NA_real_,
      test        = test,
      statistic   = if (is.null(tt)) NA_real_ else unname(tt$statistic),
      df_welch    = if (is.null(tt) || is.null(tt$parameter)) NA_real_ else unname(tt$parameter),
      p           = if (is.null(tt)) NA_real_ else tt$p.value,
      skipped     = skip,
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, rows)

  # Two multiplicity corrections, both reported:
  #   p_adj_within  Benjamini-Hochberg within each strain (what kirby.R did)
  #   p_adj_all     Benjamini-Hochberg across every comparison in the panel
  out$p_adj_all <- stats::p.adjust(out$p, method = "BH")
  out$p_adj_within <- NA_real_
  if (by_strain) {
    for (s in unique(out$Strain)) {
      idx <- which(out$Strain == s)
      out$p_adj_within[idx] <- stats::p.adjust(out$p[idx], method = "BH")
    }
  } else {
    out$p_adj_within <- out$p_adj_all
  }

  out$disc_mm <- attr(long, "disc_mm")
  # Report in display order (strain facet order, then antibiotic order), not
  # alphabetically.
  ord <- order(match(out$Strain, levels(long$Strain)),
               match(out$Antibiotic, levels(long$Antibiotic)))
  out[ord, , drop = FALSE]
}

p_to_stars <- function(p, cutoff = 0.05) {
  ifelse(is.na(p) | p > cutoff, "",
         ifelse(p <= 0.001, "***", ifelse(p <= 0.01, "**", "*")))
}

# ---- summaries and annotation geometry ------------------------------------

summarize_kirby <- function(long, by_strain = TRUE) {
  keys <- if (by_strain) list(long$Strain, long$Antibiotic, long$Group)
          else           list(long$Antibiotic, long$Group)
  sp <- split(long$diameter_mm, keys, drop = TRUE, sep = "\r")
  parts <- strsplit(names(sp), "\r", fixed = TRUE)
  n    <- vapply(sp, function(v) sum(!is.na(v)), integer(1))
  mean_ <- vapply(sp, function(v) if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE), numeric(1))
  sd_   <- vapply(sp, function(v) if (sum(!is.na(v)) > 1L) stats::sd(v, na.rm = TRUE) else NA_real_, numeric(1))

  out <- data.frame(
    Strain     = if (by_strain) vapply(parts, `[`, "", 1L) else NA_character_,
    Antibiotic = vapply(parts, `[`, "", if (by_strain) 2L else 1L),
    Group      = vapply(parts, `[`, "", if (by_strain) 3L else 2L),
    n = as.integer(n), mean = mean_, sd = sd_,
    stringsAsFactors = FALSE, row.names = NULL
  )
  out$se <- out$sd / sqrt(out$n)
  out <- out[out$n > 0L, , drop = FALSE]
  out$Antibiotic <- factor(out$Antibiotic, levels = levels(long$Antibiotic))
  out$Group      <- factor(out$Group,      levels = levels(long$Group))
  if (by_strain) out$Strain <- factor(out$Strain, levels = levels(long$Strain))
  out
}

# Bracket + star for each comparison whose chosen p-value clears the cutoff.
# p_col selects which p-value drives the stars: "p", "p_adj_within", "p_adj_all".
kirby_brackets <- function(stats_df, long, sum_df, p_col = "p", cutoff = 0.05,
                           by_strain = TRUE, half_width = 0.375, gap = 2, tick = 0.8) {
  if (!p_col %in% names(stats_df)) {
    stop("p_col '", p_col, "' not in the statistics table")
  }
  keep <- !is.na(stats_df[[p_col]]) & stats_df[[p_col]] <= cutoff
  if (!any(keep)) return(NULL)
  s <- stats_df[keep, , drop = FALSE]

  # Bracket sits above the taller of (highest replicate point, tallest bar + SE).
  y_of <- function(strain, abx) {
    sel_pt <- long$Antibiotic == abx &
      (if (by_strain) as.character(long$Strain) == strain else TRUE)
    pt <- suppressWarnings(max(long$diameter_mm[sel_pt], na.rm = TRUE))
    sel_bar <- sum_df$Antibiotic == abx &
      (if (by_strain) as.character(sum_df$Strain) == strain else TRUE)
    bar <- suppressWarnings(max(sum_df$mean[sel_bar] +
                                  ifelse(is.na(sum_df$se[sel_bar]), 0, sum_df$se[sel_bar]),
                                na.rm = TRUE))
    v <- c(pt, bar)
    v <- v[is.finite(v)]
    if (!length(v)) 0 else max(v)
  }

  y <- vapply(seq_len(nrow(s)), function(i) y_of(s$Strain[i], s$Antibiotic[i]), numeric(1))
  xc <- match(s$Antibiotic, levels(long$Antibiotic))

  ann <- data.frame(
    Antibiotic = factor(s$Antibiotic, levels = levels(long$Antibiotic)),
    x_min = xc - half_width,
    x_max = xc + half_width,
    x_mid = xc,
    y     = y + gap,
    y_tick = y + gap - tick,
    label = p_to_stars(s[[p_col]], cutoff = cutoff),
    stringsAsFactors = FALSE
  )
  if (by_strain) ann$Strain <- factor(s$Strain, levels = levels(long$Strain))
  ann[nzchar(ann$label), , drop = FALSE]
}

# ---- plot -----------------------------------------------------------------

kirby_barplot <- function(long, sum_df, ann, fill_values,
                          facet_by_strain = TRUE, facet_ncol = 2,
                          fig_width, fig_height,
                          x_lab = "Antibiotic",
                          y_lab = "Zone of Inhibition (mm)",
                          subtitle = NULL) {
  pd <- position_dodge(width = 0.75)
  pts <- long[!is.na(long$diameter_mm), , drop = FALSE]

  p <- ggplot() +
    geom_col(data = sum_df,
             aes(x = Antibiotic, y = mean, fill = Group),
             width = 0.65, position = pd, colour = "gray20", linewidth = 0.6) +
    geom_errorbar(data = sum_df[!is.na(sum_df$se), , drop = FALSE],
                  aes(x = Antibiotic, ymin = mean - se, ymax = mean + se, group = Group),
                  width = 0.2, position = pd, colour = "gray40", linewidth = 0.5) +
    geom_point(data = pts,
               aes(x = Antibiotic, y = diameter_mm, group = Group),
               position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75),
               size = 1.8, alpha = 0.7, colour = "gray30")

  if (!is.null(ann) && nrow(ann)) {
    p <- p +
      geom_segment(data = ann, aes(x = x_min, xend = x_max, y = y, yend = y),
                   inherit.aes = FALSE, linewidth = 0.35, colour = "gray20") +
      geom_segment(data = ann, aes(x = x_min, xend = x_min, y = y, yend = y_tick),
                   inherit.aes = FALSE, linewidth = 0.35, colour = "gray20") +
      geom_segment(data = ann, aes(x = x_max, xend = x_max, y = y, yend = y_tick),
                   inherit.aes = FALSE, linewidth = 0.35, colour = "gray20") +
      geom_text(data = ann, aes(x = x_mid, y = y, label = label),
                inherit.aes = FALSE, size = 4, fontface = "bold", vjust = -0.3)
  }

  if (facet_by_strain) {
    p <- p + facet_wrap(~ Strain, ncol = facet_ncol, scales = "free_y")
  }

  p +
    scale_x_discrete(expand = expansion(add = 0.5), drop = FALSE) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    scale_fill_manual(values = fill_values, drop = FALSE) +
    labs(x = x_lab, y = y_lab, fill = NULL,
         subtitle = if (FIG_SUBTITLE) subtitle else NULL) +
    theme_ludmila(fig_width, fig_height) +
    theme(
      aspect.ratio    = NULL,
      axis.text.x     = element_text(angle = 45, hjust = 1),
      plot.subtitle   = element_text(size = 9, hjust = 0, margin = margin(b = 8, unit = "pt")),
      legend.position = "none",
      panel.spacing   = unit(24, "pt")
    ) +
    coord_cartesian(clip = "off")
}

# Legend as its own file, so the assembled figure carries one legend for both
# panels rather than one per panel.
save_group_legend <- function(fill_values, outfile, title = NULL) {
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    message("Package cowplot not installed; skipped legend file: ", outfile)
    message("Install with: install.packages(\"cowplot\")")
    return(invisible(NULL))
  }
  lv <- names(fill_values)
  d <- data.frame(Group = factor(lv, levels = lv), y = rep(1, length(lv)))
  p_leg <- ggplot(d, aes(x = Group, y = y, fill = Group)) +
    geom_col(width = 0.6, colour = "gray20", linewidth = 0.6) +
    scale_fill_manual(values = fill_values, name = title) +
    theme_bw(base_size = 13, base_family = "Times New Roman") +
    theme(
      legend.position = "bottom",
      legend.title    = element_text(face = "bold"),
      legend.text     = element_text(face = "bold"),
      plot.background = element_rect(fill = "white", colour = NA)
    )
  leg <- cowplot::get_legend(p_leg)
  if (is.null(leg)) return(invisible(NULL))
  ggsave(outfile, cowplot::ggdraw(leg), width = 3.5, height = 1.25,
         units = "in", dpi = FIG_DPI, bg = "white")
  message("Wrote legend: ", outfile)
  invisible(outfile)
}

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

# Every panel writes a PNG; the PDF is best-effort and never claimed unless written.
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
