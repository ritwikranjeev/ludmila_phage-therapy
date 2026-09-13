# Supplementary Figure 4b — burst size, pre- vs post-treatment isolate.
#
# Run with ../../run.sh, which supplies FIG_DATA, FIG_OUT and FIG_SHARED.
#
# METHOD.
#
# Plateau method: burst size = (mean titre on the plateau - mean titre before
# lysis) / infective centres. It assumes the one-step growth curve HAS a
# plateau. This script tests that assumption instead of taking it on trust: it
# fits the log slope over the tail of each curve and says so when the titre is
# still climbing at the last timepoint, in which case every burst size here is
# a lower bound.
#
# THE DENOMINATOR IS IN THE DATA, AND WAS NOT USED.
#
#   Burst size is progeny per infected cell, so the denominator is the number of
#   infected cells. In a one-step growth curve that number is measured: the
#   pre-lysis plateau, because each infected cell yields exactly one plaque until
#   it bursts. The original analysis instead divided by an assumed 1e5.
#
#   The measured plateau here is 5.8-7.0 x 10^5, seven times that assumption, and
#   it reconciles with the protocol: at MOI 0.001 with near-complete adsorption
#   the plateau equals the phage input, implying a host density of 6-7 x 10^8
#   CFU/mL in all three experiments. So the data is consistent; the constant was
#   not. Using the measured plateau also self-corrects for any difference in host
#   density or adsorption between the two isolates.
#
#   S4_DENOMINATOR=assumed restores the old behaviour for comparison.
#
#   The plateau window. Burst size roughly doubles depending on which times are
#   called the plateau, so the sweep is written out alongside the headline
#   number.
#
# Optional overrides:
#   S4_OSGC_FILE     pre/post input file inside FIG_DATA
#   S4_USA300_FILE   USA300 input file inside FIG_DATA, for the Figure 2d check
#   S4_GROUPS        comma-separated column prefixes identifying the two groups
#   S4_LABELS        display labels for those groups
#   S4_IC0           assumed infective centres per mL
#   S4_PRE_TIMES     comma-separated pre-lysis timepoints
#   S4_PLATEAU_TIMES comma-separated plateau timepoints

need <- function(v) {
  x <- Sys.getenv(v, unset = NA_character_)
  if (is.na(x)) stop(v, " must be set. Run this figure with ./run.sh")
  x
}
data_dir    <- need("FIG_DATA")
results_dir <- need("FIG_OUT")
shared      <- need("FIG_SHARED")
source(shared)

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

INFILE     <- env_or("S4_OSGC_FILE", "osgc_post-pre-raw.tsv")
USA_FILE   <- env_or("S4_USA300_FILE", "osgc_usa300.tsv")
GROUPS     <- env_csv("S4_GROUPS", c("OCPT04SA01", "OCPT04SA03"))
LABELS     <- env_csv("S4_LABELS", c("Pre", "Post"))
IC0         <- as.numeric(env_or("S4_IC0", "1e5"))
DENOMINATOR <- env_or("S4_DENOMINATOR", "measured")
if (!DENOMINATOR %in% c("measured", "assumed")) {
  stop("S4_DENOMINATOR must be 'measured' (pre-lysis plateau) or 'assumed' (S4_IC0)")
}
PRE_TIMES  <- env_num("S4_PRE_TIMES", c(10, 20, 30))

# PLATEAU AT 40 MIN, AND WHY.
#
# A one-step growth curve measures ONE round of infection. Titre rises roughly
# tenfold between 30 and 40 min in all three curves here, which is the first
# burst; it then keeps climbing to the last timepoint sampled. That continued
# rise is not part of the first burst — progeny released at 40 min re-infect,
# and later timepoints mix the first burst with subsequent cycles. Taking a
# later window therefore inflates burst size with secondary infection.
#
# 40 min is the first post-burst timepoint and is used for all three curves
# (pre-treatment, post-treatment and USA300). Applying it to some curves and
# not others would not be defensible: every curve is still rising at every
# window tried, so there is no curve for which a later window is more plateau-
# like than another.
#
# The single-timepoint window is a real limitation, not a hidden one: one
# timepoint cannot demonstrate a plateau, lysis is asynchronous so some cells
# have not burst by 40 min, and the reported values are consequently lower
# bounds. The full window sweep is written to burst_size_window_sweep.tsv so
# the sensitivity of the estimate is on record rather than buried.
PLAT_TIMES <- env_num("S4_PLATEAU_TIMES", c(40))

# Enlarged from 5.5 x 5 when FIG_TEXT_SCALE doubled the type: point sizes are
# absolute, so on the old canvas the rotated y-axis title ran out of room and
# was clipped. Both dimensions are scaled by the same factor, so the panel
# proportions (aspect.ratio = fig_height / fig_width) are unchanged.
fig_width  <- 8.25
fig_height <- 7.5

d  <- read_timecourse(file.path(data_dir, INFILE), GROUPS, LABELS)
tm <- d$time

if (DENOMINATOR == "measured") {
  message("Denominator: the measured pre-lysis plateau, one value per replicate ",
          "(infected cells, from ", paste(PRE_TIMES, collapse = "/"), " min).")
} else {
  message("Denominator: assumed IC0 = ", format(IC0, scientific = TRUE),
          " — this is the original behaviour and overstates burst size ",
          "wherever the measured plateau is larger.")
}

for (nm in c("PRE_TIMES", "PLAT_TIMES")) {
  v <- get(nm)
  if (!all(v %in% tm)) {
    stop(nm, " contains times not present in the data: ",
         paste(setdiff(v, tm), collapse = ", "))
  }
}

burst_of <- function(m) {
  pre_mean  <- colMeans(m[tm %in% PRE_TIMES, , drop = FALSE], na.rm = TRUE)
  plat_mean <- colMeans(m[tm %in% PLAT_TIMES, , drop = FALSE], na.rm = TRUE)
  denom <- if (DENOMINATOR == "measured") pre_mean else rep(IC0, length(pre_mean))
  (plat_mean - pre_mean) / denom
}

# ---- has the curve plateaued? ----------------------------------------------

tail_slope <- function(m, from = min(PLAT_TIMES)) {
  sel <- tm >= from
  if (sum(sel) < 3L) return(NA_real_)
  y <- log(rowMeans(m[sel, , drop = FALSE], na.rm = TRUE))
  unname(stats::coef(stats::lm(y ~ tm[sel]))[2L])
}

slopes <- vapply(d$groups, tail_slope, numeric(1))
message("Log slope from ", min(PLAT_TIMES), " min to the last timepoint (",
        max(tm), " min):")
for (i in seq_along(slopes)) {
  message(sprintf("  %-5s %+.4f /min  -> %s", LABELS[i], slopes[i],
                  if (slopes[i] > 0.002) "STILL RISING, not a plateau" else "flat"))
}
rising <- any(slopes > 0.002, na.rm = TRUE)
if (rising) {
  message("The curve has not plateaued by the last timepoint, so the plateau ",
          "method underestimates burst size. Report these as lower bounds, or ",
          "extend the time course.")
}

# ---- headline comparison ---------------------------------------------------

b   <- lapply(d$groups, burst_of)
res <- welch(b[[1L]], b[[2L]], sprintf("burst_%s_vs_%s", LABELS[1L], LABELS[2L]))
res$denominator    <- DENOMINATOR
res$ic0_assumed    <- IC0
res$plateau_pre    <- mean(colMeans(d$groups[[1]][tm %in% PRE_TIMES, , drop = FALSE]))
res$plateau_post   <- mean(colMeans(d$groups[[2]][tm %in% PRE_TIMES, , drop = FALSE]))
res$pre_times      <- paste(PRE_TIMES, collapse = "|")
res$plateau_times  <- paste(PLAT_TIMES, collapse = "|")
res$plateau_reached <- !rising

message(sprintf("burst %s = %.1f +/- %.1f;  %s = %.1f +/- %.1f;  ratio %.2f;  Welch p = %.4g",
                LABELS[1L], res$mean1, res$sd1, LABELS[2L], res$mean2, res$sd2,
                res$ratio, res$p))

# ---- USA300, same definition, against the Figure 2d value ------------------

usa_row <- NULL
usa_path <- file.path(data_dir, USA_FILE)
if (file.exists(usa_path)) {
  u  <- utils::read.delim(usa_path, check.names = FALSE, stringsAsFactors = FALSE)
  ut <- suppressWarnings(as.numeric(u[[1L]]))
  um <- as.matrix(u[, -1L, drop = FALSE]); storage.mode(um) <- "double"
  if (all(c(PRE_TIMES, PLAT_TIMES) %in% ut)) {
    pre_mean  <- colMeans(um[ut %in% PRE_TIMES, , drop = FALSE], na.rm = TRUE)
    plat_mean <- colMeans(um[ut %in% PLAT_TIMES, , drop = FALSE], na.rm = TRUE)
    denom_u <- if (DENOMINATOR == "measured") pre_mean else rep(IC0, length(pre_mean))
    bu <- (plat_mean - pre_mean) / denom_u
    usa_row <- data.frame(
      label = "burst_USA300", n1 = length(bu), n2 = NA_integer_,
      mean1 = mean(bu), mean2 = NA_real_, sd1 = stats::sd(bu), sd2 = NA_real_,
      t = NA_real_, df = NA_real_, p = NA_real_, ratio = NA_real_,
      denominator = DENOMINATOR, ic0_assumed = IC0,
      plateau_pre = mean(pre_mean), plateau_post = NA_real_,
      pre_times = paste(PRE_TIMES, collapse = "|"),
      plateau_times = paste(PLAT_TIMES, collapse = "|"),
      plateau_reached = !(tail_slope(um, min(PLAT_TIMES)) > 0.002),
      stringsAsFactors = FALSE
    )
    message(sprintf("USA300, same definition: %.1f +/- %.1f  (Figure 2d currently reports 99 +/- 18, which used the assumed IC0)",
                    mean(bu), stats::sd(bu)))
  } else {
    message("USA300 file lacks the configured timepoints; skipped")
  }
}

write_tsv(rbind(res, usa_row), file.path(results_dir, "burst_size_welch.tsv"))

# ---- sensitivity to the plateau window -------------------------------------

late <- tm[tm > max(PRE_TIMES)]
rows <- list()
for (i in seq_along(late)) for (j in seq_along(late)) {
  if (late[j] < late[i]) next
  win <- late[i:j]
  if (length(win) > 4L) next
  bb <- lapply(d$groups, function(m) {
    p0 <- colMeans(m[tm %in% PRE_TIMES, , drop = FALSE], na.rm = TRUE)
    p1 <- colMeans(m[tm %in% win, , drop = FALSE], na.rm = TRUE)
    (p1 - p0) / (if (DENOMINATOR == "measured") p0 else rep(IC0, length(p0)))
  })
  r <- welch(bb[[1L]], bb[[2L]], paste(win, collapse = "|"))
  r$window_start <- min(win); r$window_end <- max(win)
  rows[[length(rows) + 1L]] <- r
}
sweep <- do.call(rbind, rows)
write_tsv(sweep, file.path(results_dir, "burst_size_window_sweep.tsv"))
message("Across plateau windows: ", LABELS[1L], " burst ",
        sprintf("%.0f to %.0f", min(sweep$mean1), max(sweep$mean1)),
        "; ratio ", sprintf("%.2f to %.2f", min(sweep$ratio), max(sweep$ratio)),
        " — the ratio is stable, the absolute value is not.")

# ---- panel -----------------------------------------------------------------

p <- two_group_bar(
  b, LABELS,
  y_lab = "Burst size (PFU per infected cell)",
  fig_width = fig_width, fig_height = fig_height,
  subtitle = sprintf(
    "Plateau method, %s min vs %s min; Welch p = %.3g; ratio %.2f. Denominator: %s.%s",
    paste(PLAT_TIMES, collapse = "/"), paste(PRE_TIMES, collapse = "/"),
    res$p, res$ratio,
    if (DENOMINATOR == "measured") "measured pre-lysis plateau"
    else sprintf("assumed IC0 = %s", format(IC0, scientific = TRUE)),
    if (rising) " Curve still rising: lower bound." else ""),
  p_value = res$p
)

save_panel(p, file.path(results_dir, "burst_size"), fig_width, fig_height)
