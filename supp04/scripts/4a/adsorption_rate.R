# Supplementary Figure 4a — adsorption rate, pre- vs post-treatment isolate.
#
# Run with ../../run.sh, which supplies FIG_DATA, FIG_OUT and FIG_SHARED.
#
# METHOD, and why it differs from the original.
#
# Adsorption is first order in free phage, so ln(P) falls linearly with time
# while adsorption is under way:  ln(P_t) = ln(P_0) - k*N*t.
#
# The original analysis took two timepoints and solved that equation for k. That
# makes the answer depend entirely on which two points are chosen: across the
# possible pairs in this dataset the pre/post ratio ranges from 0.01 to 3.6 and
# the effect changes sign. This script instead fits the slope over the whole
# adsorption phase by least squares, one fit per replicate, and compares the
# slopes. It also writes the full two-point sweep so the instability is on
# record rather than hidden behind one chosen pair.
#
# TWO CONSTANTS THE DATA DOES NOT CONTAIN.
#
#   Host density N. The fitted slope is k*N, the observed rate constant, in
#   units of 1/min. Converting it to k in mL/min requires N, which was assumed
#   rather than measured. The pre/post RATIO does not involve N at all: it
#   cancels. So the comparison stands on its own; the absolute k does not.
#
#   The titre at t = 0. Sampling starts at 2 min, so phage adsorbed before the
#   first sample is invisible and every absolute rate here is a lower bound.
#
# Optional overrides:
#   S4_ROA_FILE     input file name inside FIG_DATA
#   S4_GROUPS       comma-separated column prefixes identifying the two groups
#   S4_LABELS       display labels for those groups
#   S4_HOST_DENSITY assumed host density in CFU/mL
#   S4_WINDOW       "start,end" in minutes for the adsorption phase

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

INFILE  <- env_or("S4_ROA_FILE", "roa_pre_post.tsv")
GROUPS  <- env_csv("S4_GROUPS", c("OCPT04SA01", "OCPT04SA03"))
LABELS  <- env_csv("S4_LABELS", c("Pre", "Post"))
HOST_N  <- as.numeric(env_or("S4_HOST_DENSITY", "1e8"))
WINDOW  <- env_num("S4_WINDOW", c(2, 6))

fig_width  <- 5.5
fig_height <- 5

d <- read_timecourse(file.path(data_dir, INFILE), GROUPS, LABELS)
tm <- d$time

message("Host density assumed: ", format(HOST_N, scientific = TRUE),
        " CFU/mL — CONFIRM AT BENCH. The pre/post ratio does not depend on it.")
message("First sample at ", min(tm), " min; phage adsorbed before then is not observed.")

# ---- is the chosen window still in the adsorption phase? -------------------

# Free phage relative to the first sample, per group. Once this stops falling,
# adsorption is over and any later point only adds flat tail to the fit.
frac <- lapply(d$groups, function(m) rowMeans(m, na.rm = TRUE) / mean(m[1L, ], na.rm = TRUE))
message("Free phage remaining, relative to ", min(tm), " min:")
for (i in seq_along(tm)) {
  message("  ", format(tm[i], width = 4), " min   ",
          paste(sprintf("%s %.3f", LABELS, vapply(frac, `[`, numeric(1), i)), collapse = "   "))
}

# ---- per-replicate log-linear fit ------------------------------------------

fit_slopes <- function(m, t0, t1) {
  sel <- tm >= t0 & tm <= t1
  if (sum(sel) < 3L) {
    message("NOTE: only ", sum(sel), " timepoints in the window; the fit is barely constrained")
  }
  out <- vapply(seq_len(ncol(m)), function(j) {
    y <- log(m[sel, j])
    f <- stats::lm(y ~ tm[sel])
    c(slope = unname(stats::coef(f)[2L]), r2 = summary(f)$r.squared)
  }, numeric(2))
  list(k_obs = -out["slope", ], r2 = out["r2", ])
}

fits  <- lapply(d$groups, fit_slopes, t0 = WINDOW[1L], t1 = WINDOW[2L])
k_obs <- lapply(fits, `[[`, "k_obs")          # 1/min, independent of host density
k_ml  <- lapply(k_obs, function(v) v / HOST_N) # mL/min, only as good as HOST_N

message("Window ", WINDOW[1L], "-", WINDOW[2L], " min; mean R^2 of the log-linear fit: ",
        paste(sprintf("%s %.3f", LABELS, vapply(fits, function(f) mean(f$r2), numeric(1))),
              collapse = ", "))

res <- welch(k_obs[[1L]], k_obs[[2L]],
             sprintf("k_obs_%s_vs_%s_%g_to_%g_min", LABELS[1L], LABELS[2L],
                     WINDOW[1L], WINDOW[2L]))
res$units       <- "per_min"
res$host_density <- HOST_N
res$k_mL_per_min_1 <- res$mean1 / HOST_N
res$k_mL_per_min_2 <- res$mean2 / HOST_N
write_tsv(res, file.path(results_dir, "adsorption_rate_welch.tsv"))

message(sprintf("k_obs %s = %.4f +/- %.4f /min;  %s = %.4f +/- %.4f /min",
                LABELS[1L], res$mean1, res$sd1, LABELS[2L], res$mean2, res$sd2))
message(sprintf("ratio %s/%s = %.2f, Welch p = %.4g", LABELS[1L], LABELS[2L],
                res$ratio, res$p))
message(sprintf("as k: %.3g and %.3g mL/min at N = %s",
                res$k_mL_per_min_1, res$k_mL_per_min_2, format(HOST_N, scientific = TRUE)))

# ---- sensitivity: every window, both methods -------------------------------

rows <- list()
for (a in seq_along(tm)) for (b in seq_along(tm)) {
  if (tm[b] <= tm[a]) next
  # two-point estimate, the original method
  two <- lapply(d$groups, function(m) {
    -log(m[b, ] / m[a, ]) / (tm[b] - tm[a])
  })
  r2p <- welch(two[[1L]], two[[2L]], sprintf("%g_to_%g", tm[a], tm[b]))
  r2p$method <- "two_point"
  rows[[length(rows) + 1L]] <- r2p
  # regression estimate, this script's method, where the window holds >= 3 points
  if (sum(tm >= tm[a] & tm <= tm[b]) >= 3L) {
    f <- lapply(d$groups, fit_slopes, t0 = tm[a], t1 = tm[b])
    rr <- welch(f[[1L]]$k_obs, f[[2L]]$k_obs, sprintf("%g_to_%g", tm[a], tm[b]))
    rr$method <- "regression"
    rr$mean_r2 <- mean(c(f[[1L]]$r2, f[[2L]]$r2))
    rows[[length(rows) + 1L]] <- rr
  }
}
sweep <- do.call(rbind, lapply(rows, function(r) {
  if (!"mean_r2" %in% names(r)) r$mean_r2 <- NA_real_
  r
}))
write_tsv(sweep, file.path(results_dir, "adsorption_rate_window_sweep.tsv"))

tp <- sweep[sweep$method == "two_point" & is.finite(sweep$ratio), ]
message("Two-point estimates across all ", nrow(tp), " time pairs: ratio ",
        sprintf("%.2f to %.2f", min(tp$ratio), max(tp$ratio)),
        ", p ", sprintf("%.2g to %.2g", min(tp$p, na.rm = TRUE), max(tp$p, na.rm = TRUE)),
        " — which is why a single pair is not a result.")

# ---- panel -----------------------------------------------------------------

p <- two_group_bar(
  k_ml, LABELS,
  y_lab = "Adsorption rate constant k (mL/min)",
  fig_width = fig_width, fig_height = fig_height,
  subtitle = sprintf(
    "Log-linear fit %g-%g min; Welch on k_obs, p = %.3g; ratio %.2f. k shown at assumed N = %s CFU/mL.",
    WINDOW[1L], WINDOW[2L], res$p, res$ratio, format(HOST_N, scientific = TRUE))
)

save_panel(p, file.path(results_dir, "adsorption_rate"), fig_width, fig_height)
