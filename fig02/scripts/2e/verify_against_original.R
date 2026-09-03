# Figure 2e — confirm the layout-driven script reproduces the original hardcoded one.
#
# Recomputes the plotted series two ways from the same plate-reader file:
#   OLD  wells chosen by hardcoded column indices, as in moi-gc.R
#   NEW  wells chosen by reading plate_layout.tsv, as in moi_curves.R
# and compares them. Base R only, so it runs even if tidyverse is unavailable.
#
# Exits non-zero if the two disagree.

data_dir <- Sys.getenv("FIG_DATA", unset = NA_character_)
if (is.na(data_dir)) stop("FIG_DATA must be set. Run via ../../run.sh or set it manually.")

raw    <- read.delim(file.path(data_dir, Sys.getenv("MOI_DATA",   unset = "Ludmila_MOI_GC_Data.txt")), check.names = FALSE)
layout <- read.delim(file.path(data_dir, Sys.getenv("MOI_LAYOUT", unset = "plate_layout.tsv")),        check.names = FALSE)

raw <- raw[nzchar(trimws(as.character(raw$Time))), , drop = FALSE]
hours <- vapply(strsplit(as.character(raw$Time), ":"), function(p) {
  p <- as.numeric(p); (p[1] * 3600 + p[2] * 60 + p[3]) / 3600
}, numeric(1))

mean_sd <- function(v) {
  v <- v[!is.na(v)]
  c(avg = mean(v), sd = if (length(v) > 1) stats::sd(v) else 0, n = length(v))
}
wells_of <- function(nm) as.numeric(unlist(raw[i, nm, drop = TRUE]))

STRAIN   <- "USA300"
SHOW     <- c("1", "0.1", "0.01", "0.001", "Bacteria Only", "Blank")
MAX_HOURS <- 16
MOI_BY_ROW <- c(A = "1", B = "0.1", C = "0.01", D = "0.001",
                E = "0.0001", F = "0.00001", G = "0.000001", H = "Bacteria Only")

old <- new <- list()
for (i in seq_len(nrow(raw))) {
  h <- round(hours[i], 6)

  # ---- OLD: hardcoded columns 4:6, H4/H5 for bacteria-only, H3/H6 blank ----
  b_old <- mean_sd(wells_of(c("H3", "H6")))[["avg"]]
  for (r in names(MOI_BY_ROW)) {
    w <- if (r == "H") c("H4", "H5") else paste0(r, 4:6)
    m <- mean_sd(wells_of(w))
    old[[paste(h, MOI_BY_ROW[[r]])]] <- c(m[["avg"]] - b_old, m[["sd"]], m[["n"]])
  }
  bs <- mean_sd(wells_of(c("H3", "H6")))
  old[[paste(h, "Blank")]] <- c(0, bs[["sd"]], bs[["n"]])

  # ---- NEW: wells taken from plate_layout.tsv ----
  blank_w <- layout$well[layout$role == "blank"]
  b_new <- mean_sd(wells_of(blank_w))[["avg"]]
  keep <- layout$strain == STRAIN | layout$role == "blank"
  sub  <- layout[keep, , drop = FALSE]
  series <- ifelse(sub$role == "blank", "Blank",
             ifelse(sub$moi == "bacteria_only", "Bacteria Only", sub$moi))
  for (s in unique(series)) {
    m <- mean_sd(wells_of(sub$well[series == s]))
    new[[paste(h, s)]] <- c(m[["avg"]] - b_new, m[["sd"]], m[["n"]])
  }
}

keys <- names(old)
keys <- keys[vapply(strsplit(keys, " ", fixed = TRUE),
                    function(p) as.numeric(p[1]) <= MAX_HOURS, logical(1))]
keys <- keys[vapply(keys, function(k) sub("^[^ ]+ ", "", k) %in% SHOW, logical(1))]

miss <- setdiff(keys, names(new))
cmp  <- intersect(keys, names(new))
d_avg <- vapply(cmp, function(k) abs(old[[k]][1] - new[[k]][1]), numeric(1))
d_sd  <- vapply(cmp, function(k) abs(old[[k]][2] - new[[k]][2]), numeric(1))
d_n   <- vapply(cmp, function(k) abs(old[[k]][3] - new[[k]][3]), numeric(1))

cat("\n  compared", length(cmp), "(timepoint x series) values, Hours <=", MAX_HOURS, "\n")
cat("  max difference in mean :", format(max(d_avg), scientific = TRUE), "\n")
cat("  max difference in SD   :", format(max(d_sd),  scientific = TRUE), "\n")
cat("  replicate-count mismatches:", sum(d_n > 0), "\n")
cat("  series missing from new  :", length(miss), "\n\n")

TOL <- 1e-10
if (length(miss) || max(d_avg) > TOL || max(d_sd) > TOL || any(d_n > 0)) {
  bad <- cmp[d_avg > TOL | d_sd > TOL | d_n > 0]
  for (k in utils::head(bad, 5))
    cat("    MISMATCH", k, " old=", old[[k]], " new=", new[[k]], "\n")
  stop("moi_curves.R does not reproduce moi-gc.R")
}
cat("  IDENTICAL - moi_curves.R reproduces moi-gc.R exactly\n")
