# Font preflight for the whole figure set.
#
#   Rscript _shared/check_fonts.R
#
# WHY THIS EXISTS
#
# The scripts ask for a serif face by name, and three different names are in use
# across the set: "Times New Roman", "Times" and "serif". None of the save paths
# errors when a name is unavailable — ggsave goes to cairo_pdf and to ragg/quartz
# for PNG, and all of those SUBSTITUTE a default face silently rather than
# stopping. So a missing font does not fail the run; it ships a manuscript whose
# panels are set in whatever the substitute happened to be, possibly differing
# from panel to panel.
#
# That is the failure this checks for. It reports whether each requested name
# resolves, and what actual font file each one resolves to, so that "they all
# resolve" can be distinguished from "they all resolve to the same face".
#
# The only path that genuinely errors is the pdf() fallback, reached when
# cairo_pdf is unavailable: base pdf() rejects an unknown family with
# "invalid font type". That is checked at the end.

requested <- c("Times New Roman", "Times", "serif")

cat("R           ", R.version.string, "\n")
cat("platform    ", R.version$platform, "\n\n")

# ---- 1. do the graphics packages exist -------------------------------------

for (pkg in c("systemfonts", "ragg", "ggplot2")) {
  cat(sprintf("%-12s %s\n", pkg,
              if (requireNamespace(pkg, quietly = TRUE))
                paste("present", as.character(utils::packageVersion(pkg)))
              else "MISSING"))
}
# capabilities("cairo") is a COMPILE-TIME flag. On macOS the cairo shared library
# is linked against XQuartz, so it can report TRUE and still fail to load at
# runtime. Worse, that failure is a WARNING, not an error: cairo_pdf() returns
# without writing a file, any tryCatch(error=) guard around it does not fire, and
# the calling script reports success for a PDF that does not exist. So test it by
# actually writing a file and checking the result.
cat("cairo compile flag ", capabilities("cairo"), "\n")
.cairo_ok <- local({
  f <- tempfile(fileext = ".pdf")
  warned <- NULL
  ok <- tryCatch(
    withCallingHandlers({
      grDevices::cairo_pdf(f, width = 3, height = 3); plot(1); grDevices::dev.off(); TRUE
    }, warning = function(w) { warned <<- conditionMessage(w); invokeRestart("muffleWarning") }),
    error = function(e) { warned <<- conditionMessage(e); FALSE })
  if (!is.null(grDevices::dev.list())) try(grDevices::dev.off(), silent = TRUE)
  wrote <- file.exists(f) && file.size(f) > 1000
  unlink(f)
  if (!is.null(warned)) cat("cairo runtime       FAILED: ", warned, "\n", sep = "")
  cat("cairo actually wrote a PDF: ", wrote, "\n", sep = "")
  if (!wrote) {
    cat("\n  *** No PDF output. Every .pdf in this figure set will be silently\n")
    cat("      skipped while the scripts report success. Install XQuartz:\n")
    cat("        brew install --cask xquartz\n")
    cat("      then start a NEW R session and re-run this check.\n")
    cat("      /opt/X11 present: ", dir.exists("/opt/X11"), "\n\n", sep = "")
  }
  isTRUE(ok) && wrote
})
cat("\n")

if (!requireNamespace("systemfonts", quietly = TRUE)) {
  cat("Install systemfonts to resolve font names: install.packages(\"systemfonts\")\n")
  quit(status = 1)
}

# ---- 2. what does each requested name actually resolve to -------------------

fams <- unique(systemfonts::system_fonts()$family)
cat("Fonts visible to the graphics device:", length(fams), "families\n\n")

cat("REQUESTED NAME       LISTED?  RESOLVES TO\n")
paths <- character(0)
for (nm in requested) {
  m <- systemfonts::match_fonts(nm)          # never fails; falls back silently
  file <- if (is.data.frame(m)) m$path[1] else as.character(m)[1]
  paths <- c(paths, file)
  cat(sprintf("%-20s %-8s %s\n", nm,
              if (nm %in% fams) "yes" else "NO",
              basename(file)))
}
cat("\n")

if (length(unique(paths)) == 1L) {
  cat("All three names resolve to the same font file — typography will be consistent.\n")
} else {
  cat("WARNING: the three names resolve to DIFFERENT files. Panels using different\n")
  cat("names will be set in different faces. Unify them before the final render.\n")
}

missing <- requested[!requested %in% fams]
if (length(missing)) {
  cat("\nNot registered: ", paste(missing, collapse = ", "), "\n", sep = "")
  cat("On macOS, Times New Roman normally ships at\n")
  cat("  /System/Library/Fonts/Supplemental/Times New Roman.ttf\n")
  cat("If that file exists but the name is not listed above, register it:\n")
  cat("  systemfonts::register_font(\"Times New Roman\",\n")
  cat("    plain = \"/System/Library/Fonts/Supplemental/Times New Roman.ttf\")\n")
  f <- "/System/Library/Fonts/Supplemental/Times New Roman.ttf"
  cat("  file present on this machine: ", file.exists(f), "\n", sep = "")
}

# ---- 3. does the pdf() fallback accept these names --------------------------

cat("\npdf() fallback (only reached if cairo_pdf is unavailable):\n")
for (nm in requested) {
  tmp <- tempfile(fileext = ".pdf")
  ok <- tryCatch({
    grDevices::pdf(tmp, family = nm); plot(1); grDevices::dev.off(); TRUE
  }, error = function(e) {
    if (!is.null(grDevices::dev.list())) try(grDevices::dev.off(), silent = TRUE)
    conditionMessage(e)
  })
  unlink(tmp)
  cat(sprintf("  %-20s %s\n", nm, if (isTRUE(ok)) "accepted" else paste("ERROR:", ok)))
}

# ---- 4. what actually lands in a rendered file ------------------------------

if (requireNamespace("ggplot2", quietly = TRUE)) {
  cat("\nRendering a test panel in each name -> _shared/font_check_<name>.png\n")
  cat("Open them and compare: if a name silently substituted, it will look different.\n")
  outdir <- dirname(sub("^--file=", "",
                        grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
  if (is.na(outdir) || !nzchar(outdir)) outdir <- "."
  for (nm in requested) {
    p <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
      ggplot2::geom_text(label = paste0(nm, ": Staphylococcus aureus 1234567890"),
                         family = nm, size = 5) +
      ggplot2::theme_bw(base_size = 13, base_family = nm) +
      ggplot2::labs(title = nm, x = "zone diameter (mm)", y = "PFU/mL")
    f <- file.path(outdir, paste0("font_check_", gsub("[^A-Za-z]", "", nm), ".png"))
    suppressWarnings(ggplot2::ggsave(f, p, width = 6, height = 3, dpi = 150, bg = "white"))
    cat("  wrote ", basename(f), "\n", sep = "")
  }
}
