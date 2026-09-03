# Supplementary Table 3 — Kirby-Bauer zone diameters: every comparison shown in
# Figure 6 and Supplementary Figure 6, with the full statistics behind them.
#
# Run with ../run.sh, which supplies FIG_DATA (Figure 6's data directory),
# FIG_OUT and FIG_SHARED (Figure 6's helpers).
#
# WHY THIS TABLE EXISTS
#
# Figure 6 draws brackets only where p <= 0.05, so a reader cannot tell a
# comparison that was tested and came out null from one that was never testable.
# This table reports all of them: group means and SDs, the test statistic, the
# nominal p, and Benjamini-Hochberg across the whole set, plus a reason for each
# comparison that could not be tested.
#
# The correction is applied across all comparisons in the table at once, which
# is the level a reader takes the figure in at. Correcting within each strain
# would be more generous and harder to justify: nothing in the design makes the
# five strains five separate experiments.
#
# Zone diameters are recorded CLSI-style, inclusive of the 6 mm disc, so a plate
# with no visible clearing reads 6 mm. Earlier versions of this sheet recorded
# no zone as 0 mm, which made "no zone" and "not measured" indistinguishable and
# required a zero policy to resolve; the re-measured sheet removes that ambiguity
# and the policy with it. check_disc_floor() enforces the convention.
#
# Optional overrides:
#   KB_DISC_MM, KB_TEST, KB_ABX_ORDER, KB_STRAIN_ORDER — as in Figure 6

need <- function(v) {
  x <- Sys.getenv(v, unset = NA_character_)
  if (is.na(x)) stop(v, " must be set. Run these tables with ./run.sh")
  x
}
data_dir    <- need("FIG_DATA")
results_dir <- need("FIG_OUT")
shared      <- need("FIG_SHARED")
source(shared)

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

DISC_MM      <- as.numeric(env_or("KB_DISC_MM", "6"))
TEST         <- env_or("KB_TEST", "welch")
ABX_ORDER    <- env_csv("KB_ABX_ORDER",    c("VA", "LEV", "MEM", "SXT", "MH", "LZD"))
STRAIN_ORDER <- env_csv("KB_STRAIN_ORDER", c("USA300", "SA412", "OCPT02SA01", "OCPT04SA01", "PedCF"))

mutants <- check_disc_floor(read_kirby_long(
  file.path(data_dir, "6a", env_or("KB_FILE_6A", "kirby_bauer_in_vitro_mutants.tsv")),
  group_col = "Clone", group_levels = c("Ancestor", "Resistant"),
  abx_order = ABX_ORDER, strain_order = STRAIN_ORDER), DISC_MM)

clinical <- check_disc_floor(read_kirby_long(
  file.path(data_dir, "6b", env_or("KB_FILE_6B", "kirby_bauer_pre_post.tsv")),
  group_col = "Clone", group_levels = c("Pre", "Post"),
  abx_order = ABX_ORDER, strain_order = NULL), DISC_MM)

a <- kirby_welch(mutants,  by_strain = TRUE,  test = TEST)
b <- kirby_welch(clinical, by_strain = FALSE, test = TEST)
b$Strain <- "Clinical"
common <- intersect(names(a), names(b))
out <- rbind(a[, common, drop = FALSE], b[, common, drop = FALSE])

out$fdr_bh_all <- stats::p.adjust(out$p, method = "BH")
out <- out[order(out$p, na.last = TRUE), , drop = FALSE]

utils::write.table(out, file.path(results_dir, "supp_table03_kirby.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

# ---- what the table says ---------------------------------------------------

n_tested <- sum(!is.na(out$p))
message("Supplementary Table 3: ", nrow(out), " comparisons, ", n_tested,
        " tested, test = ", TEST)

skipped <- out[is.na(out$p), , drop = FALSE]
if (nrow(skipped)) {
  message("  not testable (", nrow(skipped), "):")
  for (i in seq_len(nrow(skipped))) {
    message("    ", skipped$Strain[i], " ", skipped$Antibiotic[i], " — ", skipped$skipped[i])
  }
}

nominal <- out[!is.na(out$p) & out$p <= 0.05, , drop = FALSE]
message("  nominally p <= 0.05: ", nrow(nominal))
for (i in seq_len(nrow(nominal))) {
  message(sprintf("    %-12s %-4s p = %.4g   q = %.3g",
                  nominal$Strain[i], nominal$Antibiotic[i],
                  nominal$p[i], nominal$fdr_bh_all[i]))
}

surviving <- out[!is.na(out$fdr_bh_all) & out$fdr_bh_all < 0.05, , drop = FALSE]
if (nrow(surviving)) {
  message("  surviving BH across all ", n_tested, " comparisons (", nrow(surviving), "):")
  for (i in seq_len(nrow(surviving))) {
    message(sprintf("    %-12s %-4s  %.2f -> %.2f mm,  q = %.4g",
                    surviving$Strain[i], surviving$Antibiotic[i],
                    surviving$mean1[i], surviving$mean2[i], surviving$fdr_bh_all[i]))
  }
} else {
  message("  no comparison survives BH correction")
}
