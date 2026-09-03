# Supplementary Figure 2d — VIRIDIC intergenomic similarity heatmap, with a
# genome-length bar for each genome and OCPT-Ludmila's row and column outlined.
#
# Run with ../../run.sh, which supplies FIG_DATA and FIG_OUT.
#
# Rewritten from the original in three ways:
#   - paths come from the environment; the original assigned data_dir from
#     itself before it existed and could not run
#   - base R and ggplot2 only; the original called install.packages() at run time
#   - VIRIDIC's species and genus cluster assignments are carried into the
#     metadata output and summarised in the log, so the taxonomic claim in the
#     text has a file behind it
#
# Optional overrides:
#   S2_SIM_FILE      similarity matrix inside FIG_DATA
#   S2_CLUSTER_FILE  VIRIDIC cluster table inside FIG_DATA
#   S2_FASTA_FILE    FASTA used as VIRIDIC input inside FIG_DATA
#   S2_SUMMARY_FILE  NCBI data summary inside FIG_DATA
#   S2_FOCAL_ID      exact genome id of the focal phage

if (!interactive()) pdf(NULL)

suppressPackageStartupMessages({
  library(ggplot2)
})

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

SIM_FILE     <- env_or("S2_SIM_FILE",     "VIRIDIC_sim-dist_table.tsv")
CLUSTER_FILE <- env_or("S2_CLUSTER_FILE", "VIRIDIC_cluster_table.tsv")
FASTA_FILE   <- env_or("S2_FASTA_FILE",   "viridic_input_new.fasta")
SUMMARY_FILE <- env_or("S2_SUMMARY_FILE", "data_summary.tsv")
FOCAL_ID     <- env_or("S2_FOCAL_ID",     "OCPT-Ludmila-scaffold_1")

FONT_FAMILY <- "Times New Roman"
if (requireNamespace("systemfonts", quietly = TRUE)) {
  if (!FONT_FAMILY %in% unique(systemfonts::system_fonts()$family)) FONT_FAMILY <- "Times"
} else {
  FONT_FAMILY <- "Times"
}

path_of <- function(f) {
  p <- file.path(data_dir, f)
  if (!file.exists(p)) stop("Input file not found: ", p)
  p
}

# ---- similarity matrix -----------------------------------------------------

d   <- utils::read.delim(path_of(SIM_FILE), quote = "\"", check.names = FALSE,
                         stringsAsFactors = FALSE)
ids <- gsub("^\"|\"$", "", as.character(d[[1L]]))
m   <- as.matrix(d[, -1L, drop = FALSE])
storage.mode(m) <- "double"
colnames(m) <- gsub("^\"|\"$", "", colnames(m))
rownames(m) <- ids
if (!identical(rownames(m), colnames(m))) {
  stop("Similarity matrix is not square on the same ids")
}

# VIRIDIC writes similarity above and below the diagonal; they should agree.
asym <- max(abs(m - t(m)), na.rm = TRUE)
if (asym > 0.5) {
  message("WARNING: matrix halves disagree by up to ", signif(asym, 3),
          " — check that this file is the similarity table, not sim/dist combined")
} else if (asym > 0) {
  message("Matrix halves agree to ", signif(asym, 3), "; symmetrised by averaging")
}
m <- (m + t(m)) / 2
diag(m) <- 100

if (!FOCAL_ID %in% rownames(m)) {
  stop("Focal genome '", FOCAL_ID, "' not in the matrix. Present: ",
       paste(rownames(m), collapse = ", "))
}

# ---- names and lengths from the FASTA --------------------------------------

headers <- grep("^>", readLines(path_of(FASTA_FILE)), value = TRUE)
acc <- org <- character(0)
len <- integer(0)
for (h in headers) {
  rest  <- sub("^>", "", h)
  parts <- strsplit(rest, "\\s+", perl = TRUE)[[1L]]
  a     <- parts[1L]
  acc[a] <- a
  org[a] <- if (length(parts) > 1L) paste(parts[-1L], collapse = " ") else ""
  hit <- regmatches(rest, regexpr("length=(\\d+)", rest))
  if (length(hit)) len[a] <- as.integer(sub("length=", "", hit))
}

display_name <- function(a) {
  o <- unname(org[a])                     # NA when the id is absent from the FASTA
  if (is.na(o) || !nzchar(o)) return(a)
  short <- sub("^Staphylococcus (phage|virus) ", "", o, ignore.case = TRUE)
  short <- trimws(strsplit(short, ",")[[1L]][1L])
  if (a == FOCAL_ID) return("OCPT-Ludmila")
  if (nzchar(short)) short else a
}
id_to_display <- vapply(rownames(m), display_name, character(1))

# ---- genome length ---------------------------------------------------------

# NCBI's summary is keyed by organism name, the matrix by accession; match on
# the phage name that follows "phage " in both.
summ <- utils::read.delim(path_of(SUMMARY_FILE), check.names = FALSE,
                          stringsAsFactors = FALSE)
name_col <- grep("scientific", names(summ), ignore.case = TRUE, value = TRUE)[1L]
if (is.na(name_col)) name_col <- names(summ)[1L]
size_col <- names(summ)[tolower(names(summ)) == "size"][1L]
if (is.na(size_col)) size_col <- grep("size|length", names(summ), ignore.case = TRUE, value = TRUE)[1L]

key_of <- function(x) {
  x <- unname(x)
  if (is.na(x) || !nzchar(x)) return("")
  trimws(strsplit(sub(".*phage\\s+", "", x), ",")[[1L]][1L])
}

lengths_bp <- setNames(rep(NA_real_, nrow(m)), rownames(m))
for (a in rownames(m)) {
  if (a %in% names(len)) lengths_bp[[a]] <- len[[a]]
}
for (i in seq_len(nrow(summ))) {
  size <- suppressWarnings(as.numeric(summ[[size_col]][i]))
  if (is.na(size) || size <= 0) next
  ks <- key_of(as.character(summ[[name_col]][i]))
  for (a in rownames(m)) {
    if (startsWith(a, "OCPT") || !is.na(lengths_bp[[a]])) next
    kf <- key_of(org[a])
    if (nzchar(ks) && nzchar(kf) && (grepl(ks, kf, fixed = TRUE) || grepl(kf, ks, fixed = TRUE))) {
      lengths_bp[[a]] <- size
      break
    }
  }
}
missing_len <- names(lengths_bp)[is.na(lengths_bp)]
if (length(missing_len)) {
  message("No genome length matched for: ", paste(missing_len, collapse = ", "),
          " — their bars are drawn at zero")
}

# ---- VIRIDIC clusters ------------------------------------------------------

clust <- utils::read.delim(path_of(CLUSTER_FILE), check.names = FALSE,
                           stringsAsFactors = FALSE)
names(clust)[1L] <- "genome"
ci <- match(rownames(m), clust$genome)
species_cluster <- if ("species_cluster" %in% names(clust)) clust$species_cluster[ci] else NA
genus_cluster   <- if ("genus_cluster"   %in% names(clust)) clust$genus_cluster[ci]   else NA

focal_sp <- species_cluster[match(FOCAL_ID, rownames(m))]
focal_gn <- genus_cluster[match(FOCAL_ID, rownames(m))]
same_sp  <- setdiff(rownames(m)[!is.na(species_cluster) & species_cluster == focal_sp], FOCAL_ID)
same_gn  <- setdiff(rownames(m)[!is.na(genus_cluster)   & genus_cluster   == focal_gn], FOCAL_ID)
message("VIRIDIC species cluster ", focal_sp, ": ",
        if (length(same_sp)) paste(id_to_display[same_sp], collapse = ", ")
        else "OCPT-Ludmila alone (no genome at the 95% species threshold)")
message("VIRIDIC genus cluster ", focal_gn, ": ",
        if (length(same_gn)) paste(id_to_display[same_gn], collapse = ", ")
        else "OCPT-Ludmila alone (no genome at the 70% genus threshold)")

sim_to_focal <- m[FOCAL_ID, setdiff(rownames(m), FOCAL_ID)]
message("Highest similarity to OCPT-Ludmila: ",
        id_to_display[[names(which.max(sim_to_focal))]], " at ",
        signif(max(sim_to_focal), 4), "%")

# ---- order and layout ------------------------------------------------------

ord <- rownames(m)[stats::hclust(stats::as.dist(100 - m), method = "average")$order]
n   <- length(ord)
idx <- setNames(seq_len(n), ord)

grid <- expand.grid(row = seq_len(n), col = seq_len(n))
grid$value <- m[ord, ord][cbind(grid$row, grid$col)]
grid$x <- grid$col
grid$y <- n - grid$row + 1L

bar_start       <- n + 1.2
bar_width_cells <- 6.0
len_ord <- lengths_bp[ord]
ratio   <- len_ord / max(len_ord, na.rm = TRUE)
ratio[!is.finite(ratio)] <- 0

bars <- data.frame(
  genome_id = ord,
  y         = n - idx[ord] + 1L,
  xmin      = bar_start,
  xmax      = bar_start + ratio * bar_width_cells,
  is_focal  = ord == FOCAL_ID,
  stringsAsFactors = FALSE
)

focal_col <- match(FOCAL_ID, ord)
focal_row <- n - focal_col + 1L
outline <- data.frame(
  xmin = c(focal_col - 0.5, 0.5),
  xmax = c(focal_col + 0.5, n + 0.5),
  ymin = c(0.5, focal_row - 0.5),
  ymax = c(n + 0.5, focal_row + 0.5)
)

p <- ggplot() +
  geom_tile(data = grid, aes(x = x, y = y, fill = value)) +
  scale_fill_viridis_c(name = "Intergenomic\nsimilarity (%)", limits = c(0, 100),
                       na.value = "transparent", option = "cividis") +
  geom_segment(data = bars,
               aes(x = xmin, xend = xmax, y = y, yend = y, colour = is_focal),
               linewidth = 6, lineend = "butt", inherit.aes = FALSE) +
  scale_colour_manual(values = c("TRUE" = "black", "FALSE" = "gray70"), guide = "none") +
  geom_rect(data = outline,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, colour = "black", linewidth = 1.2, inherit.aes = FALSE) +
  scale_x_continuous(position = "top", breaks = seq_len(n),
                     labels = unname(id_to_display[ord]),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(breaks = seq_len(n),
                     labels = unname(id_to_display[rev(ord)]),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "Genome", y = "Genome") +
  theme_minimal(base_size = 10, base_family = FONT_FAMILY) +
  theme(
    text         = element_text(colour = "black", family = FONT_FAMILY),
    axis.text    = element_text(colour = "black", family = FONT_FAMILY),
    axis.text.x  = element_text(angle = 90, hjust = 0, vjust = 0.5),
    axis.title   = element_text(colour = "black", family = FONT_FAMILY),
    legend.text  = element_text(colour = "black", family = FONT_FAMILY),
    legend.title = element_text(colour = "black", family = FONT_FAMILY),
    panel.grid   = element_blank(),
    legend.position = "right"
  ) +
  coord_fixed(xlim = c(0.5, bar_start + bar_width_cells + 0.5), clip = "off")

fig_width  <- max(10, n * 0.35)
fig_height <- max(8,  n * 0.30)

out_base <- file.path(results_dir, "viridic_similarity_heatmap")
ggsave(paste0(out_base, ".png"), p, width = fig_width, height = fig_height,
       dpi = FIG_DPI, bg = "white")
pdf_ok <- save_pdf_best_effort(p, paste0(out_base, ".pdf"), width = fig_width,
                               height = fig_height, bg = "white")

meta_out <- data.frame(
  genome_id        = rownames(m),
  display_name     = unname(id_to_display),
  genome_length_bp = unname(lengths_bp),
  species_cluster  = species_cluster,
  genus_cluster    = genus_cluster,
  similarity_to_focal = unname(m[FOCAL_ID, ]),
  is_focal         = rownames(m) == FOCAL_ID,
  stringsAsFactors = FALSE
)
meta_out <- meta_out[order(-meta_out$similarity_to_focal), , drop = FALSE]
utils::write.table(meta_out, file.path(results_dir, "viridic_genome_metadata.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

if (pdf_ok) {
  message("Wrote: ", out_base, ".png, .pdf and viridic_genome_metadata.tsv")
} else {
  message("Wrote: ", out_base, ".png and viridic_genome_metadata.tsv",
          "   (no .pdf: cairo_pdf wrote no file)")
}
