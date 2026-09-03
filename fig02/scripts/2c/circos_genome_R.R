#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Ludmila phage genome – circular plot in R (run locally, not on HPC)
# Reads annotation from Pharokka (GFF or GenBank) and draws one contig + gene track.
# Install: install.packages("circlize"); BiocManager::install("rtracklayer")
# ------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1L && args[1L] %in% c("-h", "--help")) {
  cat("Usage: Rscript circos_genome_R.R [project_root] [output_dir]\n")
  cat("  project_root: path to 2a_Ludmila_Genome (default: current directory)\n")
  cat("  output_dir:   optional; folder for ludmila_circos.pdf (default: project results/figures)\n")
  cat("Output: results/figures/ludmila_circos.pdf or output_dir/ludmila_circos.pdf\n")
  quit(save = "no", status = 0)
}

# Project root: first argument, or directory containing this script, or getwd()
if (length(args) >= 1L && !args[1L] %in% c("-h", "--help")) {
  project_root <- normalizePath(args[1L], mustWork = TRUE)
} else {
  tryCatch({
    script_path <- if (exists("sys.frame") && exists("ofile")) {
      dirname(get("ofile", envir = sys.frame(1)))
    } else {
      dirname(sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())]))
    }
    if (length(script_path) && nzchar(script_path)) project_root <- normalizePath(file.path(script_path, ".."), mustWork = TRUE)
    else project_root <- getwd()
  }, error = function(e) project_root <- getwd())
}

# Optional output directory (e.g. your local folder for figures)
if (length(args) >= 2L) {
  out_dir <- normalizePath(args[2L], mustWork = FALSE)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  out_dir <- normalizePath(out_dir, mustWork = TRUE)
} else {
  out_dir <- file.path(project_root, "results", "figures")
}
out_pdf <- file.path(out_dir, "ludmila_circos.pdf")

# Paths under project (same layout as HPC pipeline; copy results/ from HPC if needed)
assembly_fasta  <- file.path(project_root, "results", "assembly", "Ludmila.fasta")
pharokka_dir    <- file.path(project_root, "results", "annotation", "pharokka")
annotation_gff  <- file.path(pharokka_dir, "pharokka.gff")
annotation_gbk  <- file.path(pharokka_dir, "pharokka.gbk")

# Find annotation file (GFF preferred, then GenBank)
annotation_file <- if (file.exists(annotation_gff)) annotation_gff else annotation_gbk
if (!file.exists(annotation_file)) {
  stop("Annotation not found. Copy results/annotation/pharokka/ from HPC, or run pipeline.\n",
       "  Looked for: ", annotation_gff, "\n  and: ", annotation_gbk)
}

# Optional: genome length from FASTA (for single contig); else from annotation
genome_length <- NA_integer_
contig_name   <- NA_character_
if (file.exists(assembly_fasta)) {
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    message("Install Biostrings for FASTA length: BiocManager::install('Biostrings')")
  } else {
    suppressPackageStartupMessages(library(Biostrings))
    fa <- Biostrings::readDNAStringSet(assembly_fasta)
    if (length(fa) > 0L) {
      genome_length <- width(fa)[1L]
      contig_name   <- names(fa)[1L]
    }
  }
}

# Load packages
suppressPackageStartupMessages({
  if (!requireNamespace("circlize", quietly = TRUE))
    stop("Install circlize: install.packages('circlize')")
  if (!requireNamespace("rtracklayer", quietly = TRUE))
    stop("Install rtracklayer: BiocManager::install('rtracklayer')")
  library(circlize)
  library(rtracklayer)
})

# Import features (GFF or GenBank)
gr <- rtracklayer::import(annotation_file)
if (length(gr) == 0L) stop("No features in ", annotation_file)

# Ensure we have one contig
seqlevels_in_use <- as.character(unique(seqnames(gr)))
if (length(seqlevels_in_use) > 1L) {
  message("Multiple contigs in annotation; using first: ", seqlevels_in_use[1L])
  gr <- gr[seqnames(gr) == seqlevels_in_use[1L]]
}
contig_name   <- if (is.na(contig_name)) seqlevels_in_use[1L] else contig_name
genome_length <- if (is.na(genome_length)) max(end(gr), na.rm = TRUE) else genome_length

# Restrict to gene/CDS-type features for the track (optional: keep all)
type_col <- NULL
if ("type" %in% names(mcols(gr))) type_col <- mcols(gr)$type
if (is.null(type_col) && "source" %in% names(mcols(gr))) type_col <- mcols(gr)$source
gene_like <- rep(TRUE, length(gr))
if (!is.null(type_col)) {
  type_col <- tolower(as.character(type_col))
  gene_like <- type_col %in% c("gene", "cds", "cds_from_genomic", "exon", "mrna")
}
gr_genes <- gr[gene_like]
if (length(gr_genes) == 0L) gr_genes <- gr

# BED-like data frame for circlize: chr, start, end, then optional columns
df_genes <- data.frame(
  chr    = as.character(seqnames(gr_genes)),
  start  = start(gr_genes),
  end    = end(gr_genes),
  strand = as.character(strand(gr_genes)),
  stringsAsFactors = FALSE
)
if ("type" %in% names(mcols(gr_genes))) df_genes$type <- as.character(mcols(gr_genes)$type)
if ("product" %in% names(mcols(gr_genes))) df_genes$product <- as.character(mcols(gr_genes)$product)

# Initialize circular layout with one sector (single contig)
circos_df <- data.frame(
  name  = contig_name,
  start = 0,
  end   = genome_length
)
ylim_genes <- c(0, 1)

# Write PDF (all drawing inside the PDF device)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
pdf(out_pdf, width = 8, height = 8)

circos.genomicInitialize(circos_df, plotType = NULL)
circos.track(
  ylim = c(0, 1),
  track.height = 0.08,
  panel.fun = function(x, y) {
    circos.axis(h = "top", labels.cex = 0.6)
    circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(2),
                CELL_META$sector.index, cex = 0.8, niceFacing = TRUE)
  },
  bg.border = NA
)
circos.genomicTrack(
  df_genes,
  ylim = ylim_genes,
  track.height = 0.25,
  panel.fun = function(region, value, ...) {
    strand <- value$strand
    y_bottom <- ifelse(strand == "+", 0.5, 0)
    y_top    <- ifelse(strand == "+", 1, 0.5)
    col_fwd  <- "#4A90D9"
    col_rev  <- "#E24A4A"
    col_vec  <- ifelse(strand == "+", col_fwd, col_rev)
    circos.genomicRect(region, value, ytop = y_top, ybottom = y_bottom,
                       col = col_vec, border = NA, ...)
  },
  bg.border = NA
)
circos.clear()

dev.off()
message("Wrote: ", out_pdf)
