# Figure 5a: Mutation table visualization
# Rows = genes or events, columns = resistant mutants, cell color = mutation class.
# Usage: Rscript plot_figure5a.R --input <tsv> --output-png <path> --output-pdf <path> --output-csv <path>

if (!interactive()) pdf(NULL)

args <- commandArgs(trailingOnly = TRUE)
`%||%` <- function(a, b) if (is.null(a)) b else a
get_arg <- function(name) {
  i <- match(paste0("--", name), args)
  if (is.na(i) || i >= length(args)) return(NULL)
  args[i + 1L]
}
input_tsv   <- get_arg("input")       %||% stop("Need --input <tsv>")
output_png  <- get_arg("output-png")  %||% stop("Need --output-png <path>")
output_pdf  <- get_arg("output-pdf")  %||% stop("Need --output-pdf <path>")
output_csv  <- get_arg("output-csv")  %||% stop("Need --output-csv <path>")

# Read normalized mutation table (base R only, no package installs on cluster nodes)
df <- read.delim(input_tsv, sep = "\t", stringsAsFactors = FALSE, na.strings = c("", "NA", "."))
if (nrow(df) == 0) {
  df <- data.frame(mutant_id = character(), gene = character(), mutation_type = character(), stringsAsFactors = FALSE)
}

# Ensure required columns exist
for (c in c("mutant_id", "gene", "mutation_type")) {
  if (!c %in% names(df)) df[[c]] <- NA_character_
}
df$gene[is.na(df$gene) | df$gene == ""] <- "unknown"
df$mutation_type[is.na(df$mutation_type) | df$mutation_type == ""] <- "unknown"

genes <- sort(unique(df$gene))
mutants <- sort(unique(df$mutant_id))
if (length(genes) == 0) genes <- "no_genes"
if (length(mutants) == 0) mutants <- "no_mutants"

mat <- matrix("unknown", nrow = length(genes), ncol = length(mutants), dimnames = list(genes, mutants))
for (i in seq_len(nrow(df))) {
  g <- df$gene[i]
  m <- df$mutant_id[i]
  if (!is.na(g) && !is.na(m) && g %in% genes && m %in% mutants) {
    mat[g, m] <- df$mutation_type[i]
  }
}

base_colors <- c(
  "SNP" = "#2e86ab",
  "ins" = "#a23b72",
  "del" = "#f18f01",
  "MOB" = "#c73e1d",
  "JC" = "#6a994e",
  "MC" = "#bc4b51",
  "unknown" = "gray85"
)
all_types <- sort(unique(c(as.vector(mat), "unknown")))
extra_types <- setdiff(all_types, names(base_colors))
if (length(extra_types) > 0) {
  base_colors <- c(base_colors, setNames(rep("gray70", length(extra_types)), extra_types))
}
type_levels <- names(base_colors)
idx <- matrix(match(mat, type_levels), nrow = nrow(mat), ncol = ncol(mat))

draw_plot <- function() {
  par(mar = c(8, 10, 4, 10))
  image(
    x = seq_len(ncol(idx)),
    y = seq_len(nrow(idx)),
    z = t(idx[nrow(idx):1, , drop = FALSE]),
    col = unname(base_colors),
    xaxt = "n", yaxt = "n",
    xlab = "Resistant mutant", ylab = "Gene",
    main = "Figure 5a: Mutation table"
  )
  axis(1, at = seq_len(ncol(idx)), labels = colnames(mat), las = 2, cex.axis = 0.7)
  axis(2, at = seq_len(nrow(idx)), labels = rev(rownames(mat)), las = 2, cex.axis = 0.6)
  legend("topright", inset = c(-0.35, 0), legend = names(base_colors), fill = unname(base_colors), cex = 0.8, bty = "n")
}

dir.create(dirname(output_png), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)

png(output_png, width = max(2400, 90 * length(mutants)), height = max(1800, 30 * length(genes)), res = 300)
draw_plot()
dev.off()

pdf(output_pdf, width = max(8, 0.5 * length(mutants)), height = max(6, 0.25 * length(genes)))
draw_plot()
dev.off()

write.csv(df, output_csv, row.names = FALSE)
message("Wrote ", output_png, ", ", output_pdf, ", ", output_csv)
