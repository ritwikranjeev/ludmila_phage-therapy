# Figure 5C: Kirby Bauer — zone diameter by strain, antibiotic, Ancestor vs Resistant
# Box-and-whisker plot; Wilcoxon test per antibiotic, bracket + star for significant differences.
# Aesthetic aligned with Figure 3 & 4 (theme_ludmila, Times New Roman, panel tag C).

if (!interactive()) pdf(NULL)

library(tidyverse)
library(rstatix)

# Figure dimensions: one row of 5 strain panels, wide layout
fig_width   <- 22
fig_height  <- 5.5
plot_margin <- margin(t = 12, r = 14, b = 10, l = 36, unit = "pt")
theme_ludmila <- function() {
  theme_bw(base_size = 13, base_family = "Times New Roman") +
  theme(
    axis.title       = element_text(face = "bold"),
    axis.title.x     = element_text(face = "bold", size = 14, margin = margin(t = 12, unit = "pt")),
    axis.title.y     = element_text(face = "bold", size = 14, margin = margin(r = 12, unit = "pt")),
    axis.text        = element_text(face = "bold"),
    axis.text.x      = element_text(face = "bold", margin = margin(t = 6, unit = "pt")),
    axis.text.y      = element_text(face = "bold", margin = margin(r = 6, unit = "pt")),
    panel.grid.minor = element_blank(),
    aspect.ratio     = fig_height / fig_width
  )
}

# Project root
project_dir <- Sys.getenv("PROJECT_DIR", "")
if (!nzchar(project_dir)) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[grepl("^--file=", args)]
  if (length(file_arg)) {
    script_path <- sub("^--file=", "", file_arg[1L])
    project_dir <- normalizePath(file.path(dirname(script_path), "..", "..", ".."), mustWork = FALSE)
  } else {
    project_dir <- normalizePath(file.path(getwd(), "..", "..", ".."), mustWork = FALSE)
  }
}
project_dir <- normalizePath(project_dir, mustWork = FALSE)

data_dir    <- file.path(project_dir, "data", "Figure5_Phage_resistant_Mutant_Analysis", "kirby_bauers")
results_dir <- file.path(project_dir, "analyses", "Figure5_Phage_resistant_Mutant_Analysis", "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 1) Read TSV ----
infile <- file.path(data_dir, "5c_KirbyBauer_Input.tsv")
df_wide <- readr::read_tsv(infile, show_col_types = FALSE)

# ---- 2) Pivot long ----
df <- df_wide %>%
  pivot_longer(
    cols = ends_with("_Diam"),
    names_to = "Antibiotic",
    values_to = "diameter_mm"
  ) %>%
  mutate(
    Antibiotic = str_replace(Antibiotic, "_Diam$", ""),
    Clone = factor(Clone, levels = c("Ancestor", "Resistant")),
    Antibiotic = factor(Antibiotic, levels = c("VA", "LEV", "MEM", "SXT", "MH")),
    Strain = factor(Strain, levels = c("USA300", "SA412", "OCPT02SA01", "OCPT04SA01", "PedCF"))
  )

df <- df %>% mutate(diameter_mm = na_if(diameter_mm, 0))

# ---- 3) Stats: Wilcoxon (Ancestor vs Resistant) within each Strain x Antibiotic ----
# Non-parametric; compares the two groups (Clone) for each antibiotic within each strain.
stats_df <- df %>%
  filter(!is.na(diameter_mm)) %>%
  group_by(Strain, Antibiotic) %>%
  wilcox_test(diameter_mm ~ Clone) %>%
  ungroup() %>%
  group_by(Strain) %>%
  mutate(p.adj = p.adjust(p, method = "BH")) %>%
  ungroup() %>%
  mutate(p_label = case_when(
    p.adj <= 0.001 ~ "***",
    p.adj <= 0.01  ~ "**",
    p.adj <= 0.05  ~ "*",
    TRUE           ~ "ns"
  ))

# ---- 4) Bracket positions: line from Ancestor to Resistant box, star above ----
# position_dodge(width = 0.75) places boxes at (x - 0.375) and (x + 0.375) for factor x.
y_pos <- df %>%
  group_by(Strain, Antibiotic) %>%
  summarise(y_max = max(diameter_mm, na.rm = TRUE), .groups = "drop") %>%
  mutate(y_max = if_else(is.finite(y_max), y_max, 0)) %>%
  mutate(y.position = y_max + 2) %>%
  select(Strain, Antibiotic, y.position)

# Annotation data: bracket from Ancestor to Resistant bar, star above (significant only)
stats_annot <- stats_df %>%
  left_join(y_pos, by = c("Strain", "Antibiotic")) %>%
  filter(p.adj <= 0.05) %>%
  mutate(
    x_min = as.numeric(Antibiotic) - 0.375,
    x_max = as.numeric(Antibiotic) + 0.375,
    label = p_label
  )

# ---- 5) Plot (box-and-whisker) ----
# Box: median, IQR; whiskers: 1.5*IQR (default); raw data overlaid as points (outliers hidden on box to avoid double-plotting)
pd <- position_dodge(width = 0.75)

p <- ggplot(df, aes(x = Antibiotic, y = diameter_mm, fill = Clone)) +
  geom_boxplot(
    position = pd,
    width = 0.65,
    outlier.shape = NA,
    coef = 1.5,
    color = "gray20",
    linewidth = 0.6
  ) +
  geom_point(
    data = df %>% filter(!is.na(diameter_mm)),
    position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75),
    size = 1.8,
    alpha = 0.7,
    color = "gray30"
  ) +
  # Bracket line (Ancestor bar to Resistant bar) + star for significant comparisons
  geom_segment(
    data = stats_annot,
    aes(x = x_min, xend = x_max, y = y.position, yend = y.position),
    inherit.aes = FALSE,
    linewidth = 0.35,
    color = "gray20"
  ) +
  geom_segment(
    data = stats_annot,
    aes(x = x_min, xend = x_min, y = y.position, yend = y.position - 0.8),
    inherit.aes = FALSE,
    linewidth = 0.35,
    color = "gray20"
  ) +
  geom_segment(
    data = stats_annot,
    aes(x = x_max, xend = x_max, y = y.position, yend = y.position - 0.8),
    inherit.aes = FALSE,
    linewidth = 0.35,
    color = "gray20"
  ) +
  geom_text(
    data = stats_annot,
    aes(x = (x_min + x_max) / 2, y = y.position, label = label),
    inherit.aes = FALSE,
    size = 4,
    fontface = "bold",
    vjust = -0.3
  ) +
  facet_wrap(~ Strain, ncol = 5, scales = "free_y") +
  scale_x_discrete(expand = expansion(add = 0.5)) +
  labs(
    x = "Antibiotic",
    y = "Zone of Inhibition (mm)",
    fill = "Clone",
    tag = "C"
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_fill_manual(
    values = c(
      "Ancestor" = "white",
      "Resistant" = "gray20"
    )
  ) +
  theme_ludmila() +
  theme(
    aspect.ratio     = NULL,
    axis.text.x      = element_text(angle = 45, hjust = 1),
    strip.text       = element_text(face = "bold"),
    legend.position  = "none",
    panel.spacing    = unit(24, "pt"),
    plot.margin      = plot_margin,
    plot.tag         = element_text(size = 18, face = "bold", margin = margin(r = 12, unit = "pt")),
    plot.tag.position = "topleft"
  ) +
  coord_cartesian(clip = "off")

p <- p + theme(
  axis.title.y      = element_text(face = "bold"),
  axis.title.y.left = element_text(face = "bold"),
  axis.text.y       = element_text(face = "bold"),
  axis.text.y.left  = element_text(face = "bold")
)

if (interactive()) print(p)
# One row of 5 strain panels, wider figure
ggsave(file.path(results_dir, "5c_KirbyBauer.png"), plot = p, width = fig_width, height = fig_height, units = "in", dpi = 300, bg = "white")
