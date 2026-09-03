# Analyses (by figure)

Analyses are organized by manuscript figure. Project root: `LudmilaCharacterization/`.

## Analyses structure

| Figure | Folder | Panels / content |
|--------|--------|-------------------|
| **Figure 2** | `Figure2_Genomic_Characterization/` | 2a Ludmila_Genome, 2b Ludmila_Phylogeny, 2c Ludmila_VIRIDIC_Analysis (each with `scripts/`, `results/`, `config/`) |
| **Figure 3** | `Figure3_Stability_Testing/` | 3a pH_Stability, 3b Temp_Stability (scripts: `pH.R`, `temp.R`; panel PNGs in `results/`) |
| **Figure 4** | `Figure4_Physiological_Characterization/` | 4a Ludmila_OSGC, 4b Ludmila_MOI (scripts: `osgc.R`, `moi-gc.R`; panel PNGs in `results/`) |
| **Figure 5** | `Figure5_Phage_resistant_Mutant_Analysis/` | Panel C: Kirby Bauer (`scripts/kirby.R`, `results/`) |
| **Figure 6** | `Figure6_Host_Range/` | Host-range: `scripts/`, `results/`, `config/` |

**Paths:** `analyses/FigureN_<title>/` then either `<panel>/scripts|results|config` (Fig 2–4), or `scripts/` and `results/` at figure level (Fig 5), or `scripts|results|config` (Fig 6).

## Data structure (figure-based)

Input data mirrors the figure layout under `data/`:

- **Figure 2** – `data/Figure2_Genomic_Characterization/` → 2a_Ludmila_Genome, 2b_Ludmila_Phylogeny, 2c_Ludmila_VIRIDIC
- **Figure 3** – `data/Figure3_Stability_Testing/` → 3a_pH_Stability, 3b_Temp_Stability
- **Figure 4** – `data/Figure4_Physiological_Characterization/` → 4a_Ludmila_OSGC, 4b_Ludmila_MOI
- **Figure 5** – `data/Figure5_Phage_resistant_Mutant_Analysis/` → kirby_bauers/
- **Figure 6** – `data/Figure6_Host_Range/` → hosts/ (faa, fna, gff), meta/, phage/ (faa, fna, gff)
