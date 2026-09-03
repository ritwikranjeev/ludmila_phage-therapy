# Figure 3 – Stability Testing of Ludmila

Panels:
- **3a_pH_Stability** – pH stability (`scripts/pH.R`, `results/Ludmila_pH.png`)
- **3b_Temp_Stability** – Temperature stability (`scripts/temp.R`, `results/Ludmila_temp.png`)

Each panel has `scripts/`, `results/`, and optionally `config/`. Panel plots use shared aesthetics (8.5×6 in, Times New Roman, bold axes, consistent margins).

**Run panels:** From project root, run each script with `PROJECT_DIR` set, e.g. `Rscript analyses/Figure3_Stability_Testing/3a_pH_Stability/scripts/pH.R` (or set `PROJECT_DIR` to the project root).

**Build compiled figure:** From project root, run `Rscript manuscript/figures/build_fig3_fig4.R` to generate panel PNGs and `manuscript/figures/Figure3/Fig3.png`.
