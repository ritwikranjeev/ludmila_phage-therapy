# Figure 4 – Physiological Characterization of Ludmila

Panels:
- **4a_Ludmila_OSGC** – One-step growth curve (`scripts/osgc.R`, `results/Ludmila_OSGC.png`)
- **4b_Ludmila_MOI** – MOI curves (`scripts/moi-gc.R`, `results/Ludmila_MOI_OCPT04SA01.png`, `Ludmila_MOI_SA412.png`)

Each panel has `scripts/`, `results/`, and optionally `config/`. Panel plots use shared aesthetics (8.5×6 in, Times New Roman, bold axes, consistent margins).

**Run panels:** From project root, run each script with `PROJECT_DIR` set (or the build script sets it).

**Build compiled figure:** From project root, run `Rscript manuscript/figures/build_fig3_fig4.R` to generate panel PNGs and `manuscript/figures/Figure4/Fig4.png`.
