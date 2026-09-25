# Reciprocal Transplant Microbiome Analysis

R Markdown analyses of microbial community composition, alpha and beta diversity, PERMANOVA, convergence, and NMDS for the reciprocal transplant study.

## Project contents

- `01_data_load.Rmd` loads and normalizes the input tables and prepares intermediate data.
- `02_alpha_div.Rmd` analyzes alpha diversity.
- `03_beta_div_permanova.Rmd` analyzes beta diversity and PERMANOVA.
- `04_convergence.Rmd` analyzes community convergence.
- `05_nmds_plots.Rmd` generates NMDS figures.
- `06_combined_figures.Rmd` assembles combined figures.
- `R/utils.R` contains shared settings and helper functions.
- `data/` contains the count and abundance tables plus the metadata and reference inputs listed below.
- `results/PERMANOVAS_Primer/` contains two PRIMER summary tables used by the beta-diversity report.
- `plots/` and `results/` contain figures, tables, and analysis outputs. `results/rds/` contains regenerable intermediate files and is not tracked.
- `RTs.Rmd` is the original, consolidated analysis report. Its inputs are project-local, but some output paths remain machine-specific; use the numbered reports for the staged workflow.

## Requirements

Use R and RStudio, with the R packages imported by the analysis reports installed. All file inputs are included in the repository; no sibling `CA_curto` directory or `CA_CURTO_ROOT` environment variable is needed. The additional inputs bundled from the surrounding study folders are:

- `data/external/curto_core_genesDB/masterMD4DB.txt`
- `data/external/curto_core_genesDB/gene_lengths_by_genome.csv`
- `data/external/isolates/master_md_all_sequenced_Isolates_with_extendedMD.txt`
- `data/external/isolates/clade_colors.txt`
- `data/external/RTs/RTs_sample_metadata.txt`
- `data/external/RTs/LOMA_PRISM_ppt_tmin_tmean_tmax_30yr_normal_800m_monthly_normals.csv`
- `data/external/mg/PRISM_ppt_tmin_tmean_tmax_30yr_normal_800m_monthly_normals.csv`
- `data/external/mg/litter_chem/litter_chem.txt`
- `results/PERMANOVAS_Primer/PERMANOVAS_primer.csv`
- `results/PERMANOVAS_Primer/PERMANOVAS_primer_genus.csv`

Confirm that you have permission to publish and redistribute the bundled study data.

## Run the analysis

Open `RTs_analysis.Rproj` in RStudio. Render the numbered reports in order, starting with `01_data_load.Rmd`; later reports use intermediate RDS files produced by earlier steps. Render a report with `rmarkdown::render("01_data_load.Rmd")`, replacing the filename for each subsequent report. Outputs are written under the project’s `results/` and `plots/` directories.
