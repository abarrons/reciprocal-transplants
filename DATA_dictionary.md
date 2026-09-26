# Data Dictionary

Input files required to run the reciprocal transplant microbiome analysis. These are the raw files read directly by the numbered reports. Everything under `results/rds/` is a regenerable intermediate produced by the pipeline (starting with `01_data_load.Rmd`) and is not listed here.

## Where inputs are read

| File | Read by | Purpose |
|------|---------|---------|
| `data/count_table.txt` | `01_data_load.Rmd` | Genome-level count table (whole community) |
| `data/GTDB_Genus_abundance_table.tsv` | `01_data_load.Rmd` | Genus-level abundance table (GTDB taxonomy) |
| `data/external/curto_core_genesDB/masterMD4DB.txt` | `01_data_load.Rmd` | Maps internal genome IDs to public IDs / genera |
| `data/external/curto_core_genesDB/gene_lengths_by_genome.csv` | `01_data_load.Rmd` | Total gene lengths per genome, for RPK normalization |
| `data/external/isolates/master_md_all_sequenced_Isolates_with_extendedMD.txt` | `01_data_load.Rmd` | Isolate metadata; source of clade / subclade taxonomy |
| `data/external/isolates/clade_colors.txt` | `01_data_load.Rmd` | Clade and subclade color palette |
| `data/external/RTs/RTs_sample_metadata.txt` | `01_data_load.Rmd` | Sample metadata for all reciprocal transplant samples |
| `data/external/mg/litter_chem/litter_chem.txt` | `01_data_load.Rmd` | Litter chemistry assay results |
| `data/external/mg/PRISM_ppt_tmin_tmean_tmax_30yr_normal_800m_monthly_normals.csv` | `04_convergence.Rmd` | PRISM 30-yr climate normals for reserve sites |
| `data/external/RTs/LOMA_PRISM_ppt_tmin_tmean_tmax_30yr_normal_800m_monthly_normals.csv` | `04_convergence.Rmd` | PRISM 30-yr climate normals for the LOMA reference site |
| `results/PERMANOVAS_Primer/PERMANOVAS_primer.csv` | `03_beta_div_permanova.Rmd` | PRIMER PERMANOVA summary (Curtobacterium) |
| `results/PERMANOVAS_Primer/PERMANOVAS_primer_genus.csv` | `03_beta_div_permanova.Rmd` | PRIMER PERMANOVA summary (genus level) |

`R/utils.R` sets `CA_CURTO_ROOT <- here::here("data", "external")`, so all `file.path(CA_CURTO_ROOT, ...)` paths resolve inside `data/external/`.

---

## File descriptions

### `data/count_table.txt`
Tab-separated genome-by-sample count matrix (the whole-community table).

- **Shape:** 470 genomes × ~325 samples (471 lines incl. header).
- **First column:** `genomeID` — internal genome identifier; remapped to public IDs via `masterMD4DB.txt` during loading.
- **Remaining columns:** one per sample, named with the transplant sample code (e.g. `L00S00I00-1`). Column-name dots are converted to hyphens on load.
- **Values:** raw read counts. Normalized to reads-per-kilobase (RPK) using `gene_lengths_by_genome.csv`.

### `data/GTDB_Genus_abundance_table.tsv`
Tab-separated genus-by-sample abundance matrix using GTDB taxonomy. Used to build the genus-level phyloseq object (`rt_g`).

- **Shape:** 28,990 genera × ~325 samples.
- **First column:** `genus` — genus name (row identifier).
- **Remaining columns:** one per sample; same sample-code naming as the count table.
- **Values:** genus abundances.

### `data/external/curto_core_genesDB/masterMD4DB.txt`
Whitespace-delimited database metadata table, no header. Maps genome IDs across naming schemes.

- **V1:** genus (e.g. `Curtobacterium`)
- **V2:** filename (e.g. `02D05.fna`)
- **V3:** genome ID matching `gene_lengths_by_genome.csv` columns
- **V4:** public genome ID used downstream (rownames of the count table become these)
- **V5:** internal genome ID matching `count_table.txt` `genomeID`
- **Shape:** 470 rows.

### `data/external/curto_core_genesDB/gene_lengths_by_genome.csv`
Comma-separated gene-length table used to normalize counts to RPK.

- **First column (unnamed, row names):** gene identifier.
- **Remaining columns:** one per genome, named by the `V3` genome ID (remapped to public IDs via `masterMD4DB.txt` on load).
- **Values:** gene lengths; column sums give each genome's total gene length.
- **Shape:** 1,075 genes × ~470 genomes.

### `data/external/isolates/master_md_all_sequenced_Isolates_with_extendedMD.txt`
Tab-separated isolate metadata with extended fields. Supplies the clade/subclade taxonomy joined into the taxonomy table.

- **Key columns used:** `Name` (genome ID join key), `clade`, `subclade` (columns 2, 3, 9 in the merge).
- **Other columns:** many additional isolate attributes (alternative names, seq IDs, filenames, reserve, isolation source/media, etc.).
- **Shape:** 608 isolates.

### `data/external/isolates/clade_colors.txt`
Tab-separated color mapping for clades and subclades.

- **`variable`:** clade or subclade label (e.g. `Clade.I`, `IA`, `IIB`).
- **`color`:** color value (name or hex).
- **Note:** the first 7 rows (clade-level `Clade.*` entries) are dropped on load; remaining rows define clade colors (`I`–`V`, `UND`) and subclade colors.
- **Shape:** 43 rows.

### `data/external/RTs/RTs_sample_metadata.txt`
Tab-separated sample metadata for the reciprocal transplant design. Row identifier is `original_label`.

- **`original_label`:** sample ID (matches count/abundance column names).
- **`Site`, `Reserve`, `Reserve_number`, `Reserve_ecosystem`:** site/reserve identifiers and ecosystem.
- **`Substrate_litter`, `Substrate_reserve`, `Litter_type`, `Litter_species`:** litter substrate origin and type.
- **`Inoculum`, `Inoculum_reserve`:** inoculum source (`IX` = no inoculum control).
- **`Bag_type`:** `Closed` or `Open`.
- **Derived on load:** reserve/substrate/inoculum code names (`reserve_code`, `s_code`, `i_code`), home/away flags (`r_home`, `s_home`, `i_home`), and composite `treatment` / `treatment_rep` / `treatment_2` labels.
- **Shape:** 332 samples.

### `data/external/mg/litter_chem/litter_chem.txt`
Tab-separated litter chemistry assay output (commercial forage panel). Heavily reshaped on load.

- **`desc_1`:** sample description; parsed into `sampleID`.
- **Columns used:** `ADF`, `NDF`, `CP` (crude protein), `Lignin`, `NFC`, `NSC`.
- **Derived:** `cellulose` (ADF − Lignin), `hemicellulose` (NDF − ADF), `Structural_Carbohydrates` (NFC − NSC); aggregated per reserve and used for a PCA of litter chemistry and Euclidean distances to LOMA.
- **Shape:** 95 rows (many extra assay columns are dropped).

### `data/external/mg/PRISM_ppt_tmin_tmean_tmax_30yr_normal_800m_monthly_normals.csv`
Comma-separated PRISM 30-year monthly climate normals for the reserve sites.

- **Columns:** `Name`, `Longitude`, `Latitude`, `Elevation (m)`, `Date` (month), `ppt (mm)`, `tmin (degrees C)`, `tmean (degrees C)`, `tmax (degrees C)`.
- **Shape:** 657 rows (12 months × sites).

### `data/external/RTs/LOMA_PRISM_ppt_tmin_tmean_tmax_30yr_normal_800m_monthly_normals.csv`
Same schema as the reserve PRISM file, but for the LOMA reference site only.

- **Shape:** 13 rows (12 monthly normals + header).

### `results/PERMANOVAS_Primer/PERMANOVAS_primer.csv`
Comma-separated PERMANOVA summary exported from PRIMER (Curtobacterium community). Read to build the variance-explained plots.

- **Columns:** `Source`, `df`, `SS`, `MS`, `Pseudo-F`, `P(perm)`, `% var. explained`, `pair`.
- **`Source`:** effect term (e.g. `Climate`, `Litter`).
- **`pair`:** reserve-vs-LOMA comparison label.
- **Note:** the header has a UTF-8 BOM on the first column.
- **Shape:** 48 rows.

### `results/PERMANOVAS_Primer/PERMANOVAS_primer_genus.csv`
Same schema as `PERMANOVAS_primer.csv`, for the genus-level PERMANOVA. Shape: 48 rows.

---

## Related files not read by the pipeline

- **`data/coreDB_Genus_abundance_table.tsv`** — a genus abundance table (coreDB taxonomy, ~3,136 genera) present in `data/` but not referenced by any of the numbered reports. Likely an alternate to the GTDB table.
- **`results/PERMANOVAS_Primer/PERMANOVA_*.txt` and `PERMANOVAS_primer.xlsx`** — raw PRIMER per-comparison outputs. Only the two `PERMANOVAS_primer*.csv` summaries are consumed by the reports.
