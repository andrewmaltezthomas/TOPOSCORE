# Custom scoring based on ecological topology of gut microbiota associated with cancer immunotherapy outcome

The gut microbiota influences clinical responses of cancer patients to immune checkpoint inhibitors. However, there is no consensus definition of detrimental dysbiosis. Based on metagenomics (MG) sequencing of 245 non-small cell lung cancer (NSCLC) patient feces, we constructed species-level co-abundance networks that were clustered into species interacting groups (SIG) correlating with overall survival. Thirty-seven and 45 MG species were associated with resistance (SIG1) and response (SIG2) to ICI, respectively. If combined with the quantification of Akkermansia species, this procedure allowed a person-based calculation of a topological score (TOPOSCORE) that was validated in additional 254 NSCLC patients and in 216 genitourinary cancers. Finally, this TOPOSCORE was translated into a 21 bacterial probe set-based qPCR-scoring that was validated in a prospective cohort of NSCLC patients, as well as in colorectal and melanoma patients. This approach could represent a dynamic diagnosis tool of intestinal dysbiosis to guide personalized microbiota-centered interventions.

[Original article](https://doi.org/10.1016/j.cell.2024.05.029)

If you find it useful, please cite:
Derosa L, Iebba V, Alves Costa Silva C, Piccinno G, Wu G, Lordello L, Routy B, et al. Custom scoring based on ecological topology of gut microbiota associated with cancer immunotherapy outcome. *Cell* (2024). https://doi.org/10.1016/j.cell.2024.05.029

# TOPOSCORE Calculator

A Shiny application that calculates TOPOSCORE from taxonomic profiles produced by several profilers and databases.

Web app: https://oncobiome-atlas.shinyapps.io/TOPOSCORE/

## Supported profilers and databases

| Option | Input | Identifiers used |
|---|---|---|
| Species - MetaPhlAn 4 vJan21 | `merge_metaphlan_tables.py` output | `s__` names (original TOPOSCORE species) |
| SGB - MetaPhlAn 4 vJan21 / vJun23 / vJan25 | `merge_metaphlan_tables.py` output with `t__SGB` rows | SGB IDs |
| Species - GTDB r207 / r220 | `clade_name` with `d__;p__;...;s__` lineages (sylph, mOTUs4, Kraken2/Bracken + GTDB) | GTDB species |
| Species - UHGG v2.0.2 | Kraken2/Bracken + UHGG, same lineage format | UHGG (GTDB-labelled) species |
| MSP - Meteor2 | `msp`, `msp_taxonomy` and sample columns | MSP IDs |

Taxa are in rows and samples in columns. Relative abundances can be percentages or fractions: fractions are detected and converted to percent. Annotation columns (`NCBI_tax_id`, `additional_species`, `bracken_taxid`, `genome_size`, ...) are ignored.

For Kraken2/Bracken, correct the Bracken read fractions for genome length before uploading. The *Akkermansia* High cut-off is defined on taxonomic (cell) abundance, and uncorrected read fractions overestimate it.

The non-MetaPhlAn options use a genome-based crosswalk (`app/data/toposcore_crosswalk.tsv`) that maps each of the 82 signature species to the identifier(s) used by each database. A species counts as detected if any of its identifiers has a non-zero abundance. Signature species that a table does not contain are listed after each run, because they always count as not detected.

## How the score is calculated

1. S-score = 0.5 × (SIG2 detected / 45 − SIG1 detected / 37 + 1), rounded to 3 decimals.
2. SIG class: SIG1 if S ≤ 0.5351, SIG2 if S ≥ 0.7911, otherwise grey zone.
3. *Akkermansia muciniphila*: Zero (0), High (≥ 4.799%; ≥ 4.52% for Meteor2), otherwise Normal.
4. TOPOSCORE: SIG1 → SIG1+, SIG2 → SIG2+; grey zone → SIG1+ if *Akkermansia* is Zero or High, SIG2+ if Normal.

## Using the app

1. Upload a merged profile table (tab-separated). The app recognises the table format and warns if it does not match the selected option.
2. Choose the profiler and database.
3. Press **Calculate TOPOSCORE**.

Results:
- Summary boxes: SIG1+ / SIG2+ counts and SIG class distribution.
- Run details: database, *Akkermansia* identifier and cut-off, abundance units, and the signature species missing from the table.
- **Results table**: per-sample S-score, SIG class, TOPOSCORE, *Akkermansia* status and abundance, and the detected SIG1/SIG2 species. Download it as CSV; the CSV includes the database and app version.
- **Score distribution**: your samples against the discovery cohort, split by overall survival. Highlight samples by name or by selecting table rows.
- **Signature species**: a heatmap of which SIG1/SIG2 species were detected in each sample, ordered by S-score, with species missing from the table marked.

The **Signature species** tab lists the crosswalk for the selected database.

Example input files are available in the sidebar.

## Running locally

Requires R ≥ 4.1.

```r
install.packages(c("shiny", "bslib", "DT", "ggplot2", "ggrepel", "tibble"))
shiny::runApp("app")
```

`source("app/dependencies.R")` lists any missing packages without installing them. For a reproducible setup, use `renv::init()` in the repository and commit `renv.lock`.

### Using the function without the app

```r
source("app/R/toposcore_calc.R")
res <- calculate_toposcore("my_profiles.tsv", taxonomy_database = "GTDB_r220")
attr(res, "diagnostics")   # species not found, Akkermansia identifier and cut-off, scaling
attr(res, "presence")      # samples x signature species, TRUE if detected
```

`toposcore_databases` lists the valid `taxonomy_database` values.

### Tests

```sh
cd app && Rscript tests/run_tests.R
```

### Deployment

```r
rsconnect::deployApp("app")
```

## Project structure

```
app/
├── app.R                      # Shiny UI and server
├── R/
│   ├── toposcore_calc.R       # calculation, input parsing, format detection
│   └── plots.R                # score distribution and signature heatmap
├── data/
│   ├── toposcore_crosswalk.tsv          # 82 signature species x database identifiers
│   ├── scores_disc.csv, clin_disc.csv   # discovery cohort reference
│   └── *_test_data.tsv                  # example inputs
├── www/                       # logos
├── tests/run_tests.R
└── dependencies.R
```

Uploaded files are processed in memory and not stored. For research use only.

## Contact

Valerio Iebba: valerio.iebba@gmail.com

Abhilash Dhal: adhalbiophysics@gmail.com
