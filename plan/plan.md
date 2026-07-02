# Directory Standardization Plan

## Goal

Impose a consistent folder skeleton across all 21 research experiment directories in this
repository, without renaming the top-level directories or altering supplementary-file
naming conventions tied to published papers.

---

## Target Standard Structure

Each experiment directory should contain:

```
ExperimentDir/
├── code/          # R scripts, Rmd notebooks, shell scripts
├── datasets/      # raw + processed data (RDS, RData, CSV, FASTA, TXT, XLSX, BIOM, etc.)
├── functions/     # helper R functions sourced by scripts in code/
├── misc/          # scratch analyses, old versions, loose figures, miscellaneous outputs
├── README.md      # project description, data sources, how to run the code
└── LICENSE        # MIT license
```

**Optional folders** (kept where they already exist, not created elsewhere):
- `results/` — statistical outputs, model summaries, tables
- `figures/` — publication-ready figures

**Project-specific folders in Liu_etal_2026** (left in place, not removed):
`phylogeny/`, `sequence_alignments/`, `HPCC/`, `SRA/`, `renv/`, `demo_output/`

---

## Current State of All 21 Directories

### Category A — Partially structured (folders exist but gaps remain)

| Directory | Missing folders | Missing LICENSE |
|---|---|---|
| Benucci_etal_2018_FijiSoilMicrobiome | `misc/` | no |
| Benucci_etal_2019_Lycopods | `functions/`, `misc/` | yes |
| Benucci_etal_2019_MorchellaMicrobiome | `misc/` | yes |
| Benucci_etal_2025_CompetitionReleaseTruffle | `misc/` | yes |
| Garcia-Barreda_et_al_TruffleNestBacteria | `functions/`, `misc/` | no |
| Kelley_etal_2024_ArabidopsisDrought | `functions/`, `misc/` | yes |
| Liu_etal_2026_AMF_FertUnfertSwitchgrass | already well-structured | yes |

**Action:** create missing folders + add MIT `LICENSE` where absent. No file moves.

### Category B — Flat structure (all files at root, need sorting)

| Directory | Notes |
|---|---|
| Bell-Dereske_etal_2023_FungiBactSwitchgrass | 2 RData + 1 R script |
| Bellucci_etal_2025_IsopreneEffectOnMicrobiome | 2 RDS + 1 R script + .Rproj |
| Benucci_etal_2022_DeepSoilMicrobiome | 2 RDS + 1 RData + 1 R script |
| Benucci_etal_2023_LeafDecompSoilMicrocosm | 4 RDS + 1 R script |
| daCosta_and_Benucci_etal_2021_SwithgrassMicrobiome | 2 RDS + 2 FASTA + 1 R script |
| Haan_etal_2022_BiodiversityKBS | 10 FASTA + 1 R script + 1 shell script |
| Liu_etal_2024_TransgenicSwitchgrassMicrobiome | 2 RData + 1 R script |
| Longley_etal_2023_CoralMicrobiomeBleaching | 1 R script |
| Marozzi_etal_2022_WhiteTruffleMicrobiome | OTU tables, FASTA, metadata, 1 R script |
| Martin-Pinto_etal_2022_CistusFire | OTU tables, metadata CSV, FASTA, 1 R script |
| Noel_etal_2022_FungicidePulseDisturbance | 1 R script |
| VanWallendael_etal_2021_SwitchgrassLeafFungalMicrobiome | 3 helper R files + 1 main R script |

**Action:** create skeleton + move loose files per the rules below.

### Category C — Supplementary-file naming (leave flat)

| Directory | Notes |
|---|---|
| Benucci_etal_2020_Mao-Tofu_Microbiome | Files named `Suppl_File_X_*` — matches published paper |
| Benucci_etal_2020_Patient_Propagules | Files named `Suppl_File_X_*` — matches published paper |

**Action:** create skeleton folders with `.gitkeep` only. No file moves.

---

## File Sorting Rules

| Extension(s) | Destination |
|---|---|
| `.R`, `.r`, `.Rmd`, `.rmd`, `.sh` | `code/` |
| `.rds`, `.RDS`, `.RData`, `.Rdata`, `.rda` | `datasets/` |
| `.csv`, `.tsv`, `.txt`, `.xlsx`, `.xls`, `.biom` | `datasets/` |
| `.fasta`, `.fa`, `.fastq` | `datasets/` |
| `.pdf`, `.eps`, `.png`, `.jpg`, `.jpeg` | `misc/` |
| `.Rproj`, `README.md`, `LICENSE*`, `renv.lock` | root (never moved) |

Files with unrecognised extensions are flagged as `[skip]` in the log — nothing is
silently lost.

---

## How to Run

```bash
# 1. Preview everything first (no changes made)
bash plan/standardize_dirs.sh --dry-run

# 2. Review the output, then apply
bash plan/standardize_dirs.sh
```

The script prints a log of every action: `[mkdir]`, `[mv]`, `[LICENSE]`, `[skip]`.

---

## Post-run Verification

```bash
find /home/gian/Dropbox/6_PROJETCS/Published-R-Code \
  -maxdepth 2 -type d -not -path '*/\.*' | sort
```

Every experiment directory should now show `code/`, `datasets/`, `functions/`, `misc/`
alongside `README.md` and `LICENSE`.

---

## .gitignore additions (to apply after the script)

Add to the repo root `.gitignore`:

```
# Generated output — not tracked
**/results/
**/figures/
```
