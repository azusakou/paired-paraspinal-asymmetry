# Cross-etiology paraspinal metabolic asymmetry in scoliosis

Reproducible project for the study concept:

> Paired transcriptomes identify concave-side mitochondrial expression asymmetry
> across scoliosis etiologies

## Current status

The primary analysis now separates a cross-etiology concave-versus-convex
paraspinal response from AIS-associated systemic/mesenchymal differences. It
uses formal shared-effect and disease-by-side interaction contrasts rather than
intersections of separately thresholded differential-expression lists.

The plasma severity analysis remains intentionally gated. Phase 0 audits data availability, sample
annotations, clinical-label completeness, and overlap with published studies.
Prediction modelling is **not yet authorised by the data**: the public
PXD024246 files contain 84 raw runs and an 84-column protein-abundance table,
but the public accession XML does not map run IDs to clinical groups or expose
individual Cobb angles/covariates. The Wiley article page confirms the formal
supplement filename, but automated retrieval is blocked by publisher security
verification; the supplement has therefore not been misrepresented as local
data.

All analyses that do not require plasma group labels are complete: formal
paired bulk contrasts, leave-one-donor-out and exhaustive permutation
sensitivity analyses, donor-aware single-cell localization, BM-MSC comparison,
unsupervised plasma QC, and HPA/UniProt source discrimination. 

## Data

data will be here after paper published

## Environment restoration

R 4.4.2 and Bioconductor 3.20 dependencies are frozen in `renv.lock`. In a
clean R installation, restore them with:

```r
install.packages("renv")
renv::restore()
```

Python scripts use Python 3.9.6 and the standard library only. The compact
human-readable package summary is retained in `environment-lock-phase0.txt`.

## Phase-0 reproduction

```bash
python3 scripts/03_checksum_reused_inputs.py

python3 scripts/01_parse_pxd_proteins.py \
  --proteins-xml /Volumes/wang/paper16/ais_sc_plasma/raw/proteomics/PXD024246/quant/Proteins.xml \
  --px-xml /Volumes/wang/paper16/ais_sc_plasma/raw/proteomics/PXD024246/metadata/ProteomeXchange_PXD024246.xml \
  --outdir results/audit/pxd024246

python3 scripts/02_validate_phase0.py --project-root .
```

The parser uses only the Python standard library. It preserves SpreadsheetML
column indices, distinguishes blank cells from numerical zero, and does not
impute missing abundances.

## Phase-1 reproduction

```bash
Rscript scripts/10_build_sc_origin_reference.R
Rscript scripts/11_analyze_gse254300_bulk.R
Rscript scripts/12_proteomics_unsupervised_qc.R
Rscript scripts/13_build_tissue_support_matrix.R

python3 scripts/14_prepare_uniprot_requests.py
# Run each JSON payload through the bundled UniProt REST client documented in
# metadata/uniprot_requests/manifest.tsv, then:
Rscript scripts/15_integrate_external_origin_annotations.R
python3 scripts/16_validate_phase1.py
```

## Conserved-adaptation pivot reproduction

```bash
python3 scripts/20_parse_gse110359_platform.py \
  --raw-dir /path/to/gse110359 \
  --out-dir metadata/gse110359
Rscript scripts/21_analyze_gse110359_msc.R
Rscript scripts/22_analyze_orthogonal_evidence_axes.R
Rscript scripts/23_map_shared_pathways_to_celltypes.R
Rscript scripts/24_exact_paired_permutation_bulk.R
Rscript scripts/26_quantify_cross_etiology_equivalence.R
Rscript scripts/27_test_composition_and_precision_sensitivity.R
Rscript scripts/28_plan_validation_power.R
Rscript scripts/32_validate_model_dependence.R
python3 scripts/25_validate_pivot.py
Rscript scripts/31_build_nature_main_figures.R .
python3 scripts/29_build_source_data.py
```

The primary bulk analysis is produced by the updated
`scripts/11_analyze_gse254300_bulk.R`; rerun it before the pivot commands when
starting from raw counts.

## Analysis guardrails

- No group label is inferred from run-number order.
- No continuous Cobb-angle model is run without individual Cobb angles.
- Concave and convex libraries from the same donor are paired throughout.
- PRJNA722100 has only two AIS donors; cell-level counts are never treated as
  biological replicate counts.
- GSE254300 has five paired AIS donors and five paired congenital-scoliosis
  donors; it is a side-asymmetry dataset, not a healthy-vs-AIS cohort.
- Claims of paraspinal cell origin require expression, specificity, prevalence,
  donor consistency, secretion evidence, and exclusion of dominant liver/blood
  sources.

## Plasma prediction gate

Phase 1 may start differential proteomics only after an auditable mapping from
all 84 PXD024246 abundance columns to healthy/follow-up/surgery groups is
available. Nested cross-validation additionally requires an individual-level
outcome table and a stable subject identifier that keeps technical replicates
within folds.

The exact publisher supplement needed to resolve this gate is
`prca2169-sup-0006-TablesS1-S6.docx` (115 KB). A valid copy must be a ZIP-based
DOCX (`PK` magic bytes), not an HTML security-verification response.

