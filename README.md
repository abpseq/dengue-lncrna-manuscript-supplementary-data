# Supplementary data tables (for manuscript)

This repository contains supplementary tables generated as part of the analysis described in the manuscript.

## Quick start

- All tables are available as:
  - tab-separated values (TSV) under `data/` (recommended for reproducibility)
  - a single Excel workbook: `data/Supplementary_Tables_All.xlsx` (reviewer-friendly)

## File index

See `docs/Supplementary_File_Index.tsv` for a one-line description of each supplementary table.

## Contents

### Differential expression result tables

Location: `data/differential_expression/`

- `Supplementary_Table_res_table_DI.tsv`
- `Supplementary_Table_res_table_DS.tsv`
- `Supplementary_Table_res_table_DWS.tsv`

These tables include (at minimum) identifiers/locus, group labels, effect size (`log2_fold_change`), test statistic, p-value, FDR-adjusted q-value, and a `significant` flag.

### Coding potential outputs

Location: `data/coding_potential/`

**CPC** (Coding Potential Calculator) tables:
- `Supplementary_Table_DI_CPC.tsv`
- `Supplementary_Table_DS_CPC.tsv`
- `Supplementary_Table_DWS_CPC.tsv`

**CPAT** (Coding Potential Assessment Tool) tables:
- `Supplementary_Table_DI_CPAT.tsv`
- `Supplementary_Table_DS_CPAT.tsv`
- `Supplementary_Table_DWS_CPAT.tsv`

## Notes

- Column names were normalized for consistency (trimmed whitespace, spaces converted to underscores, and minor punctuation cleaned).
- The biological meaning of group labels (e.g., DI/DS/DWS) match the manuscript definitions.

## Citation

## Soon
