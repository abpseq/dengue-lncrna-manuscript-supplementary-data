# Data dictionary

This document briefly describes columns present in the supplementary tables.

## Differential expression tables (`res_table_*.tsv`)

Typical columns:

- `test_id`: feature tested (e.g., transcript/gene locus ID)
- `gene_id`: feature identifier (often same as `test_id` in these tables)
- `gene`: gene symbol (may be `-` when not annotated)
- `locus`: genomic coordinates
- `sample_1`, `sample_2`: comparison labels (group names)
- `status`: tool status indicator (e.g., OK)
- `value_1`, `value_2`: group-level expression/abundance summaries as produced by the DE tool
- `log2_fold_change`: effect size (log2 scale) for `sample_2` vs `sample_1` (confirm with your Methods)
- `test_stat`: test statistic
- `p_value`: raw p-value
- `q_value`: multiple-testing adjusted p-value (FDR)
- `significant`: yes/no flag as provided by the tool

## CPC tables (`*_CPC.tsv`)

Typical columns:

- `id`: transcript ID (e.g., TCONS_*)
- `c_nc`: coding or noncoding label from CPC
- `coding_potential_score`: CPC score
- `evidence`: tool link/flag field from CPC output
- `utr_db_hits`, `rna_db_hits`: CPC link/flag fields from CPC output

## CPAT tables (`*_CPAT.tsv`)

Typical columns:

- `data_id`: row index in CPAT output
- `sequence_name`: transcript ID (e.g., TCONS_*)
- `rna_size`: transcript length
- `orf_size`: predicted ORF length
- `ficket_score`: CPAT feature score
- `hexamer_score`: CPAT feature score
- `coding_probability`: predicted coding probability
- `coding_label`: yes/no label
