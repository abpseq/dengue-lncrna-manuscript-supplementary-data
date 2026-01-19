# Data dictionary (supplementary tables)

## Differential expression result tables (`res_table_*`)

Typical columns:

- `test_id`, `gene_id`: transcript/locus identifier used in the analysis
- `gene`: gene symbol (may be `-` if not annotated)
- `locus`: genomic locus (e.g., `chr:start-end`)
- `sample_1`, `sample_2`: compared groups/conditions
- `value_1`, `value_2`: normalized expression/abundance values used by the DE tool
- `log2_fold_change`: log2(value_2 / value_1)
- `test_stat`: test statistic reported by the DE method
- `p_value`: nominal p-value
- `q_value`: multiple-testing adjusted p-value (FDR)
- `significant`: yes/no flag based on the manuscript threshold

## CPC tables (`*_CPC`)

- `ID`: transcript identifier
- `C/NC`: coding/noncoding call
- `CODING_POTENTIAL_SCORE`: CPC score
- Remaining columns are CPC report links/fields as produced by CPC

## CPAT tables (`*_CPAT`)

- `Data_ID`: row number in CPAT output
- `Sequence_Name`: transcript identifier
- `RNA_Size`, `ORF_Size`
- `Ficket_Score`, `Hexamer_Score`
- `Coding_Probability`
- `Coding_Label`: yes/no call
