#!/usr/bin/env Rscript
# =============================================================================
# 01_filter_novel_lnc_results.R
# -----------------------------------------------------------------------------
# Purpose
#   - Read 3 differential analysis result tables (DI/DWS/DS or A/B/C)
#   - Apply significance filters (p, q, |log2FC|)
#   - Export selected lncRNA IDs and (optionally) subset circos input files
#
# Why this exists
#   Your legacy script repeated the same filtering blocks many times with
#   hard-coded paths. This refactor makes it reproducible and CLI-driven.
#
# Inputs
#   - Result tables (TSV) with at least: gene_id, p_value, q_value, log2FC column
#   - Optional circos input tables containing a gene/ID column
#
# Outputs
#   - out_dir/summary_counts.tsv
#   - out_dir/selected_ids/<group>_sig_q_lfc_ids.txt
#   - out_dir/selected_tables/<group>_sig_q_lfc_table.tsv
#   - Optional: out_dir/circos_subset/<group>_circos_subset.tsv
#
# Example
#   Rscript code/01_filter_novel_lnc_results.R \
#     --tables DI=/path/res_table_DI.txt,DWS=/path/res_table_DWS.txt,DS=/path/res_table_DS.txt \
#     --out_dir results_filtered \
#     --p_cut 0.05 --q_cut 0.05 --lfc_cut 2 \
#     --circos DI=/path/DI_circos_input.txt,DWS=/path/DWS_circos_input.txt,DS=/path/DS_circos_input.txt \
#     --circos_id_col Gene
# =============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("optparse", quietly = TRUE)) {
    install.packages("optparse", repos = "https://cloud.r-project.org")
  }
  library(optparse)
})

# ----------------------------- helpers ---------------------------------
parse_kv_list <- function(x) {
  # "DI=a,DS=b" -> named vector c(DI="a", DS="b")
  parts <- unlist(strsplit(x, ","))
  parts <- trimws(parts)
  kv <- strsplit(parts, "=", fixed = TRUE)
  keys <- vapply(kv, function(z) z[[1]], character(1))
  vals <- vapply(kv, function(z) paste(z[-1], collapse="="), character(1))
  vals <- trimws(vals)
  stats::setNames(vals, keys)
}

detect_log2fc_col <- function(df) {
  candidates <- c("log2.fold_change.", "log2_fold_change", "log2_fold_change.", "log2.fold_change",
                  "log2FoldChange", "log2_fc", "log2FC")
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) > 0) return(hit[[1]])
  # last resort: try partial match
  pm <- grep("log2", colnames(df), ignore.case = TRUE, value = TRUE)
  if (length(pm) == 1) return(pm[[1]])
  stop("Could not find a log2 fold-change column. Found columns: ", paste(colnames(df), collapse=", "))
}

read_table_safe <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
  read.table(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = TRUE, quote = "", comment.char = "")
}

ensure_cols <- function(df, needed) {
  miss <- setdiff(needed, colnames(df))
  if (length(miss) > 0) stop("Missing required columns: ", paste(miss, collapse=", "))
}

# ----------------------------- CLI -------------------------------------
option_list <- list(
  make_option("--tables", type="character", default=NULL,
              help="Comma-separated KEY=path pairs, e.g., DI=...txt,DWS=...txt,DS=...txt"),
  make_option("--out_dir", type="character", default="filtered_outputs",
              help="Output directory [default: %default]"),
  make_option("--p_cut", type="double", default=0.05, help="p-value cutoff [default: %default]"),
  make_option("--q_cut", type="double", default=0.05, help="q-value cutoff [default: %default]"),
  make_option("--lfc_cut", type="double", default=2, help="|log2FC| cutoff [default: %default]"),
  make_option("--circos", type="character", default=NULL,
              help="Optional circos inputs as KEY=path pairs (same keys as --tables)"),
  make_option("--circos_id_col", type="character", default="Gene",
              help="Column name in circos input containing gene/ID [default: %default]")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$tables)) stop("Provide --tables KEY=path pairs.")
tables <- parse_kv_list(opt$tables)

out_dir <- normalizePath(opt$out_dir, winslash = "/", mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(file.path(out_dir, "selected_ids"), showWarnings = FALSE)
dir.create(file.path(out_dir, "selected_tables"), showWarnings = FALSE)

circos <- NULL
if (!is.null(opt$circos)) {
  circos <- parse_kv_list(opt$circos)
  dir.create(file.path(out_dir, "circos_subset"), showWarnings = FALSE)
}

# ----------------------------- processing --------------------------------
summary_rows <- list()

filter_one <- function(df, p_cut, q_cut, lfc_cut) {
  ensure_cols(df, c("gene_id", "p_value", "q_value"))
  lfc_col <- detect_log2fc_col(df)
  df$.__log2fc__ <- as.numeric(df[[lfc_col]])

  p_sig   <- df[df$p_value < p_cut, ]
  q_sig   <- df[df$q_value < q_cut, ]
  lfc_sig <- df[abs(df$.__log2fc__) >= lfc_cut, ]
  q_lfc   <- df[(df$q_value < q_cut) & (abs(df$.__log2fc__) >= lfc_cut), ]

  list(
    lfc_col = lfc_col,
    p_sig = p_sig,
    q_sig = q_sig,
    lfc_sig = lfc_sig,
    q_lfc = q_lfc
  )
}

for (grp in names(tables)) {
  df <- read_table_safe(tables[[grp]])
  res <- filter_one(df, opt$p_cut, opt$q_cut, opt$lfc_cut)

  # write selected q+lfc
  out_ids <- file.path(out_dir, "selected_ids", paste0(grp, "_sig_q_lfc_ids.txt"))
  write.table(as.character(res$q_lfc$gene_id), file = out_ids,
              row.names = FALSE, col.names = FALSE, quote = FALSE)

  out_tab <- file.path(out_dir, "selected_tables", paste0(grp, "_sig_q_lfc_table.tsv"))
  df_out <- res$q_lfc
  df_out$log2FC_used_column <- res$lfc_col
  write.table(df_out, file = out_tab, sep = "\t", row.names = FALSE, quote = FALSE)

  # circos subset (optional)
  if (!is.null(circos) && grp %in% names(circos)) {
    cdf <- read_table_safe(circos[[grp]])
    ensure_cols(cdf, c(opt$circos_id_col))
    keep <- cdf[cdf[[opt$circos_id_col]] %in% res$q_lfc$gene_id, ]
    out_circos <- file.path(out_dir, "circos_subset", paste0(grp, "_circos_subset.tsv"))
    write.table(keep, file = out_circos, sep = "\t", row.names = FALSE, quote = FALSE)
  }

  summary_rows[[grp]] <- data.frame(
    group = grp,
    n_total = nrow(df),
    n_p_sig = nrow(res$p_sig),
    n_q_sig = nrow(res$q_sig),
    n_lfc_sig = nrow(res$lfc_sig),
    n_q_lfc = nrow(res$q_lfc),
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summary_rows)
write.table(summary_df, file = file.path(out_dir, "summary_counts.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

message("Done. Outputs: ", out_dir)
