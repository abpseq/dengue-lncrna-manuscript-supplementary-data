#!/usr/bin/env Rscript
# =============================================================================
# 02_cummeRbund_heatmaps.R
# -----------------------------------------------------------------------------
# Purpose
#   Generate cummeRbund csHeatmap() plots for gene/lncRNA ID sets per group.
#
# Inputs
#   - Cuffdiff output folders (cummeRbund readable): e.g., cuffdiff_DI, cuffdiff_DWS, cuffdiff_DS
#   - Gene/lncRNA IDs list per group:
#       * either a simple text file (one ID per line), OR
#       * a result table containing a 'gene_id' column
#
# Outputs
#   - out_dir/<group>_csHeatmap.png  (600 dpi by default)
#
# Example
#   Rscript code/02_cummeRbund_heatmaps.R \
#     --cuffdirs DI=/path/cuffdiff_DI,DWS=/path/cuffdiff_DWS,DS=/path/cuffdiff_DS \
#     --gene_inputs DI=/path/res_table_DI.txt,DWS=/path/res_table_DWS.txt,DS=/path/res_table_DS.txt \
#     --gene_input_type table \
#     --out_dir heatmaps --dpi 600
# =============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("optparse", quietly = TRUE)) {
    install.packages("optparse", repos = "https://cloud.r-project.org")
  }
  library(optparse)

  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
  if (!requireNamespace("cummeRbund", quietly = TRUE)) BiocManager::install("cummeRbund", ask = FALSE, update = FALSE)
  library(cummeRbund)
})

parse_kv_list <- function(x) {
  parts <- unlist(strsplit(x, ","))
  parts <- trimws(parts)
  kv <- strsplit(parts, "=", fixed = TRUE)
  keys <- vapply(kv, function(z) z[[1]], character(1))
  vals <- vapply(kv, function(z) paste(z[-1], collapse="="), character(1))
  vals <- trimws(vals)
  stats::setNames(vals, keys)
}

read_gene_ids <- function(path, input_type=c("list","table")) {
  input_type <- match.arg(input_type)
  if (!file.exists(path)) stop("Missing file: ", path)
  if (input_type == "list") {
    ids <- scan(path, what = "character", quiet = TRUE)
    return(unique(ids))
  }
  df <- read.table(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = TRUE, quote = "", comment.char = "")
  if (!("gene_id" %in% colnames(df))) stop("Table missing gene_id column: ", path)
  unique(as.character(df$gene_id))
}

option_list <- list(
  make_option("--cuffdirs", type="character", default=NULL,
              help="KEY=dir pairs, e.g., DI=/path/cuffdiff_DI,DWS=/path/cuffdiff_DWS,DS=/path/cuffdiff_DS"),
  make_option("--gene_inputs", type="character", default=NULL,
              help="KEY=path pairs for ID lists or tables"),
  make_option("--gene_input_type", type="character", default="table",
              help="Either 'table' (expects gene_id column) or 'list' (1 ID per line) [default: %default]"),
  make_option("--out_dir", type="character", default="heatmaps",
              help="Output directory [default: %default]"),
  make_option("--dpi", type="integer", default=600, help="PNG dpi [default: %default]"),
  make_option("--width_in", type="double", default=5, help="Figure width in inches [default: %default]"),
  make_option("--height_in", type="double", default=8, help="Figure height in inches [default: %default]"),
  make_option("--cluster", type="character", default="both", help="csHeatmap cluster arg [default: %default]")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$cuffdirs) || is.null(opt$gene_inputs)) stop("Provide --cuffdirs and --gene_inputs")

cuffdirs <- parse_kv_list(opt$cuffdirs)
gene_inputs <- parse_kv_list(opt$gene_inputs)

if (!all(names(cuffdirs) %in% names(gene_inputs))) {
  stop("Keys in --cuffdirs must also exist in --gene_inputs. cuffdirs keys: ",
       paste(names(cuffdirs), collapse=","), " gene_input keys: ", paste(names(gene_inputs), collapse=","))
}

out_dir <- normalizePath(opt$out_dir, winslash = "/", mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (grp in names(cuffdirs)) {
  cuff_dir <- normalizePath(cuffdirs[[grp]], winslash = "/", mustWork = FALSE)
  ids_path <- normalizePath(gene_inputs[[grp]], winslash = "/", mustWork = FALSE)

  ids <- read_gene_ids(ids_path, opt$gene_input_type)
  message("Group ", grp, ": loaded ", length(ids), " IDs")

  cuff <- readCufflinks(cuff_dir)
  myGenes <- getGenes(cuff, ids)

  out_png <- file.path(out_dir, paste0(grp, "_csHeatmap.png"))
  png(out_png, units = "in", width = opt$width_in, height = opt$height_in, res = opt$dpi)
  csHeatmap(myGenes, cluster = opt$cluster)
  dev.off()

  message("Saved: ", out_png)
}

message("Done.")
