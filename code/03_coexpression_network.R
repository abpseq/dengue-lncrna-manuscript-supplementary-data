#!/usr/bin/env Rscript
# =============================================================================
# 03_coexpression_network.R
# -----------------------------------------------------------------------------
# Purpose
#   Build a co-expression edge list between:
#     - a set of mRNAs (gene symbols) and
#     - a set of novel lncRNA IDs (e.g., XLOC_* or transcript IDs)
#   using repFpkmMatrix from cummeRbund + Pearson correlation thresholding.
#
# Outputs
#   - edge list (TSV): Start, End
#   - optional edge list filtered to XLOC nodes
#
# Example
#   Rscript code/03_coexpression_network.R \
#     --cuffdiff_dir /path/cuffdiff_DS \
#     --mrna_list /path/sig_DEG_list.txt \
#     --novel_ids /path/group_DS_sel_result.fasta \
#     --out_edges out/groupDS_edges.tsv \
#     --out_edges_xloc out/groupDS_edges_XLOC.tsv \
#     --cutoff 0.7
# =============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("optparse", quietly = TRUE)) {
    install.packages("optparse", repos = "https://cloud.r-project.org")
  }
  library(optparse)

  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
  if (!requireNamespace("cummeRbund", quietly = TRUE)) BiocManager::install("cummeRbund", ask = FALSE, update = FALSE)
  if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph", repos = "https://cloud.r-project.org")

  library(cummeRbund)
  library(igraph)
})

read_mrna_list <- function(path) {
  if (!file.exists(path)) stop("Missing mRNA list: ", path)
  unique(scan(path, what = "character", quiet = TRUE))
}

read_novel_ids <- function(path) {
  if (!file.exists(path)) stop("Missing novel IDs file: ", path)
  if (grepl("\\.fasta$|\\.fa$|\\.fna$", path, ignore.case = TRUE)) {
    # Parse FASTA headers: take first token after '>'
    x <- readLines(path, warn = FALSE)
    hdr <- x[grepl("^>", x)]
    ids <- sub("^>", "", hdr)
    ids <- sub("\\s.*$", "", ids)
    return(unique(ids))
  }
  # Otherwise: try table/csv
  df <- tryCatch(
    read.table(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE),
    error = function(e) read.csv(path, header = FALSE, stringsAsFactors = FALSE)
  )
  # use 2nd col if present (legacy), else first
  col_to_use <- if (ncol(df) >= 2) 2 else 1
  unique(as.character(df[[col_to_use]]))
}

dedup_by_iqr <- function(expr_mat, gene_short_names) {
  # expr_mat rows align to gene_id rows subset; gene_short_names aligns to rows
  b <- as.numeric(apply(expr_mat, 1, IQR, na.rm = TRUE))
  cx <- unique(gene_short_names)
  maxin <- integer(length(cx))
  for (i in seq_along(cx)) {
    idx <- which(cx[i] == gene_short_names)
    maxin[i] <- idx[which.max(b[idx])]
  }
  expr_mat[maxin, , drop = FALSE]
}

option_list <- list(
  make_option("--cuffdiff_dir", type="character", default=NULL, help="Cuffdiff output dir (cummeRbund readable)"),
  make_option("--mrna_list", type="character", default=NULL, help="mRNA gene symbols list (1 per line)"),
  make_option("--novel_ids", type="character", default=NULL, help="Novel IDs (FASTA or 1/2-column table)"),
  make_option("--cutoff", type="double", default=0.7, help="Pearson correlation cutoff [default: %default]"),
  make_option("--method", type="character", default="pearson", help="Correlation method [default: %default]"),
  make_option("--out_edges", type="character", default="edges.tsv", help="Output edge list [default: %default]"),
  make_option("--out_edges_xloc", type="character", default=NULL, help="Optional output filtered to XLOC edges")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$cuffdiff_dir) || is.null(opt$mrna_list) || is.null(opt$novel_ids)) {
  stop("Provide --cuffdiff_dir, --mrna_list, --novel_ids")
}

cuff <- readCufflinks(normalizePath(opt$cuffdiff_dir, winslash = "/", mustWork = FALSE))
mrna_syms <- read_mrna_list(opt$mrna_list)
novel_ids <- read_novel_ids(opt$novel_ids)

# Map gene symbols -> gene_id via annotation
anno <- annotation(object = genes(cuff))
anno <- anno[, c("gene_id", "gene_short_name")]
anno$gene_short_name <- gsub("\\,.*", "", anno$gene_short_name)

# subset gene_ids for mRNAs in list
f1 <- anno[anno$gene_short_name %in% mrna_syms, , drop = FALSE]
x_loc <- f1$gene_id

expr <- repFpkmMatrix(genes(cuff))

# mRNA expression by gene_id, then deduplicate by IQR across duplicated gene symbols
expr_mrna <- expr[x_loc, , drop = FALSE]
expr_mrna_dedup <- dedup_by_iqr(expr_mrna, f1$gene_short_name)

# novel expression by ID (direct rownames)
novel_present <- novel_ids[novel_ids %in% rownames(expr)]
if (length(novel_present) == 0) stop("None of the novel IDs were found in repFpkmMatrix rownames.")
expr_novel <- expr[novel_present, , drop = FALSE]

# combine
expr_all <- rbind(expr_mrna_dedup, expr_novel)

# correlation between rows (genes), based on sample expression columns
cor_mat <- cor(t(expr_all), method = opt$method, use = "pairwise.complete.obs")
diag(cor_mat) <- 0

# threshold to adjacency
adj <- abs(cor_mat)
adj[adj < opt$cutoff] <- 0
adj[adj >= opt$cutoff] <- 1
adj[lower.tri(adj)] <- 0

# build graph and edge list
g <- graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE, weighted = NULL)
edgelist <- as.data.frame(get.edgelist(g), stringsAsFactors = FALSE)
colnames(edgelist) <- c("Start", "End")

# write outputs
dir.create(dirname(opt$out_edges), recursive = TRUE, showWarnings = FALSE)
write.table(edgelist, file = opt$out_edges, sep = "\t", row.names = FALSE, quote = FALSE)

if (!is.null(opt$out_edges_xloc)) {
  xloc_idx <- grepl("XLOC", edgelist$Start, ignore.case = TRUE) | grepl("XLOC", edgelist$End, ignore.case = TRUE)
  ed2 <- edgelist[xloc_idx, , drop = FALSE]
  dir.create(dirname(opt$out_edges_xloc), recursive = TRUE, showWarnings = FALSE)
  write.table(ed2, file = opt$out_edges_xloc, sep = "\t", row.names = FALSE, quote = FALSE)
}

message("Done. Wrote edges: ", opt$out_edges)
