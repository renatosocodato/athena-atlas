#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(limma))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: s06_edgeR_ql_contrasts.R <counts.csv> <metadata.csv> <route> <output.csv>")
}

counts_path <- args[[1]]
metadata_path <- args[[2]]
route_name <- args[[3]]
output_path <- args[[4]]

raw_counts <- read.csv(counts_path, check.names = FALSE, colClasses = "character", na.strings = NULL)
genes <- raw_counts[[1L]]
counts <- as.matrix(raw_counts[-1L])
storage.mode(counts) <- "double"
rownames(counts) <- genes
meta <- read.csv(metadata_path, check.names = FALSE, colClasses = "character", na.strings = NULL)

if (anyDuplicated(genes) || anyDuplicated(colnames(counts)) || any(!is.finite(counts)) ||
    any(counts < 0) || any(counts != floor(counts)) || any(colSums(counts) <= 0)) {
  stop("Unique genes/samples, finite nonnegative integer counts and positive libraries required")
}

if (!"sample_id" %in% colnames(meta)) {
  stop("Metadata must include sample_id")
}
if (!"sex_route" %in% colnames(meta)) {
  stop("Metadata must include sex_route")
}
if (!"diagnosis" %in% colnames(meta)) {
  stop("Metadata must include diagnosis")
}
if (!"donor_id" %in% colnames(meta) || anyDuplicated(meta$sample_id) ||
    anyDuplicated(meta$donor_id) || !setequal(meta$sample_id, colnames(counts))) {
  stop("Counts must match one sample per independent donor")
}

meta <- meta[meta$sex_route == route_name, , drop = FALSE]
if (nrow(meta) < 4) {
  stop(sprintf("Route %s does not contain enough samples for edgeR", route_name))
}

keep_cols <- intersect(meta$sample_id, colnames(counts))
if (length(keep_cols) < 4) {
  stop(sprintf("Route %s does not retain enough count columns", route_name))
}
meta <- meta[match(keep_cols, meta$sample_id), , drop = FALSE]
counts <- counts[, keep_cols, drop = FALSE]

group <- factor(meta$diagnosis, levels = c("control", "ad"))
if (anyNA(group) || any(table(group) < 2L)) {
  stop(sprintf("Route %s requires at least two independent donors per diagnosis", route_name))
}

dge <- DGEList(counts = counts, samples = meta, group = group)
keep <- filterByExpr(dge, group = group)
if (sum(keep) < 1) {
  stop(sprintf("Route %s retained no genes after filterByExpr", route_name))
}
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

if ("braak_numeric" %in% colnames(meta)) {
  braak_numeric <- suppressWarnings(as.numeric(meta$braak_numeric))
  if (any(!is.finite(braak_numeric))) stop("Braak covariate must be finite")
}
if ("braak_numeric" %in% colnames(meta) && length(unique(braak_numeric)) > 1) {
  design <- model.matrix(~ group + braak_numeric)
} else {
  design <- model.matrix(~ group)
}
if (qr(design)$rank < ncol(design) || nrow(design) <= ncol(design)) {
  stop("Diagnosis contrast is not estimable with this design")
}

dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)
coef_name <- grep("^groupad$", colnames(design), value = TRUE)
if (length(coef_name) != 1) {
  stop(sprintf("Route %s could not resolve the groupad coefficient", route_name))
}
qlf <- glmQLFTest(fit, coef = which(colnames(design) == coef_name))
tab <- topTags(qlf, n = Inf, sort.by = "none")$table
tab$gene <- rownames(tab)
tab$sex_route <- route_name
tab$n_ad_donors <- sum(meta$diagnosis == "ad")
tab$n_control_donors <- sum(meta$diagnosis == "control")
write.csv(tab, output_path, row.names = FALSE)
