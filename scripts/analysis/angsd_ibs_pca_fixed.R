args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop(
    "Usage: Rscript angsd_ibs_pca_fixed.R ",
    "<ibsMat_file> [bam_list] [inds2pops] [out_prefix]"
  )
}

ibs_file <- args[1]
bam_file <- if (length(args) >= 2) args[2] else "bams_ref"
inds_file <- if (length(args) >= 3) args[3] else ""

out_prefix <- if (length(args) >= 4) {
  args[4]
} else {
  sub("\\.ibsMat$", "", basename(ibs_file))
}

suppressPackageStartupMessages({
  library(vegan)
})

# Sample order
bams <- read.table(
  bam_file,
  stringsAsFactors = FALSE
)[, 1]

samples <- basename(bams)
samples <- sub("\\.bam$", "", samples)

# IBS matrix
ma <- as.matrix(
  read.table(
    ibs_file,
    check.names = FALSE
  )
)

if (
  nrow(ma) != length(samples) ||
  ncol(ma) != length(samples)
) {
  stop("IBS matrix dimensions do not match the sample list.")
}

dimnames(ma) <- list(samples, samples)

if (max(abs(ma - t(ma))) > 1e-10) {
  stop("IBS matrix is not symmetric.")
}

diag(ma) <- 0

# Population assignments, if provided
has_pops <- FALSE

if (
  nzchar(inds_file) &&
  file.exists(inds_file) &&
  file.info(inds_file)$size > 0
) {
  i2p <- read.table(
    inds_file,
    sep = "\t",
    stringsAsFactors = FALSE,
    header = FALSE
  )

  if (ncol(i2p) < 2) {
    stop("inds2pops must have at least 2 columns: sample<TAB>population")
  }

  rownames(i2p) <- i2p[, 1]
  pop <- i2p[samples, 2]
  has_pops <- TRUE
} else {
  pop <- factor(
    rep("all_samples", length(samples))
  )
}

pop <- as.factor(pop)
colors <- as.numeric(pop)

# Hierarchical clustering
hc <- hclust(
  as.dist(ma),
  method = "average"
)

# PCoA
pcoa <- cmdscale(
  as.dist(ma),
  k = 2,
  eig = TRUE
)

pc_percent <- round(
  100 * pcoa$eig /
    sum(pcoa$eig[pcoa$eig > 0]),
  2
)

pts <- as.data.frame(pcoa$points)
colnames(pts) <- c("PC1", "PC2")
pts$sample <- samples
pts$pop <- pop

pdf(
  paste0(out_prefix, "_ibs_pca.pdf"),
  width = 10,
  height = 8
)

par(mfrow = c(2, 2))

plot(
  hc,
  cex = 0.6,
  main = "Hierarchical clustering from IBS matrix",
  xlab = "",
  sub = ""
)

plot(
  pcoa$eig,
  type = "b",
  pch = 16,
  main = "Eigenvalues",
  xlab = "Axis",
  ylab = "Eigenvalue"
)

plot(
  pts$PC1,
  pts$PC2,
  pch = 19,
  col = colors,
  xlab = paste0(
    "PC1 (",
    ifelse(length(pc_percent) >= 1, pc_percent[1], NA),
    "%)"
  ),
  ylab = paste0(
    "PC2 (",
    ifelse(length(pc_percent) >= 2, pc_percent[2], NA),
    "%)"
  ),
  main = "PCoA from IBS matrix"
)

text(
  pts$PC1,
  pts$PC2,
  labels = pts$sample,
  pos = 3,
  cex = 0.55
)

if (has_pops) {
  ord <- capscale(
    ma ~ pop,
    data = data.frame(pop = pop)
  )

  plot(
    ord,
    choices = c(1, 2),
    type = "n",
    main = "PCoA/CAP by population"
  )

  points(
    ord,
    choices = c(1, 2),
    pch = 19,
    col = colors
  )

  ordispider(
    ord,
    choices = c(1, 2),
    groups = pop,
    col = "grey80"
  )

  ordiellipse(
    ord,
    choices = c(1, 2),
    groups = pop,
    draw = "polygon",
    col = seq_along(levels(pop)),
    alpha = 60,
    label = TRUE
  )
} else {
  plot.new()

  text(
    0.5,
    0.5,
    "No population file provided",
    cex = 1.2
  )
}

dev.off()

write.table(
  pts,
  file = paste0(
    out_prefix,
    "_ibs_pca_coordinates.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Wrote:\n")
cat("  ", paste0(out_prefix, "_ibs_pca.pdf"), "\n")
cat(
  "  ",
  paste0(out_prefix, "_ibs_pca_coordinates.tsv"),
  "\n"
)
