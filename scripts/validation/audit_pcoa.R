# Check PCoA results from the ANGSD IBS matrix

# IBS matrix and sample order
D <- as.matrix(
  read.table(
    "processed_data/ref_main.ibsMat",
    header = FALSE
  )
)

ids <- readLines(
  "processed_data/ref_main_sample_order.txt"
)

ids <- sub("\\.bam$", "", ids)

cat("Matrix dimensions:", dim(D), "\n")
cat("Number of samples:", length(ids), "\n")
cat(
  "Symmetric:",
  isTRUE(all.equal(D, t(D), tolerance = 1e-8)),
  "\n"
)
cat("Diagonal range:", range(diag(D)), "\n\n")

if (
  nrow(D) != length(ids) ||
  ncol(D) != length(ids)
) {
  stop("IBS matrix dimensions do not match the sample list.")
}

rownames(D) <- ids
colnames(D) <- ids

# PCoA
pcoa <- cmdscale(
  as.dist(D),
  k = 10,
  eig = TRUE,
  add = FALSE
)

eig <- pcoa$eig

cat("First 10 eigenvalues:\n")
print(eig[1:10], digits = 10)

cat("\nNumber positive:", sum(eig > 0), "\n")
cat("Number negative:", sum(eig < 0), "\n")
cat("Most negative eigenvalue:", min(eig), "\n\n")

# Variation explained using positive eigenvalues
positive <- eig[eig > 0]
pct_positive <- 100 * positive / sum(positive)

cat("VARIATION EXPLAINED - POSITIVE EIGENVALUES\n")
cat("PCoA1 =", pct_positive[1], "%\n")
cat("PCoA2 =", pct_positive[2], "%\n")
cat(
  "PCoA1 + PCoA2 =",
  sum(pct_positive[1:2]),
  "%\n\n"
)

# Compare with percentages based on all eigenvalues
pct_all <- 100 * eig / sum(eig)

cat("RELATIVE TO SUM OF ALL EIGENVALUES\n")
cat("PCoA1 =", pct_all[1], "%\n")
cat("PCoA2 =", pct_all[2], "%\n\n")

# Compare reconstructed coordinates with saved coordinates
old <- read.delim(
  "processed_data/ref_main_ibs_pca_coordinates.tsv",
  stringsAsFactors = FALSE
)

old <- old[match(ids, old$sample), ]

if (anyNA(old$sample)) {
  stop("Some samples are missing from the saved PCoA coordinates.")
}

cat(
  "Samples correctly matched:",
  all(old$sample == ids),
  "\n\n"
)

# Axis signs can reverse without changing the ordination
cor1 <- cor(
  pcoa$points[, 1],
  old$PC1
)

cor2 <- cor(
  pcoa$points[, 2],
  old$PC2
)

cat("Coordinate correlations:\n")
cat("Axis 1 correlation =", cor1, "\n")
cat("Axis 2 correlation =", cor2, "\n")
cat("Absolute Axis 1 correlation =", abs(cor1), "\n")
cat("Absolute Axis 2 correlation =", abs(cor2), "\n")
