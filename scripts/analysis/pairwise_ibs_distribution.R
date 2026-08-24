# Pairwise IBS distances within and between genets

ibs_file <- "processed_data/ref_main.ibsMat"
sample_file <- "processed_data/ref_main_sample_order.txt"
genet_file <- "metadata/genotype_assignments.tsv"

# IBS matrix and sample order
ibs <- as.matrix(
  read.table(
    ibs_file,
    header = FALSE,
    check.names = FALSE
  )
)

samples <- readLines(sample_file)
samples <- sub("\\.bam$", "", samples)

if (
  nrow(ibs) != length(samples) ||
  ncol(ibs) != length(samples)
) {
  stop("IBS matrix dimensions do not match the sample list.")
}

if (max(abs(ibs - t(ibs))) > 1e-10) {
  stop("IBS matrix is not symmetric.")
}

diag(ibs) <- 0
dimnames(ibs) <- list(samples, samples)

# Genet assignments
genets <- read.delim(
  genet_file,
  stringsAsFactors = FALSE
)

if (!all(samples %in% genets$sample)) {
  stop("Some samples are missing genet assignments.")
}

genet_lookup <- setNames(
  genets$genotype,
  genets$sample
)

# Unique pairwise comparisons
idx <- which(
  upper.tri(ibs),
  arr.ind = TRUE
)

pairs <- data.frame(
  sample1 = rownames(ibs)[idx[, 1]],
  sample2 = colnames(ibs)[idx[, 2]],
  IBS_distance = ibs[idx],
  stringsAsFactors = FALSE
)

pairs$genotype1 <- unname(
  genet_lookup[pairs$sample1]
)

pairs$genotype2 <- unname(
  genet_lookup[pairs$sample2]
)

if (
  anyNA(pairs$genotype1) ||
  anyNA(pairs$genotype2)
) {
  stop("Missing genet assignments in pairwise comparisons.")
}

# Comparison groups
pairs$comparison <- mapply(
  function(g1, g2) {
    if (g1 == g2) {
      paste0("within_", g1)
    } else {
      paste(
        sort(c(g1, g2)),
        collapse = "_vs_"
      )
    }
  },
  pairs$genotype1,
  pairs$genotype2
)

comparison_order <- c(
  "within_AC1",
  "within_AC2",
  "within_AC3",
  "AC1_vs_AC2",
  "AC1_vs_AC3",
  "AC2_vs_AC3"
)

pairs$comparison <- factor(
  pairs$comparison,
  levels = comparison_order
)

if (anyNA(pairs$comparison)) {
  stop("Unexpected genet comparison category.")
}

if (nrow(pairs) != choose(45, 2)) {
  stop(
    "Expected 990 pairwise comparisons but found ",
    nrow(pairs)
  )
}

cat("\nTotal pairwise comparisons:", nrow(pairs), "\n\n")

cat("Pairs by comparison:\n")
print(table(pairs$comparison))

# Overall within- and between-genet ranges
within <- pairs$IBS_distance[
  pairs$genotype1 == pairs$genotype2
]

between <- pairs$IBS_distance[
  pairs$genotype1 != pairs$genotype2
]

cat("\nWithin-genet IBS range:\n")
cat(min(within), "to", max(within), "\n")

cat("\nBetween-genet IBS range:\n")
cat(min(between), "to", max(between), "\n")

cat(
  "\nDistributions overlap:",
  ifelse(
    max(within) >= min(between),
    "YES",
    "NO"
  ),
  "\n"
)

# Save pairwise dataset
write.table(
  pairs,
  "processed_data/all_pairwise_IBS_distances.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
