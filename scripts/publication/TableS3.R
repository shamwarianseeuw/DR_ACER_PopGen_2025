###############################################################################
# Supplementary Table S3: pairwise IBS distances within and between genets
###############################################################################

ibs <- read.delim(
  "processed_data/all_pairwise_IBS_distances.tsv",
  stringsAsFactors = FALSE
)

cat("Total pairwise comparisons:", nrow(ibs), "\n")
cat("Expected for 45 colonies:", choose(45, 2), "\n\n")

if (nrow(ibs) != choose(45, 2)) {
  stop("Expected 990 pairwise comparisons.")
}

comparison_order <- c(
  "within_AC1",
  "within_AC2",
  "within_AC3",
  "AC1_vs_AC2",
  "AC1_vs_AC3",
  "AC2_vs_AC3"
)

ibs$comparison <- factor(
  ibs$comparison,
  levels = comparison_order
)

summarize_group <- function(x) {
  c(
    n_pairs = length(x),
    mean = mean(x),
    sd = sd(x),
    median = median(x),
    min = min(x),
    max = max(x),
    Q25 = unname(quantile(x, 0.25)),
    Q75 = unname(quantile(x, 0.75))
  )
}

# Summary by comparison
results <- do.call(
  rbind,
  lapply(
    comparison_order,
    function(g) {
      summarize_group(
        ibs$IBS_distance[ibs$comparison == g]
      )
    }
  )
)

results <- data.frame(
  comparison = c(
    "Within AC1",
    "Within AC2",
    "Within AC3",
    "AC1-AC2",
    "AC1-AC3",
    "AC2-AC3"
  ),
  results,
  row.names = NULL
)

# Overall within- vs between-genet ranges
within <- ibs$IBS_distance[
  ibs$genotype1 == ibs$genotype2
]

between <- ibs$IBS_distance[
  ibs$genotype1 != ibs$genotype2
]

cat("WITHIN-GENET COMPARISONS\n")
cat("n =", length(within), "\n")
cat("min =", min(within), "\n")
cat("max =", max(within), "\n\n")

cat("BETWEEN-GENET COMPARISONS\n")
cat("n =", length(between), "\n")
cat("min =", min(between), "\n")
cat("max =", max(between), "\n\n")

cat(
  "Gap between maximum within-genet and minimum between-genet distance =",
  min(between) - max(within),
  "\n\n"
)

cat(
  "Do distributions overlap?",
  ifelse(max(within) >= min(between), "YES", "NO"),
  "\n\n"
)

print(results, digits = 8)

table_S3 <- results

for (j in c(
  "mean", "sd", "median",
  "min", "max", "Q25", "Q75"
)) {
  table_S3[[j]] <- round(table_S3[[j]], 4)
}

table_S3$n_pairs <- as.integer(table_S3$n_pairs)

cat("\nTABLE S3\n\n")
print(table_S3, row.names = FALSE)

if (interactive()) {
  View(table_S3)
}

# Save table
dir.create("supplementary", showWarnings = FALSE)

write.table(
  table_S3,
  "supplementary/Supplementary_Table_S3.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
