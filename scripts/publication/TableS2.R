##################################################################################
# Supplementary Table S2: sample metadata, genet assignment, and sequencing depth
##################################################################################

meta <- read.delim(
  "metadata/sample_nursery.tsv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

geno <- read.delim(
  "metadata/genotype_assignments.tsv",
  stringsAsFactors = FALSE
)

depth <- read.delim(
  "processed_data/bam_depth_all_2515_sites.tsv",
  header = FALSE,
  stringsAsFactors = FALSE
)

colnames(depth) <- c(
  "sample",
  "SNP_sites_with_coverage",
  "mean_depth_all_2515_sites"
)

# Standardize IDs
if ("Tube_ID" %in% names(meta)) {
  names(meta)[names(meta) == "Tube_ID"] <- "sample"
}

geno$sample <- sub("\\.trim\\.bt2$", "", geno$sample)
depth$sample <- sub("\\.trim\\.bt2$", "", depth$sample)

# Join
tab <- merge(
  meta[, c("sample", "Nursery")],
  geno[, c("sample", "genotype")],
  by = "sample"
)

tab <- merge(
  tab,
  depth[, c(
    "sample",
    "SNP_sites_with_coverage",
    "mean_depth_all_2515_sites"
  )],
  by = "sample"
)

if (nrow(tab) != 45) {
  stop("Expected 45 joined samples; found ", nrow(tab))
}

if (anyNA(tab)) {
  stop("Missing values detected after joining Table S2 inputs.")
}

# Sort N99-N143 numerically
tab$sample_number <- as.numeric(sub("^N", "", tab$sample))
tab <- tab[order(tab$sample_number), ]

#Depth summary
cat("Number of joined samples:", nrow(tab), "\n\n")

cat("Nursery x genet table:\n")
print(table(tab$Nursery, tab$genotype))
cat("\n")

x <- tab$mean_depth_all_2515_sites

cat("DEPTH ACROSS ALL 2,515 RETAINED SNP SITES\n")
cat("n =", length(x), "\n")
cat("mean =", mean(x), "\n")
cat("SD =", sd(x), "\n")
cat("median =", median(x), "\n")
cat("min =", min(x), "\n")
cat("max =", max(x), "\n\n")

table_S2 <- data.frame(
  Tube_ID = tab$sample,
  Nursery = tab$Nursery,
  Genet = tab$genotype,
  SNP_sites_with_coverage = tab$SNP_sites_with_coverage,
  Mean_depth_all_2515_sites = round(tab$mean_depth_all_2515_sites, 2)
)

print(table_S2, row.names = FALSE)

if (interactive()) {
  View(table_S2)
}

#Save table
dir.create("supplementary", showWarnings = FALSE)

write.table(
  table_S2,
  "supplementary/Supplementary_Table_S2.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
