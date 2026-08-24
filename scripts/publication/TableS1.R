###############################################################################
# Supplementary Table S1: SNP filtering summary
###############################################################################

infile <- "processed_data/SNP_filter_counts.txt"

x <- readLines(infile)

extract_number <- function(pattern) {
  z <- grep(pattern, x, value = TRUE)

  if (length(z) != 1) {
    stop(
      "Expected exactly one match for: ",
      pattern
    )
  }

  as.numeric(
    sub(".*: *([0-9.]+)%?.*", "\\1", z)
  )
}

candidate <- extract_number(
  "^Candidate SNPs before final SNP-level QC:"
)

minind <- extract_number(
  "^Pass minInd >= 34:"
)

maf <- extract_number(
  "^Pass MAF >= 0.05:"
)

both <- extract_number(
  "^Pass minInd >= 34 AND MAF >= 0.05:"
)

final <- extract_number(
  "^Final QC-filtered SNPs:"
)

retention <- extract_number(
  "^Final retention:"
)

table_S1 <- data.frame(
  Filter_or_stage = c(
    "Candidate SNPs before final SNP-level QC",
    "Pass minInd >= 34",
    "Pass MAF >= 0.05",
    "Pass minInd >= 34 and MAF >= 0.05",
    "Final QC-filtered SNPs"
  ),
  SNPs_retained = c(
    candidate,
    minind,
    maf,
    both,
    final
  ),
  Percent_of_candidate_SNPs = round(
    100 * c(
      candidate,
      minind,
      maf,
      both,
      final
    ) / candidate,
    1
  ),
  stringsAsFactors = FALSE
)

# Check final retention
calculated_retention <- round(
  100 * final / candidate,
  1
)

if (calculated_retention != retention) {
  stop(
    "Recorded final retention (",
    retention,
    "%) does not match calculated retention (",
    calculated_retention,
    "%)."
  )
}

cat("\nTABLE S1\n\n")
print(table_S1, row.names = FALSE)

View(table_S1)

# Save table
dir.create("supplementary", showWarnings = FALSE)

write.table(
  table_S1,
  "supplementary/Supplementary_Table_S1.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)