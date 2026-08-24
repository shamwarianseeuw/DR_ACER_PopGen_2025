# Genotypic richness and Simpson/Fager evenness

nursery <- read.delim(
  "metadata/sample_nursery.tsv",
  stringsAsFactors = FALSE
)

genets <- read.delim(
  "metadata/genotype_assignments.tsv",
  stringsAsFactors = FALSE
)

# Standardize sample IDs
if ("Tube_ID" %in% names(nursery)) {
  names(nursery)[names(nursery) == "Tube_ID"] <- "sample"
}

if ("genotype" %in% names(genets)) {
  names(genets)[names(genets) == "genotype"] <- "genet"
}

genets$sample <- sub(
  "\\.trim\\.bt2$",
  "",
  genets$sample
)

# Join nursery and genet assignments
dat <- merge(
  nursery[, c("sample", "Nursery")],
  genets[, c("sample", "genet")],
  by = "sample"
)

if (nrow(dat) != 45) {
  stop(
    "Expected 45 matched samples but found ",
    nrow(dat)
  )
}

if (anyNA(dat$Nursery) || anyNA(dat$genet)) {
  stop("Missing nursery or genet assignment after join.")
}

calculate_R <- function(counts) {
  counts <- counts[counts > 0]

  N <- sum(counts)
  G <- length(counts)

  (G - 1) / (N - 1)
}

calculate_V <- function(counts) {
  counts <- counts[counts > 0]

  N <- sum(counts)
  G <- length(counts)

  D <- 1 -
    sum(counts * (counts - 1)) /
    (N * (N - 1))

  Dmin <- ((G - 1) * (2 * N - G)) /
    (N * (N - 1))

  Dmax <- (N * (G - 1)) /
    (G * (N - 1))

  (D - Dmin) / (Dmax - Dmin)
}

summarize_group <- function(x, label) {
  counts <- table(
    factor(
      x$genet,
      levels = c("AC1", "AC2", "AC3")
    )
  )

  data.frame(
    Nursery = label,
    AC1 = unname(counts["AC1"]),
    AC2 = unname(counts["AC2"]),
    AC3 = unname(counts["AC3"]),
    N = sum(counts),
    G = sum(counts > 0),
    R = calculate_R(counts),
    V = calculate_V(counts),
    stringsAsFactors = FALSE
  )
}

# Calculate metrics by nursery and overall
cap_cana <- summarize_group(
  dat[dat$Nursery == "Cap Cana", ],
  "Cap Cana"
)

acuario <- summarize_group(
  dat[dat$Nursery == "Acuario", ],
  "Acuario"
)

overall <- summarize_group(
  dat,
  "Overall"
)

results <- rbind(
  cap_cana,
  acuario,
  overall
)

# Check sample totals
stopifnot(
  cap_cana$N == 19,
  acuario$N == 26,
  overall$N == 45
)

cat("\nGENOTYPIC DIVERSITY RESULTS\n\n")
print(results, digits = 8)

# Save results
write.table(
  results,
  "processed_data/genotypic_diversity_metrics.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
