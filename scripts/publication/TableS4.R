###############################################################################
# Supplementary Table S4: relatedness within and among multilocus genets
###############################################################################

# Full 45-colony dataset
full <- read.delim(
  "processed_data/full_relatedness_labeled.tsv",
  stringsAsFactors = FALSE
)

# Check number of pairwise comparisons
expected_pairs <- choose(45, 2)

if (nrow(full) != expected_pairs) {
  stop(
    "Expected ", expected_pairs,
    " full-dataset pairs but found ",
    nrow(full)
  )
}

comparison_order <- c(
  "Within AC1",
  "Within AC2",
  "Within AC3",
  "AC1-AC2",
  "AC1-AC3",
  "AC2-AC3"
)

full$comparison <- factor(
  full$comparison,
  levels = comparison_order
)

# Relatedness by comparison
summaries <- lapply(
  comparison_order,
  function(g) {
    
    z <- full$rab[
      full$comparison == g
    ]
    
    data.frame(
      comparison = g,
      n_pairs = length(z),
      mean = mean(z),
      sd = sd(z),
      min = min(z),
      max = max(z)
    )
  }
)

table_S4A <- do.call(
  rbind,
  summaries
)

rownames(table_S4A) <- NULL

# One representative per genet
reduced <- read.delim(
  "processed_data/ngsrelate_reduced_3genets.txt",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Sample order in the reduced BCF:
#   0 = N120 = AC3
#   1 = N128 = AC2
#   2 = N130 = AC1

representative_lookup <- data.frame(
  index = c(0, 1, 2),
  sample = c(
    "N120",
    "N128",
    "N130"
  ),
  genet = c(
    "AC3",
    "AC2",
    "AC1"
  ),
  stringsAsFactors = FALSE
)

# Check number of comparisons
if (nrow(reduced) != choose(3, 2)) {
  stop(
    "Expected 3 reduced-dataset comparisons but found ",
    nrow(reduced)
  )
}

get_sample <- function(i) {
  representative_lookup$sample[
    match(i, representative_lookup$index)
  ]
}

get_genet <- function(i) {
  representative_lookup$genet[
    match(i, representative_lookup$index)
  ]
}

s1 <- get_sample(reduced$a)
s2 <- get_sample(reduced$b)

g1 <- get_genet(reduced$a)
g2 <- get_genet(reduced$b)

# Keep genet labels in a consistent order
comparison <- vapply(
  seq_len(nrow(reduced)),
  function(i) {
    paste(
      sort(c(g1[i], g2[i])),
      collapse = "-"
    )
  },
  character(1)
)

# Order representatives by sample ID
sample_number <- function(x) {
  as.numeric(sub("^N", "", x))
}

rep1 <- character(nrow(reduced))
rep2 <- character(nrow(reduced))

for (i in seq_len(nrow(reduced))) {
  
  if (sample_number(s1[i]) <= sample_number(s2[i])) {
    rep1[i] <- s1[i]
    rep2[i] <- s2[i]
  } else {
    rep1[i] <- s2[i]
    rep2[i] <- s1[i]
  }
}

table_S4B <- data.frame(
  comparison = comparison,
  representative_1 = rep1,
  representative_2 = rep2,
  nSites = reduced$nSites,
  rab = reduced$rab,
  stringsAsFactors = FALSE
)

table_S4B <- table_S4B[
  match(
    c(
      "AC1-AC2",
      "AC1-AC3",
      "AC2-AC3"
    ),
    table_S4B$comparison
  ),
]

rownames(table_S4B) <- NULL

# Check outputs
cat("\nTABLE S4A — FULL DATASET\n\n")
print(
  table_S4A,
  digits = 10
)

cat(
  "\nTotal pairwise comparisons:",
  sum(table_S4A$n_pairs),
  "\n"
)

cat("\nTABLE S4B — REDUCED DATASET\n\n")
print(
  table_S4B,
  row.names = FALSE
)

# Final checks
stopifnot(
  sum(table_S4A$n_pairs) == 990,
  nrow(table_S4B) == 3,
  table_S4B$comparison ==
    c("AC1-AC2", "AC1-AC3", "AC2-AC3")
)

View(table_S4A)
View(table_S4B)

# Save tables
dir.create("supplementary", showWarnings = FALSE)

write.table(
  table_S4A,
  "supplementary/Supplementary_Table_S4A_full_relatedness.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  table_S4B,
  "supplementary/Supplementary_Table_S4B_reduced_relatedness.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)