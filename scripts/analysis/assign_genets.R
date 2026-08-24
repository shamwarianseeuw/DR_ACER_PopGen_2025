# Assign multilocus genets from the IBS distance matrix

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop(
    "Usage: Rscript assign_genets.R ",
    "<ibsMat_file> <sample_order_file> <output_file>"
  )
}

ibs_file <- args[1]
sample_file <- args[2]
output_file <- args[3]

# IBS matrix and sample order
ibs <- as.matrix(
  read.table(
    ibs_file,
    header = FALSE,
    check.names = FALSE
  )
)

samples <- readLines(sample_file)
samples <- basename(samples)
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

# Average-linkage clustering
hc <- hclust(
  as.dist(ibs),
  method = "average"
)

clusters <- cutree(
  hc,
  k = 3
)

cluster_sizes <- table(clusters)

if (!identical(
  sort(as.integer(cluster_sizes)),
  c(8L, 14L, 23L)
)) {
  stop(
    "Unexpected cluster sizes: ",
    paste(sort(as.integer(cluster_sizes)), collapse = ", ")
  )
}

# Match cluster sizes to genet labels
size_to_genet <- c(
  "8" = "AC1",
  "14" = "AC3",
  "23" = "AC2"
)

cluster_to_genet <- setNames(
  size_to_genet[
    as.character(cluster_sizes)
  ],
  names(cluster_sizes)
)

genet <- unname(
  cluster_to_genet[
    as.character(clusters)
  ]
)

assignments <- data.frame(
  sample = names(clusters),
  genotype = genet,
  stringsAsFactors = FALSE
)

assignments <- assignments[
  match(samples, assignments$sample),
]

if (anyNA(assignments$genotype)) {
  stop("Some clusters could not be assigned to a genet.")
}

cat("\nGenet counts:\n")
print(table(assignments$genotype))

write.table(
  assignments,
  output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
