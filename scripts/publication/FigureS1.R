# Supplementary Figure S1: sensitivity of genet assignment to clustering cut height

library(ggplot2)

# Data
ibs <- read.delim(
  "processed_data/all_pairwise_IBS_distances.tsv",
  stringsAsFactors = FALSE
)

assign <- read.delim(
  "metadata/genotype_assignments.tsv",
  stringsAsFactors = FALSE
)

# Reconstruct IBS distance matrix
samples <- sort(unique(c(ibs$sample1, ibs$sample2)))
n <- length(samples)

if (nrow(ibs) != choose(n, 2)) {
  stop("Unexpected number of pairwise IBS comparisons.")
}

D <- matrix(
  0,
  nrow = n,
  ncol = n,
  dimnames = list(samples, samples)
)

for (i in seq_len(nrow(ibs))) {
  s1 <- ibs$sample1[i]
  s2 <- ibs$sample2[i]
  d <- ibs$IBS_distance[i]
  
  D[s1, s2] <- d
  D[s2, s1] <- d
}

# Hierarchical clustering
hc <- hclust(
  as.dist(D),
  method = "average"
)

clusters3 <- cutree(hc, k = 3)

cluster_check <- data.frame(
  sample = names(clusters3),
  cluster = as.integer(clusters3),
  stringsAsFactors = FALSE
)

cluster_check <- merge(
  cluster_check,
  assign,
  by = "sample"
)

cat("\nCluster vs genet assignment:\n")
print(table(cluster_check$cluster, cluster_check$genotype))

# Cut-height interval producing three clusters
lower <- hc$height[n - 3]
upper <- hc$height[n - 2]
midpoint <- (lower + upper) / 2

cat("\nThree-cluster interval:\n")
cat("lower =", lower, "\n")
cat("upper =", upper, "\n")

mid_clusters <- cutree(hc, h = midpoint)

if (length(unique(mid_clusters)) != 3) {
  stop("Midpoint does not produce three clusters.")
}

# Number of clusters across cut heights
heights <- seq(
  0,
  max(hc$height) * 1.02,
  length.out = 1000
)

nclusters <- sapply(
  heights,
  function(h) {
    length(unique(cutree(hc, h = h)))
  }
)

sensitivity <- data.frame(
  cut_height = heights,
  n_clusters = nclusters
)

# Figure
figure_S1 <- ggplot(
  sensitivity,
  aes(
    x = cut_height,
    y = n_clusters
  )
) +
  annotate(
    "rect",
    xmin = lower,
    xmax = upper,
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.07
  ) +
  geom_step(
    linewidth = 0.8
  ) +
  geom_hline(
    yintercept = 3,
    linetype = "dashed",
    linewidth = 0.45
  ) +
  geom_vline(
    xintercept = c(lower, upper),
    linetype = "dotted",
    linewidth = 0.45
  ) +
  annotate(
    "text",
    x = midpoint,
    y = 5,
    label = sprintf(
      "3-genet interval\n%.3f-%.3f",
      lower,
      upper
    ),
    size = 3.5
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_continuous(
    breaks = 1:10
  ) +
  coord_cartesian(
    ylim = c(1, 10)
  ) +
  labs(
    x = "Dendrogram cut height",
    y = "Number of clusters"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9)
  )

figure_S1

# Save figure
dir.create(
  "figures",
  showWarnings = FALSE
)

ggsave(
  filename = "figures/Figure_S1_clustering_sensitivity.pdf",
  plot = figure_S1,
  width = 6.5,
  height = 4.5,
  units = "in"
)

ggsave(
  filename = "figures/Figure_S1_clustering_sensitivity.png",
  plot = figure_S1,
  width = 6.5,
  height = 4.5,
  units = "in",
  dpi = 600
)