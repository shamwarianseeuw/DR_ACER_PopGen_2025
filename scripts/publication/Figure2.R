###############################################################################
# Figure 2: genomic evidence for multilocus genet assignment
###############################################################################

library(ggplot2)
library(ggdendro)
library(patchwork)

# Plot settings
genet_levels <- c("AC1", "AC2", "AC3")

genet_cols <- c(
  AC1 = "#0072B2",
  AC2 = "#D55E00",
  AC3 = "#009E73"
)

genet_shapes <- c(
  AC1 = 21,
  AC2 = 22,
  AC3 = 24
)

base_theme <- theme_classic(base_size = 11) +
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.tag = element_text(
      face = "bold",
      size = 14
    )
  )

# Genet assignments
geno <- read.delim(
  "metadata/genotype_assignments.tsv",
  stringsAsFactors = FALSE
)

geno$genotype <- factor(
  geno$genotype,
  levels = genet_levels
)

genet_lookup <- setNames(
  as.character(geno$genotype),
  geno$sample
)

# Panel A: PCoA
pcoa <- read.delim(
  "processed_data/ref_main_ibs_pca_coordinates.tsv",
  stringsAsFactors = FALSE
)

pcoa$genet <- factor(
  genet_lookup[pcoa$sample],
  levels = genet_levels
)

if (anyNA(pcoa$genet)) {
  stop("Some PCoA samples are missing genet assignments.")
}

# % variation explained
pcoa1_pct <- 46.86
pcoa2_pct <- 29.78

panel_A <- ggplot(
  pcoa,
  aes(
    x = PC1,
    y = PC2,
    fill = genet,
    shape = genet
  )
) +
  geom_point(
    size = 3.2,
    stroke = 0.55,
    color = "black"
  ) +
  scale_fill_manual(
    values = genet_cols,
    drop = FALSE
  ) +
  scale_shape_manual(values = genet_shapes,
    drop = FALSE
  ) +
  labs(
    x = sprintf("PCoA1 (%.2f%%)", pcoa1_pct),
    y = sprintf("PCoA2 (%.2f%%)", pcoa2_pct),
    fill = "Genet",
    shape = "Genet"
  ) +
  base_theme

# Panel B: hierarchical clustering
D <- as.matrix(
  read.table(
    "processed_data/ref_main.ibsMat",
    header = FALSE,
    check.names = FALSE
  )
)

bam_order <- read.table(
  "processed_data/ref_main_sample_order.txt",
  stringsAsFactors = FALSE
)[,1]

samples <- basename(bam_order)
samples <- sub("\\.bam$", "", samples)

if (
  nrow(D) != length(samples) ||
  ncol(D) != length(samples)
) {
  stop("IBS matrix dimensions do not match sample order.")
}

if (max(abs(D - t(D))) > 1e-10) {
  stop("IBS matrix is not symmetric.")
}

diag(D) <- 0
dimnames(D) <- list(samples, samples)

hc <- hclust(
  as.dist(D),
  method = "average"
)

clusters3 <- cutree(
  hc,
  k = 3
)

cluster_check <- data.frame(
  sample = names(clusters3),
  cluster = clusters3,
  genet = genet_lookup[names(clusters3)],
  stringsAsFactors = FALSE
)

cat("\nCluster vs genet assignment:\n")
print(
  table(
    cluster_check$cluster,
    cluster_check$genet
  )
)

# Dendrogram data
dend <- as.dendrogram(hc)

ddata <- dendro_data(
  dend,
  type = "rectangle"
)

segments <- segment(ddata)
labels_df <- label(ddata)

# Short sample labels
labels_df$label_short <- sub(
  "\\.trim\\.bt2$",
  "",
  labels_df$label
)

# Genet for each tip
labels_df$genet <- factor(
  genet_lookup[labels_df$label],
  levels = genet_levels
)

panel_B <- ggplot() +
  geom_segment(
    data = segments,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend
    ),
    linewidth = 0.35
  ) +
  geom_text(
    data = labels_df,
    aes(
      x = x,
      y = y - 0.006,
      label = label_short,
      color = genet
    ),
    angle = 90,
    hjust = 1,
    size = 2.6
  ) +
  scale_color_manual(
    values = genet_cols,
    drop = FALSE
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(0.18, 0.03)
    )
  ) +
  labs(
    x = NULL,
    y = "IBS distance"
  ) +
  guides(
    color = "none"
  ) +
  base_theme +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank()
  )

# Panel C: pairwise IBS distances
ibs <- read.delim(
  "processed_data/all_pairwise_IBS_distances.tsv",
  stringsAsFactors = FALSE
)

if (nrow(ibs) != choose(45, 2)) {
  stop("Expected 990 pairwise IBS comparisons.")
}

comparison_levels <- c(
  "within_AC1",
  "within_AC2",
  "within_AC3",
  "AC1_vs_AC2",
  "AC1_vs_AC3",
  "AC2_vs_AC3"
)

comparison_labels <- c(
  "AC1",
  "AC2",
  "AC3",
  "AC1-AC2",
  "AC1-AC3",
  "AC2-AC3"
)

# Plot labels
ibs$comparison_plot <- factor(
  ibs$comparison,
  levels = comparison_levels,
  labels = comparison_labels
)

# Check comparison labels
if (anyNA(ibs$comparison_plot)) {
  stop("Some IBS comparisons could not be assigned to plotting categories.")
}

panel_C <- ggplot(
  ibs,
  aes(
    x = comparison_plot,
    y = IBS_distance
  )
) +
  geom_boxplot(
    width = 0.58,
    outlier.shape = NA,
    linewidth = 0.45
  ) +
  geom_jitter(
    width = 0.13,
    height = 0,
    size = 0.7,
    alpha = 0.45
  ) +
  geom_vline(
    xintercept = 3.5,
    linetype = "dashed",
    linewidth = 0.35
  ) +
  annotate(
    "text",
    x = 2,
    y = 0.505,
    label = "Within genet",
    size = 3.2
  ) +
  annotate(
    "text",
    x = 5,
    y = 0.505,
    label = "Between genets",
    size = 3.2
  ) +
  coord_cartesian(
    ylim = c(0.14, 0.51)
  ) +
  labs(
    x = NULL,
    y = "Pairwise IBS distance"
  ) +
  base_theme +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    )
  )

# Combine panels
figure_2 <- (
  panel_A |
    panel_B |
    panel_C
) +
  plot_layout(
    widths = c(1, 1.65, 1.2),
    guides = "collect"
  ) +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    legend.position = "top"
  )

figure_2

# Save figure

dir.create(
  "figures",
  showWarnings = FALSE
)
 
ggsave(
   filename = "figures/Figure_2_genomic_evidence.pdf",
   plot = figure_2,
   width = 12,
   height = 4.8,
   units = "in"
 )

 ggsave(
   filename = "figures/Figure_2_genomic_evidence.png",
   plot = figure_2,
   width = 12,
   height = 4.8,
   units = "in",
   dpi = 600
 )