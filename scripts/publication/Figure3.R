###############################################################################
# Figure 3: distribution of multilocus genets among nurseries
###############################################################################

library(ggplot2)

# Plot settings
genet_levels <- c("AC1", "AC2", "AC3")

genet_cols <- c(
  AC1 = "#0072B2",
  AC2 = "#D55E00",
  AC3 = "#009E73"
)

base_theme <- theme_classic(base_size = 11) +
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

# Sample data
nursery <- read.delim(
  "metadata/sample_nursery.tsv",
  stringsAsFactors = FALSE
)

genets <- read.delim(
  "metadata/genotype_assignments.tsv",
  stringsAsFactors = FALSE
)

# Standardize IDs
if ("Tube_ID" %in% names(nursery)) {
  names(nursery)[names(nursery) == "Tube_ID"] <- "sample"
}

if ("genotype" %in% names(genets)) {
  names(genets)[names(genets) == "genotype"] <- "genet"
}

# Convert N100.trim.bt2 -> N100
genets$sample <- sub(
  "\\.trim\\.bt2$",
  "",
  genets$sample
)

# Join metadata and genet assignments
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
  stop("Missing nursery or genet assignments after join.")
}

dat$genet <- factor(
  dat$genet,
  levels = genet_levels
)

dat$Nursery <- factor(
  dat$Nursery,
  levels = c("Cap Cana", "Acuario")
)

# Genet counts by nursery
genet_counts <- as.data.frame(
  table(
    Nursery = dat$Nursery,
    Genet = dat$genet
  ),
  stringsAsFactors = FALSE
)

names(genet_counts)[3] <- "n"

genet_counts$Nursery <- factor(
  genet_counts$Nursery,
  levels = c("Cap Cana", "Acuario")
)

genet_counts$Genet <- factor(
  genet_counts$Genet,
  levels = genet_levels
)

nursery_totals <- aggregate(
  n ~ Nursery,
  data = genet_counts,
  FUN = sum
)

names(nursery_totals)[2] <- "N"

genet_counts <- merge(
  genet_counts,
  nursery_totals,
  by = "Nursery"
)

genet_counts$percent <- 100 * genet_counts$n / genet_counts$N

# Restore plotting order after merge
genet_counts$Nursery <- factor(
  genet_counts$Nursery,
  levels = c("Cap Cana", "Acuario")
)

genet_counts$Genet <- factor(
  genet_counts$Genet,
  levels = genet_levels
)

cat("\nFIGURE 3 DATA\n\n")
print(genet_counts)

# Figure 3
figure_3 <- ggplot(
  genet_counts,
  aes(
    x = Nursery,
    y = n,
    fill = Genet
  )
) +
  geom_col(
    width = 0.68,
    color = "black",
    linewidth = 0.35,
    position = position_stack(reverse = TRUE)
  ) +
  geom_text(
    aes(
      label = ifelse(
        n == 0,
        "",
        paste0(Genet, "\n", n)
      )
    ),
    position = position_stack(vjust = 0.5, reverse = TRUE),
    size = 3.5
  ) +
  scale_fill_manual(
    values = genet_cols,
    drop = FALSE
  ) +
  scale_y_continuous(
    breaks = seq(0, 30, 5),
    expand = expansion(
      mult = c(0, 0.05)
    )
  ) +
  labs(
    x = NULL,
    y = "Number of sampled colonies",
    fill = "Genet"
  ) +
  base_theme +
  theme(
    axis.text.x = element_text(size = 11),
    legend.position = "top"
  )

figure_3

# Save figure
dir.create(
  "figures",
  showWarnings = FALSE
)

ggsave(
  filename = "figures/Figure_3_genet_distribution.pdf",
  plot = figure_3,
  width = 6,
  height = 4.5,
  units = "in"
)

ggsave(
  filename = "figures/Figure_3_genet_distribution.png",
  plot = figure_3,
  width = 6,
  height = 5,
  units = "in",
  dpi = 600
)