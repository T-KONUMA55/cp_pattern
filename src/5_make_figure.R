library(phytools)
library(tidyverse)
library(ggplot2)


# Ensure directory exists -------------------------------------------------

dir.create("./output/plot", recursive = TRUE, showWarnings = FALSE)

# Write phylogenetic tree of sample ---------------------------------------------

# Filter the phylogenetic tree to include only species present in `input_glm`
zanne_filtered = drop.tip(zanne, setdiff(zanne$tip.label, rownames(input_glm)))

# Save the plot as a PNG file
png("./output/plot/phylo_tree.png",
    width = 1840,
    height = 1840,
    )

# Extract the species names with a pattern_cp value of 1
selected_species = rownames(input_glm)[input_glm$pattern_cp == 1]

# Identify the tips in the filtered tree that correspond to the selected species
co_cp = zanne_filtered$tip.label %in% selected_species

# Set the default color of all tip labels to black
tip_colors = rep("black", length(zanne_filtered$tip.label))

# Change the color of the selected species to red
tip_colors[co_cp] = "red"

# Plot the phylogenetic tree
plot.phylo(zanne_filtered,
           type = "fan",
           cex = 1.2,
           tip.color = tip_colors,
           show.tip.label = TRUE
)

# Close the PNG device and save the file
dev.off()

# Plot trait data distribution -------------------------------------------------

# Set a consistent theme for all plots
theme_set(theme_bw(base_size = 18, base_family = "Helvetica"))

# Convert data to long format for visualization
input_long <- input_glm %>%
  pivot_longer(cols = c("H", "LA", "LDMC", "LeafC", "LeafN", "LeafP", "LMA", "LNA", "LNP", "SCD", "SM", "SSD"), 
               names_to = "Trait", 
               values_to = "Value") %>%
  mutate(pattern_cp = factor(pattern_cp, levels = c(0, 1), labels = c("M", "P/B")))

# Create violin plot with box plot overlay
p = ggplot(input_long, aes(x = pattern_cp, y = Value, fill = pattern_cp)) +
  geom_violin(alpha = 0.7, position = position_dodge(width = 0.9)) +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), 
               outlier.shape = NA, fill = "black", color = "black", alpha = 1) +
  stat_summary(fun = median, geom = "point", fill = "white", shape = 21, size = 3) +
  xlab("Inheritance pattern") +
  ylab("Value") +
  scale_fill_manual(values = c("lightblue", "lightgreen"), name = "Pattern CP") +
  theme(
    axis.text = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.x = element_text(size = 18),
    strip.background = element_blank(), 
  ) +
  facet_wrap(~ Trait, scales = "free_y")

# Display the plot
print(p)

# Save the plot as a PNG file
ggsave("output/plot/violin_boxplot.png", plot = p, width = 18.4, height = 8, dpi = 1000)


# Plot top5 AIC distribution ---------------------------------------------------------

# Set the theme for the plot
theme_set(theme_classic(base_size = 18, base_family = "Helvetica"))

# Create violin plot with box plot overlay
p <- ggplot(rslt_aic_top5, 
            aes(x = model, y = AIC)) +
  geom_violin(trim = FALSE, alpha = 0.7, fill = "lightgray") +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9),
               outlier.shape = NA, fill = "black", color = "black", alpha = 1) +
  stat_summary(fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  scale_x_discrete(limits = c("H + LA + LNP + SCD",
                              "H + LA + LeafN + LNP + SCD",
                              "H + LA + LeafP + LNP + SCD",
                              "H + LA + LeafN + LeafP + SCD",
                              "H + LA + LNA + LNP + SCD"), 
                   labels = c("H + LA +\n LNP + SCD",
                              "H + LA + LeafN +\n LNP + SCD",
                              "H + LA + LeafP +\n LNP + SCD",
                              "H + LA + LeafN +\n LeafP + SCD",
                              "H + LA + LNA +\n LNP + SCD")
                   ) +
  scale_y_continuous(limits = c(268, 273)) +
  ylab("AIC") +
  theme(
    axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
  ) 

# Display the plot
print(p)

# Save the plot as a PNG file
ggsave(filename = "output/plot/aic_violin_box_all.png", plot = p, width = 8.9, height = 6, dpi = 600)
 

# Plot coefficient distribution -------------------------------------------------------------

# Set the theme for the plot
theme_set(theme_classic(base_size = 18, base_family = "Helvetica"))

# Create violin plot with box plot overlay
p <- ggplot(filter(rslt_coef_top, trait != "(Intercept)"), aes(x = trait, y = coefficient)) +
  geom_violin(trim = FALSE, alpha = 0.7, fill = "lightgray") +
  geom_boxplot(width = 0.05, position = position_dodge(width = 0.9),
               outlier.shape = NA, fill = "black", color = "black", alpha = 1) +
  stat_summary(fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme(
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
  ) 

# Display the plot
print(p)

# Save the plot as a PNG file
ggsave(filename = "output/plot/coef_violin_box_all.png", plot = p, width = 8.9, height = 6, dpi = 600)
