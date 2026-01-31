library(phytools)
library(tidyverse)
library(ggplot2)


# Plot top5 AIC distribution ---------------------------------------------------------

# Set the theme for the plot
theme_set(theme_classic(base_size = 18, base_family = "Helvetica"))

# Create violin plot with box plot overlay
p <- ggplot(rslt_aic_top5_rm_gymn, 
            aes(x = model, y = AIC)) +
  geom_violin(trim = FALSE, alpha = 0.7, fill = "lightgray") +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), 
               outlier.shape = NA, fill = "black", color = "black", alpha = 1) +
  stat_summary(fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  scale_x_discrete(limits = c("H + LA + LeafN + LNA + SCD",
                              "H + LA + LNA + SCD",
                              "H + LA + LMA + LNA + SCD",
                              "H + LA + LMA + LNA + SCD + SSD",
                              "LA + LDMC + LeafN + LMA"), 
                   labels = c("H + LA + LeafN +\n LNA + SCD",
                              "H + LA +\n LNA + SCD",
                              "H + LA + LMA +\n LNA + SCD",
                              "H + LA + LMA +\n LNA + SCD + SSD",
                              "LA + LDMC +\n LeafN + LMA")
  ) +
  scale_y_continuous(limits = c(260, 270)) +
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
ggsave(filename = "output/plot/aic_violin_box_rm_gymn.png", plot = p, width = 8.9, height = 6, dpi = 600)


# Plot coefficient distribution -------------------------------------------------------------

# Set the theme for the plot
theme_set(theme_classic(base_size = 18, base_family = "Helvetica"))

# Create violin plot with box plot overlay
p <- ggplot(filter(rslt_coef_top_rm_gymn, trait != "(Intercept)"), aes(x = trait, y = coefficient)) +
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
ggsave(filename = "output/plot/coef_violin_box_rm_gymn.png", plot = p, width = 8.9, height = 6, dpi = 600)
