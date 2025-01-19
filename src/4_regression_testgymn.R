library(ape)
library(tidyverse)
source("src/utils.R")

# Read data ------------------------------------------------

# Read the processed input trait data
input_glm_rm_gymn = read_csv("output/input_glm.csv")

# Select trait columns for correlation analysis
col_tr = colnames(input_glm_rm_gymn)[!(colnames(input_glm_rm_gymn) %in% c("AccSpeciesName", "pattern_cp", "pollen", "Genus", "Family", "Group", "PG", "d15N"))]

# Read the phylogenetic tree data (Zanne et al., 2014)
zanne = read.tree("input/zanne2014.new")

# Filter out gymnosperm species from the trait data
input_glm_rm_gymn = input_glm_rm_gymn %>% filter(Group != "Gymnosperm")


# Create a color gradient for the correlation plot ----------------------------

# Define a color palette for visualizing the correlation coefficients
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

# Generate and customize the correlation plot
corrplot(cor(input_glm_rm_gymn %>% select(col_tr)),
                       method = "shade",
                       shade.col = NA,
                       tl.col = "black",
                       tl.srt = 45, 
                       col = col(200),
                       addCoef.col = "black", 
                       addcolorlabel = "no", 
                       order = "AOE")

# Perform Regression Analysis ------------------------------------------------


# Get the number of physical cores (ignoring logical/virtual cores)
detectCores(logical = FALSE)

# Run model fitting with multiple seeds
rslt_aic_all_rm_gymn <- pbmclapply(1:100, function(i) {
  phy <- edit_new_rep(input_glm_rm_gymn, zanne, seed = i)
  phyloglm_exh_test(input_glm_rm_gymn, phy, col_tr)
}, mc.cores = 4)

# Calculate summary statistics (mean and median AIC) for each model
rslt_aic_sum_rm_gymn = mkdf(rslt_aic_all_rm_gymn) %>% 
  group_by(model) %>% 
  summarise(mean = mean(AIC), median = median(AIC))

# Identify the top 5 models with the lowest median AIC values
top5_rm_gymn = rslt_aic_sum_rm_gymn %>% arrange(median) %>% head(5)

# Filter results to include only the top 5 models
rslt_aic_top5_rm_gymn = mkdf(rslt_aic_all_rm_gymn) %>% 
  filter(model %in% top5_rm_gymn$model)


# Calculate coefficients ------------------------------------------------------------


# Arrange the results by median AIC and extract the top model's name
arrange(rslt_aic_sum_rm_gymn, median)[[1,1]]
# Replace spaces in species names with underscores to match tree tip labels
rownames(input_glm_rm_gymn) = str_replace_all(input_glm_rm_gymn$AccSpeciesName, pattern = " ", replacement = "_")

# Perform phylogenetic logistic regression for 100 iterations with parallel processing
rslt_coef_all_rm_gymn = pbmclapply(1:100, function(i) {
  # Create a new phylogenetic tree by randomly selecting one species per genus
  phy = edit_new_rep(input_glm_rm_gymn, zanne, seed = i)
  
  # Fit a phylogenetic logistic regression model using the selected traits
  rslt = phyloglm(input_glm_rm_gymn$pattern_cp ~ H + LA + LeafN + LMA + SCD,
                  data = input_glm_rm_gymn,
                  phy = phy,
                  method = c("logistic_MPLE","logistic_IG10", "poisson_GEE"),
                  btol = 10,
                  log.alpha.bound = 4,
                  start.beta = NULL,
                  start.alpha = NULL,
                  boot = 0,
                  full.matrix = TRUE)
  
  # Extract coefficients from the model and store them with their respective traits
  return(data.frame(trait = names(rslt$coefficients), coefficient = rslt$coefficients, row.names = NULL))
}, mc.cores = 4) # Use 4 CPU cores for parallel computation

# Summarize coefficients across iterations ---------------------------------

# Summarize the coefficients for each trait by calculating mean and median
rslt_coef_sum_rm_gymn = mkdf(rslt_coef_all_rm_gymn) %>%
  group_by(trait) %>% 
  summarise(mean = mean(coefficient), median = median(coefficient))

# Create a data frame of all coefficients across iterations
rslt_coef_top_rm_gymn = mkdf(rslt_coef_all_rm_gymn)

# Summarize the median coefficient for each trait
result_rm_gymn <- rslt_coef_top_rm_gymn %>%
  group_by(trait) %>%
  summarise(median_coefficient = median(coefficient))