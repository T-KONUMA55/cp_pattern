library(ape)
library(tidyverse)
source("src/utils.R")

# Read data ------------------------------------------------

# Read the processed input trait data
input_glm = read_csv("output/input_glm.csv")

# Select trait columns for correlation analysis
col_tr = colnames(input_glm)[!(colnames(input_glm) %in% c("AccSpeciesName", "pattern_cp", "pollen", "Genus", "Family", "Group", "PG"))]

# Read the phylogenetic tree data (Zanne et al., 2014)
zanne = read.tree("input/zanne2014.new")


# Create a color gradient for the correlation plot ----------------------------

# Define a color palette for visualizing the correlation coefficients
col = colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

# Generate and customize the correlation plot
corrplot(cor(input_glm %>% select(col_tr)),
                       method = "shade",
                       shade.col = NA,
                       tl.col = "black",
                       tl.srt = 45, 
                       col = col(200),
                       addCoef.col = "black", 
                       addcolorlabel = "no", 
                       order = "AOE") 

# Remove traits whose correlation coefficient is 0.8 or higher
col_tr_filtered <- col_tr[col_tr != "DL"]

# Perform Regression Analysis ------------------------------------------------


# Get the number of physical cores (ignoring logical/virtual cores)
detectCores(logical = FALSE)
dir.create("./output/tree", recursive = TRUE, showWarnings = FALSE)

# Run model fitting with multiple seeds
make_tree = TRUE
# If make_tree = TRUE: generate a new phylogenetic tree.
# If make_tree = FALSE: load a saved tree from tree_path.
rslt_aic_all <- pbmclapply(1:100, function(i) {
  tree_path <- sprintf("./output/tree/phy_%03d.rds", i)
  if (make_tree) {
    phy <- edit_new_rep(input_glm, zanne, seed = i)
    saveRDS(phy, tree_path)
  } else {
    phy <- readRDS(tree_path)
  }
  phyloglm_exh(input_glm, phy, col_tr_filtered)
}, mc.cores = 4)
saveRDS(rslt_aic_all, "./output/glm_result.rds")

# Calculate summary statistics (mean and median AIC) for each model
rslt_aic_all <- readRDS("./output/glm_result.rds")
rslt_aic_sum = mkdf(rslt_aic_all) %>% 
  group_by(model) %>% 
  summarise(mean = mean(AIC), median = median(AIC))

# Identify the top 5 models with the lowest median AIC values
top5 = rslt_aic_sum %>% arrange(median) %>% head(5)

# Filter results to include only the top 5 models
rslt_aic_top5 = mkdf(rslt_aic_all) %>% 
  filter(model %in% top5$model)


# Calculate coefficients ------------------------------------------------------------


# Arrange the results by median AIC and extract the top model's name
arrange(rslt_aic_sum, median)[[1,1]]
# Replace spaces in species names with underscores to match tree tip labels
rownames(input_glm) = str_replace_all(input_glm$AccSpeciesName, pattern = " ", replacement = "_")

# Perform phylogenetic logistic regression for 100 iterations with parallel processing
rslt_coef_all = pbmclapply(1:100, function(i) {
  # Load a saved tree from tree_path
  tree_path <- sprintf("./output/tree/phy_%03d.rds", i)
  phy <- readRDS(tree_path)
  
  # Fit a phylogenetic logistic regression model using the selected traits
  rslt = phyloglm(input_glm$pattern_cp ~ H + LA + LNP + SCD,
                  data = input_glm,
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
rslt_coef_sum = mkdf(rslt_coef_all) %>%
  group_by(trait) %>%
  summarise(mean = mean(coefficient), median = median(coefficient))

# Create a data frame of all coefficients across iterations
rslt_coef_top = mkdf(rslt_coef_all)

# Summarize the median coefficient for each trait
result <- rslt_coef_top %>%
  group_by(trait) %>%
  summarise(median_coefficient = median(coefficient))