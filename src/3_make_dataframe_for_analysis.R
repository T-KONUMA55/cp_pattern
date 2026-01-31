library(ape)
library(tidyverse)
library(BHPMF)
source("src/utils.R")

# make gap-filling input ---------------------------------------------------

# bhpmf_input = make_bhpmf_input(df_target_mean, df_phylo_info, df_group)
# saveRDS(bhpmf_input, "./output/bhpmf/bhpmf_input.rds")
bhpmf_input <- readRDS("./output/bhpmf/bhpmf_input.rds")

# gap filling by BHPMF ----------------------------------------------------

set.seed(1)
dir.create("./output/bhpmf", recursive = TRUE, showWarnings = FALSE)
# eliminate or substitute '0'
GapFilling(bhpmf_input$trait, bhpmf_input$hierarchy,
           prediction.level = 4,
           used.num.hierarchy.levels = 3,
           mean.gap.filled.output.path = "./output/bhpmf/mean_gap_filled.csv",
           std.gap.filled.output.path = "./output/bhpmf/std_gap_filled.csv",
           verbose = FALSE)

# make input file ---------------------------------------------------------

col_tr = colnames(bhpmf_input$trait)
# download tree data from "http://dx.doi.org/10.5061/dryad.63q27"
tree = read.tree("./input/zanne2014.new")
input_glm = make_input_glm(bhpmf_input$hierarchy,
                           "./output/bhpmf/mean_gap_filled.csv",
                           sample %>% filter(species_try != 'Linum usitatissimum'), # Linum usitatissimum is excluded because it has a different mode of inheritance in the paper
                           col_tr,
                           tree,
                           standardize = TRUE)
write_csv(input_glm, './output/input_glm.csv')

rm(bhpmf_input, df_target_mean, df_phylo_info, df_group)
