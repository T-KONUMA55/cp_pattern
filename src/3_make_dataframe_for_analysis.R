library(ape)
library(tidyverse)
library(BHPMF)
source("src/utils.R")


# make gap-filling input ---------------------------------------------------


bhpmf_input = make_bhpmf_input(df_target_mean, df_phylo_info, df_group)


# gap filling by BHPMF ----------------------------------------------------
set.seed(42)
# eliminate or substitute '0'
GapFilling(bhpmf_input$trait, bhpmf_input$hierarchy,
           prediction.level = 4,
           used.num.hierarchy.levels = 3,
           mean.gap.filled.output.path = "./output/bhpmf/mean_gap_filled.csv",
           std.gap.filled.output.path = "./output/bhpmf/std_gap_filled.csv")


# make input file ---------------------------------------------------------

col_tr = colnames(bhpmf_input$trait)
# download tree data from "http://dx.doi.org/10.5061/dryad.63q27"
tree = read.tree("./input/zanne2014.new")
input_glm = make_input_glm(bhpmf_input$hierarchy,
                           "./output/bhpmf/mean_gap_filled.csv",
                           sample %>% filter(AccSpeciesName != 'Linum usitatissimum'), # Linum usitatissimum is excluded because it has a different mode of inheritance in the paper
                           col_tr,
                           tree,
                           standardize = TRUE)
write_csv(input_glm, './output/input_glm.csv')
