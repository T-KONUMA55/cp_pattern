library(tidyverse)
library(purrr)
source("src/utils.R")

# make trait list (trait_name = TraitID) ----------------------------------


df_rm_Dataname = read_csv("input/remove_Dataname.csv")
df_unit_conv = read_csv("input/unit_convert_factor.csv")
trait_list = list( 
  LA = 3109, # Leaf area (in case of compound leaves: leaflet, petiole excluded)
  LMA = 3115, # Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole excluded
  LFM = 163, # Leaf fresh mass
  LDMC = 47, # Leaf dry mass per leaf fresh mass (leaf dry matter content, LDMC)
  LeafC = 13, # Leaf carbon (C) content per leaf dry mass
  LeafN = 14, # Leaf nitrogen (N) content per leaf dry mass
  LeafP = 15, # Leaf phosphorus (P) content per leaf dry mass
  LNA = 50, # Leaf nitrogen (N) content per leaf area
  LNP = 56, # Leaf nitrogen/phosphorus (N/P) ratio
  d15N = 78, # Leaf nitrogen (N) isotope signature (delta 15N)
  SM = 26, # Seed dry mass
  SL = 27, # Seed length
  SNRU = 1103, # Seed number per reproductive unit
  DL = 237, # Dispersal unit length
  H = 3106, # Plant height vegetative
  SSD = 4, # Stem specific density (SSD) or wood density (stem dry mass per stem fresh volume)
  SCD = 169, # stem conduit density is the number of vessels and tracheids per unit area in a cross section,
  CEL = 282 # conduit element length refers to both vessels and tracheids
) # Bruelheide et al., 2018

# calculate the ratio of missing values (NA) for each trait cross species --------
filtered_trait_list = calc_na_ratio(df_numeric_trait, df_phylo_info, df_group, trait_list)


# extract target trait row ---------------------------------------------

df_target_mean = mean_target_trait(df_numeric_trait,
                                   filtered_trait_list,
                                   df_rm_Dataname, 
                                   df_unit_conv, 
                                   zero.omit = TRUE)

# rm(df_numeric_trait)
