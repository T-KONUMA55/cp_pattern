library(tidyverse)
source("src/utils.R")


# read data ---------------------------------------------------------------

# download trait data from "https://www.try-db.org/TryWeb/Home.php"
raw = read_tsv("input/TRY_DATA.txt", locale = locale(encoding = "ISO-8859-1"))
sample = read_csv("input/inheritance_pattern.csv", locale = locale(encoding = "ISO-8859-1"))
df_add_family = read_csv("input/genus_family.csv", locale = locale(encoding = "ISO-8859-1"))
df_replace_family = read_csv("input/family_replace.csv", locale = locale(encoding = "ISO-8859-1"))


# make species information (genus, family) file from TRY data -------------


df_phylo_info = make_phylo_info(raw, sample, df_add_family, df_replace_family)


# make group information (family, group) from TRY & sample data --------


family_gymnosperm = c("cycadaceae", "ginkgoaceae", "pinaceae", "taxodiaceae",
                      "cupressaceae", "araucariaceae", "podocarpaceae", "cephalotaxaceae",
                      "taxaceae", "welwitschiaceae", "ephedraceae", "gnetaceae")
df_group = make_group_info(raw, df_phylo_info, sample, df_replace_family, family_gymnosperm)


# remove na and extract numeric trait -----------------------------


df_numeric_trait = extract_numeric_trait(raw)
rm(raw)