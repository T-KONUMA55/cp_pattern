library(tidyverse)
library(ggplot2)
library(convertr)
library(BHPMF)
library(gtools)
library(tcltk)
library(ape)
library(phylolm)
library(corrplot)
library(parallel)
library(pbmcapply)

# Make species information (genus, family) file from TRY & sample  --------


make_phylo_info = function(df_TRY, df_sample, df_add_family, df_replace_family) {
  # Extract family information from TRY data and format it by genus
  family_TRY = df_TRY %>%
    filter(grepl("Family|family", DataName)) %>%
    separate(AccSpeciesName,into = c("Genus", "Species"), sep = " ", remove = FALSE) %>%
    select(Genus, OrigValueStr)
  
  # Format sample data by genus
  family_smpl = df_sample %>%
    mutate(Genus = genus, OrigValueStr = family) %>%
    select(Genus, OrigValueStr)
  
  # Combine TRY and sample data, fill in missing Family information using df_add_family
  family_all = bind_rows(family_smpl, filter(family_TRY, Genus %in% setdiff(family_TRY$Genus, family_smpl$Genus))) %>%
    left_join(df_add_family, by = c("Genus" = "Genus")) %>%
    mutate(OrigValueStr_replace = ifelse(is.na(Family), OrigValueStr, Family)) %>% 
    distinct(Genus, .keep_all = TRUE) %>% 
    select(-Family)
  
  
  # Extract unique species information from TRY data
  species_TRY = df_TRY %>%
    select(AccSpeciesID, AccSpeciesName) %>%
    distinct(AccSpeciesName, .keep_all = TRUE) %>%
    separate(AccSpeciesName,into = c("Genus", "Species"), sep = " ", remove = FALSE)
  
  # Check for duplicate species information
  print(species_TRY %>%
          group_by(AccSpeciesName) %>%
          filter(n() > 1))
  
  # Merge species information with family information, apply replacements for family names
  df_phylo = merge(species_TRY, family_all, by = "Genus", all.x = TRUE) %>% 
    drop_na(OrigValueStr) %>%
    arrange(AccSpeciesID) %>%
    mutate(Genus_lower = tolower(Genus),
           Family_lower = tolower(OrigValueStr)) %>%
    left_join(df_replace_family, by = c("Family_lower" = "Family")) %>%
    mutate(Family_replace = ifelse(is.na(replacement), Family_lower, replacement)) %>%
    select(AccSpeciesID, AccSpeciesName, Genus_lower, Family_replace)
  colnames(df_phylo)[3:4] <- c("Genus", "Family")
  
  return(df_phylo)
}

# Make group information (family, group) from TRY & sample data --------


make_group_info = function(df_TRY, df_phylo_info, df_sample, df_replace_family, family_gymnosperm) {
  # Define patterns to replace with 'Eudicot'
  pattern_to_replace = c("Eudicots",
                         "Core eudicot",
                         "Basal eudicot",
                         "Core Eudicots")
  
  # Process TRY data to extract genus and group information
  group_TRY = df_TRY %>%
    select(AccSpeciesName, OrigValueStr) %>%
    filter(grepl("Dicot|dicot|Monocot|monocot", OrigValueStr) & nchar(OrigValueStr) < 20) %>%
    separate(AccSpeciesName,into = c("Genus", "Species"), sep = " ", remove = TRUE) %>%
    distinct(Genus, .keep_all = TRUE) %>%
    mutate(Group = gsub(OrigValueStr, pattern = "Monocots", replacement = "Monocot", ignore.case = TRUE)) %>%
    mutate(Group = gsub(paste0("(", paste(pattern_to_replace, collapse = "|"), ")"), "Eudicot", Group, ignore.case = TRUE)) %>%
    mutate(Genus = tolower(Genus)) %>%
    merge(df_phylo_info, by = "Genus", all = TRUE) %>%
    select(Family, Group) %>%
    na.omit() %>%
    distinct(Family, .keep_all = TRUE) %>%
    filter(!Family %in% family_gymnosperm) %>%
    bind_rows(data.frame(Family = family_gymnosperm, Group = "Gymnosperm"))
  
  # Process sample data to extract family and group information
  group_smpl = df_sample %>%
    select(family, cotyledon) %>%
    na.omit() %>%
    mutate(Family_lower = tolower(family)) %>%
    left_join(df_replace_family, by = c("Family_lower" = "Family")) %>%
    mutate(Family_replace = ifelse(is.na(replacement), Family_lower, replacement),
           Group = cotyledon2) %>%
    select(Family_replace, Group) %>%
    distinct(Family_replace, .keep_all = TRUE)
  colnames(group_smpl)[1] = "Family"
  
  # Combine TRY and sample data for final group information
  df_group = bind_rows(group_TRY, group_smpl) %>%
    distinct(Family, .keep_all = TRUE)
  
  return(df_group)
}

# Make tidy trait-data frame from TRY & sample data -----------------------


extract_numeric_trait = function(df_TRY) {
  # Helper function to check if a value is numeric
  isNumeric <- function(st) {
    return(!is.na(suppressWarnings(as.numeric(st))))
  }
  
  # Filter and process the data frame to extract numeric traits
  df_numeric_trait = df_TRY %>%
    select(AccSpeciesName, AccSpeciesID, TraitID, TraitName, DataName, OriglName, OrigValueStr, OrigUnitStr, ValueKindName, Reference, Comment) %>%
    filter(!is.na(TraitID), !is.na(OrigValueStr), !is.na(OriglName), isNumeric(OrigValueStr)) %>%
    arrange(TraitID, OriglName, AccSpeciesName) %>% 
    mutate(OrigValueNum = as.numeric(OrigValueStr))
  
  return(df_numeric_trait)
}


# Calculate the ratio of missing values (NA) for each trait cross species --------


calc_na_ratio = function(df_numeric_trait, df_phylo_info, df_group, trait_list){
  # Merge data frames to incorporate phylogenetic and group information
  df_tmp = df_numeric_trait %>%
    inner_join(df_phylo_info, by = "AccSpeciesName") %>%
    inner_join(df_group, by = "Family") %>% 
    select(AccSpeciesName, Family, Genus, Group, TraitID, OrigValueNum) 
  
  # Filter the data to retain only the specified traits
  df_tmp = df_tmp[df_tmp$TraitID %in% unlist(trait_list), ]
  
  # Reshape the data into a wide format for NA calculation
  df_calc = df_tmp %>% 
    pivot_wider(names_from = TraitID, values_from = OrigValueNum, values_fill = NULL, values_fn = mean) 
  
  # Calculate NA ratio for each column
  colMeans(is.na(df_calc))
}


# Filter data for a specified trait and apply unit conversions --------


filter_trait_record = function(df_numeric_trait, trait_list, trait_name, df_rm_Dataname, df_unit_conv, zero.omit = TRUE) {
  # Extract TraitID based on the trait name
  trait_id = trait_list[[trait_name]]
  # Check if trait_id is NULL
  if (is.null(trait_id)) {
    stop(paste(trait_name , "not found in trait_list."))
  }
  
  # Filter the data for the specified trait
  df_trait = df_numeric_trait %>%
    filter(TraitID == trait_id)
  
  # If no rows are found, throw an error
  if (nrow(df_trait) == 0) {
    stop(paste("TraitID ", trait_id, "is not found in df_numeric_trait."))
  }
  # Handle zero.omit parameter
  if (zero.omit) {
    df_trait = df_trait[df_trait$OrigValueNum != 0, ]
  } else if (!zero.omit) {
  } else {
    stop("A non-Boolean type is specified for zero.omit")
  }
  
  # Filter and preprocess DataName and unit conversion information
  rm_filter = df_rm_Dataname %>% 
    filter(TraitID == trait_id)
  df_unit_conv_filter = df_unit_conv %>% 
    filter(trait == trait_name)
  
  # Filter and convert data
  df_trait$OrigUnitStr = ifelse(is.na(df_trait$OrigUnitStr), "nounit", df_trait$OrigUnitStr)
  df_filtered = df_trait %>%
    filter(!(DataName %in% rm_filter$DataName) & !ValueKindName %in% c("Maximum", "Minimum", "Median")) %>% 
    left_join(df_unit_conv_filter, by = c("OrigUnitStr" = "unit")) %>%
    mutate(ConvertedValue = case_when(
      TraitID == 47 & OrigUnitStr %in% c("gH2O / 100 g fresh mass", "gr H2O / FW*100") ~ (100 - OrigValueNum) / 100,
      substr(factor, 1, 3) == "inv" ~ (1 / OrigValueNum) * as.numeric(substr(factor, 4, nchar(factor))),
      TraitID == 78 ~ OrigValueNum,
      TRUE ~ OrigValueNum * as.numeric(factor)
    ))
  
  # Handle missing values in ConvertedValue
  if (any(is.na(df_filtered$ConvertedValue))) {
    na_rows = df_filtered %>% 
      filter(is.na(ConvertedValue))
    print(na_rows)
    print(na_rows[['OrigUnitStr']] %>% unique())
    stop(paste("NA values found in column ConvertedValue of TraitID ", trait_id))
  }
  
  return(df_filtered)
}


# Make hierarchy file for BHPMF from mk_gapdf output ----------------------


mean_target_trait = function(df_numeric_trait, trait_list, trait_name, df_rm_Dataname, df_unit_conv, zero.omit = TRUE) {
  # Initialize an empty list to store results
  output_list = list()
  
  # Apply `filter_trait_record` function to each trait in `trait_list`
  output_list = map(names(trait_list), function(trait_name) {
    filter_trait_record(df_numeric_trait, trait_list, trait_name, df_rm_Dataname, df_unit_conv)
  })
  
  # Combine all DataFrames from the list into one DataFrame
  target_rows_conv = reduce(output_list, bind_rows)
  
  # Calculate the mean value for each trait by species
  df_target_mean = target_rows_conv %>%
    select(AccSpeciesName, trait, ConvertedValue) %>% 
    pivot_wider(names_from = trait, values_from = ConvertedValue, values_fill = NULL, values_fn = mean) %>%
    arrange(AccSpeciesName)
  
  return(df_target_mean)
}


# Make gap-filling input --------------------------------------------------


make_bhpmf_input = function(df_target_mean, df_phylo_info, df_group) {
  # Log-transform trait values
  df_target_mean_log = df_target_mean %>%
    mutate(across(-AccSpeciesName, ~ ifelse(.x != 0, log10(.x), NA))) %>%
    inner_join(df_phylo_info, by = "AccSpeciesName") %>%
    inner_join(df_group, by = "Family") %>%
    separate(AccSpeciesName, into = c("Genus", "Species"), sep = " ", remove = FALSE) %>%
    arrange(AccSpeciesName) %>%
    mutate(plant_id = row_number())
  
  # Check for Inf/-Inf in the log-transformed data
  if (any(is.infinite(unlist(df_target_mean_log)))) {
    cat("The log-transformed data contains Inf or -Inf values.\n")
  } else {
    cat("The log-transformed data does not contain Inf or -Inf values.\n")
  }
  
  # Create a hierarchy data frame with classification information
  hierarchy.info = df_target_mean_log %>%
    select(plant_id, AccSpeciesName, Genus, Family, Group) %>% 
    data.frame()
  
  # Prepare a matrix of trait values
  trait.info = df_target_mean_log %>%
    select(names(trait_list)) %>%
    mutate(across(everything(), ~ ifelse(.x == 0, 0.0000001, .x))) %>%
    as.matrix()
  
  # Calculate the ratio of non-missing values for each trait
  trait_ratio = apply(trait.info, 2, function(x){length(x[!is.na(x)]) / nrow(trait.info)})
  
  # Generate a warning if any trait has more than 95% missing data
  if (any(trait_ratio < 0.05)) {
    print(trait_ratio)
    warning("Traits with missing rates higher than 0.95 are included.")
  }
  
  return(list(trait = trait.info, hierarchy = hierarchy.info))
}


# Make input file ---------------------------------------------------------

# Standardizes specified columns (col_trait) in a data frame by applying Z-transformation
ztrans = function(df, col_trait){
  input_z = df
  for (c in col_trait) {
    mu = mean(df[[c]], na.rm = TRUE)
    sigma = sd(df[[c]], na.rm = TRUE)
    input_z[c] = (df[[c]] - mu)/sigma
  }
  return(input_z)
}

make_input_glm = function(df_hierarchy, path_filled_trait, sample, col_trait, tree, standardize = TRUE) {
  # Read the gap-filled trait data and bind it with the hierarchical data
  df_gapfilled = df_hierarchy %>%
    bind_cols(read_tsv(path_filled_trait))
  
  # Assign gap-filled traits to the sample and preprocess variables
  sample_gapfilled = sample %>%
    inner_join(df_gapfilled, by = "AccSpeciesName") %>% 
    select(AccSpeciesName, pattern_cp, Genus, Family, Group, col_trait) %>%
    mutate(pattern_cp = str_replace_all(pattern_cp, pattern = c("B/M" = "1", "B" = "1", "RP" = "1", "P" = "1", "M" = "0")),
           PG = str_replace_all(Group, pattern = c("Eudicot" = "1", "dicot" = "1", "Monocot" = "1", "Gymnosperm" = "0"))) %>%
    arrange(Family)
  # Convert pattern_cp and PG columns to numeric
  sample_gapfilled$pattern_cp = as.numeric(sample_gapfilled$pattern_cp)
  sample_gapfilled$PG = as.numeric(sample_gapfilled$PG)
  
  # Get a list of unique genera in the hierarchical data
  genus_sample = df_hierarchy$Genus %>% unique()
  
  # Identify genera not found in the phylogenetic tree's tip labels
  rm_genus = setdiff(genus_sample, str_replace(tree$tip.label, pattern = "_.*$", replacement = "") %>% unique())
  
  # Remove samples belonging to the excluded genera and ensure uniqueness
  sample_gapfilled = sample_gapfilled %>%
    filter(!Genus %in% rm_genus) %>%
    unique()
  
  # Replace spaces in species names with underscores to standardize row names
  rownames(sample_gapfilled) = str_replace_all(sample_gapfilled$AccSpeciesName, pattern = " ", replacement = "_")
  
  # Print the number of available and removed species
  available_species_count = nrow(sample_gapfilled)
  removed_species_count = nrow(sample) - available_species_count
  cat(paste(available_species_count, "species available\n", removed_species_count, "species removed\n"))
  
  # Standardize the trait columns if specified
  if (standardize == TRUE) {
    sample_gapfilled_trans = ztrans(sample_gapfilled, col_trait)
    return(sample_gapfilled_trans)
  } else if (standardize == FALSE) {
    return(sample_gapfilled)
  } else {
    stop('Specify a bool type for argument "standardize"')
  }
}


# Function to edit and reassign species in a phylogenetic tree based on input -----------

edit_new_rep <- function(input, phy, seed) {
  set.seed(seed) # Set the random seed for reproducibility
  
  # Extract species names from the input
  sp_assign <- input[[1]] 
  
  # Get unique genera from the input by extracting genus names from species names
  genus_phy <- unique(str_replace(sp_assign, " .*$", "_"))
  
  # Create a data frame with genus and species extracted from sample names
  sp_phy <- data.frame(sample_name = sp_assign) %>%
    separate(sample_name, into = c("genus", "species"), sep = " ", remove = FALSE, extra = "drop")
  
  # Randomly assign one species per genus in the phylogenetic tree
  sp_zan <- sapply(genus_phy, function(gns) {
    sp <- str_subset(phy$tip.label, pattern = gns) %>% # Filter species in the tree that match the genus
      sample(1) # Randomly select one species
    return(sp)
  })
  
  # Prune the phylogenetic tree to retain only the selected species
  zanne_prune <- drop.tip(phy, phy$tip.label[!phy$tip.label %in% sp_zan])
  # Simplify tip labels by removing species names, retaining only genus
  zanne_prune$tip.label <- str_replace(zanne_prune$tip.label, "_.*$", "")
  
  # Convert the pruned tree into Newick format for editing
  zanne_prune_new <- write.tree(zanne_prune)
  
  # Replace genus names in the tree with their corresponding species names
  for (gns in zanne_prune$tip.label) {
    # Get all species in the genus that are present in the sample
    sp_in_gns <- sp_phy %>%
      filter(genus == gns) %>% # Filter rows matching the genus
      pull(sample_name) %>%    # Extract species names
      str_replace_all(" ", "_") # Replace spaces with underscores
    
    nm_sp <- length(sp_in_gns) # Count the number of species in the genus
    
    # Replace genus name in the tree with corresponding species names
    if (nm_sp == 1) {
      # Single species: replace genus with the species name
      zanne_prune_new <- str_replace(
        zanne_prune_new, 
        pattern = str_c(gns, ":", sep = ""), 
        replacement = str_c(sp_in_gns, ":", sep = "")
      )
    } else if (nm_sp > 1) {
      # Multiple species: create a sub-tree with the species
      zanne_prune_new <- str_replace(
        zanne_prune_new, 
        pattern = str_c(gns, ":", sep = ""), 
        replacement = str_c("(", str_c(sp_in_gns, collapse = ":0.000001,"), ":0.000001):", sep = "")
      )
    } else {
      # No species detected in the genus
      print("Genus which have no species in samples were detected.")
    }
  }
  
  # Return the modified phylogenetic tree by reading the edited Newick string
  return(read.tree(text = zanne_prune_new))
}

# Perform exhaustive phylogenetic logistic regression with all combinations of traits -----------------

phyloglm_exh <- function(input, phy, trait) {
  # Initialize an empty data frame to store AIC results for each model
  df_AIC <- data.frame(model = character(), AIC = numeric(), stringsAsFactors = FALSE)
  
  # Loop over the number of traits to create models with 1 to n traits
  for (r in 1:length(trait)) {
    # Generate all possible combinations of 'r' traits from the provided list
    cmb <- combinations(n = length(trait), r = r, v = trait, repeats.allowed = FALSE)
    
    # Perform phylogenetic logistic regression for each combination in parallel
    results <- mclapply(1:nrow(cmb), function(m) {
      # Select the current combination of traits and create an input data frame
      input_edit <- select(input, pattern_cp, cmb[m, ])
      # Set row names to match the species names in the phylogenetic tree
      rownames(input_edit) <- str_replace_all(input$AccSpeciesName, " ", "_")
      
      # Fit a phylogenetic logistic regression model
      result <- phyloglm(pattern_cp ~ ., 
                         data = input_edit, 
                         phy = phy, 
                         method = c("logistic_MPLE", "logistic_IG10", "poisson_GEE"),
                         btol = 10, 
                         log.alpha.bound = 4, 
                         start.beta = NULL, 
                         start.alpha = NULL, 
                         boot = 0, 
                         full.matrix = TRUE)
      
      # Create a model identifier by concatenating the trait names in the combination
      mdl <- str_flatten(cmb[m, ], collapse = " + ")
      # Return a data frame with the model identifier and its AIC value
      return(data.frame(model = mdl, AIC = result$aic, stringsAsFactors = FALSE))
    }, mc.cores = detectCores() - 1)
    
    # Combine the results from the current iteration with the overall results
    df_AIC <- bind_rows(df_AIC, do.call(rbind, results))
  }
  
  return(df_AIC)
}

# Convert a list of data frames into a single data frame ------------------------

mkdf = function(result_phylo) {
  # Initialize an empty data frame with two columns
  out = data.frame(matrix(rep(NA, 2), nrow = 1))[numeric(0), ]
  
  # Iterate over each element in the input list
  for (n in result_phylo) {
    # Append each data frame in the list to the output data frame
    out = rbind(out, n)
  }
  
  # Ensure the second column is numeric
  out[2] = as.numeric(out[[2]])
  
  # Return the combined data frame
  return(out)
}

