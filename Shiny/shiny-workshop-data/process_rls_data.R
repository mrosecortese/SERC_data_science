# Process data from https://doi.org/10.25573/serc.20344158.v1
# Data/R script accompanying Blinded by the bright publication

# This script prepares the data for inclusion in a Shiny application
# Converting to long-form and dropping 33 observations where more resolved IDs were present in the data

library(tidyverse)

# Import the survey data
rls_data_in <- read_csv("Shiny/shiny-workshop-data/blinded-by-the-bright-fighshare/CBC_RLS_2015-2019_taxonomic.csv")

# Pivot from wide to long form, clean up data
rls_data <- rls_data_in %>%
  select(-c(total_count, valid_name, -`...1`)) %>%
  pivot_longer(invert_count:`400`,
               names_to = "size_class",
               values_to = "size_count") %>%
  # Remove data not of interest to the analysis
  # Blueground and Tobacco seagrass removed because each only surveyed once
  filter(size_count > 0,
         method == 1,
         phylum == "Chordata",
         !location_name %in% c("Blueground Seagrass", "Tobacco Seagrass")) %>%
  mutate(habitat = recode_factor(habitat,
    "forereef" = "Fore Reef",
    "patch reef" = "Patch Reef",
    "mangrove" = "Mangrove",
    "sand" = "Sand",
    "seagrass" = "Seagrass"
  ))

species_df <- rls_data %>%
  filter(rank == "Species")

genus_df <- rls_data %>%
  filter(rank == "Genus")

family_df <- rls_data %>%
  filter(rank == "Family")

rls_data_out <- species_df %>%
  # Only include genus-level observations if that genus is not found in the species data
  bind_rows(
    genus_df %>%
      filter(!(genus %in% species_df$genus))
  ) %>%
  # Only include family-level observations if that family is not found in species or genus data
  bind_rows(
    family_df %>%
      filter(!(family %in% species_df$family),
             !(family %in% species_df$family))
  )

# Write out for use by the RMarkdown script
write_csv(rls_data_out, "Shiny/shiny-workshop-data/cbc_rls_shiny_workshop_data.csv")

# Write out for use by the example application
write_csv(rls_data_out, "Shiny/app-example/resources/cbc_rls_shiny_workshop_data.csv")
