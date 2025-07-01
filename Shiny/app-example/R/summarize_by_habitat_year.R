summarize_by_habitat_year <- function(df){

  # Loop over each year:
  unique_species_list <- lapply(unique(df$year), function(year){

    # One row per species observed in each habitat
    fish_year_subset <- df %>%
      filter(year == !!year) %>%
      distinct(habitat, scientific_name, rank)

    # Loop over each habitat and create a list of dataframes, one per habitat
    habitat_list <- lapply(unique(fish_year_subset$habitat), function(habitat){

      fish_year_subset %>%
        filter(habitat == !!habitat) %>%
        mutate(habitat = !!habitat)

    })

    # Calculate the total number of species per habitat
    total_species_per_habitat <- bind_rows(habitat_list) %>%
      group_by(habitat) %>%
      summarize(n_total = n()) %>%
      mutate(year = !!year)

    # Loop over the habitat list and pull out the species in each habitat
    unique_species <- lapply(habitat_list, function(df) unique(df$scientific_name))

    # Loop over the species list and determine which species are unique for a given habitat
    unique_species_list_year <- lapply(seq_along(habitat_list), function(i) {

      other_species_values <- unlist(unique_species[-i])
      unique_species_values <- setdiff(unlist(unique_species[i]),
                                       other_species_values)

      habitat_list[[i]] %>%
        filter(scientific_name %in% unique_species_values)

    })

    # Calculate the number of unique species per habitat
    unique_species_per_habitat <- bind_rows(unique_species_list_year) %>%
      group_by(habitat) %>%
      summarize(n_unique = n())

    full_join(total_species_per_habitat, unique_species_per_habitat, by = "habitat") %>%
      mutate(n_unique = replace_na(n_unique, 0),
             percent_unique = n_unique / n_total,
             year = !!year)

  })

  # Take the mean for each habitat of all years
  unique_species_df <- bind_rows(unique_species_list) %>%
    group_by(habitat) %>%
    summarize(n_unique = mean(n_unique),
              n_total = mean(n_total),
              percent_unique = mean(percent_unique))

  return(unique_species_df)

}
