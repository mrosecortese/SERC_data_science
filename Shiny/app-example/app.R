# To start the application, click the "Run App" in the script toolbar.
# By default, it will start in an RStudio window, but it is usually
# preferred to run in a browser window. Click the little white arrow facing
# down next to "Run App" and select "Run External".

# Our application has a `/R` folder. Shiny applications automatically source all
# scripts they find in the `/R` folder, making it a convenient place to store
# custom functions.

# Place package loading calls and
# static variables and data at the top of the app.R script.
# In a 3 file application, this top section would be in `global.R`.

library(shiny)
library(tidyverse)
library(bslib)
library(leaflet)

fish_data <- read_csv("resources/cbc_rls_shiny_workshop_data.csv")
coordinates <- read_csv("resources/RLS Metadata.csv")

# The UI function. In a 3 file application, this would go in the `ui.R` script.
ui <- page_sidebar(

  title = "Blinded by the Bright Dashboard",

  sidebar = sidebar(

    selectInput(inputId = "habitat",
                label = "Select habitats to include",
                choices = unique(fish_data$habitat),
                selected = unique(fish_data$habitat),
                multiple = TRUE),

    selectInput(inputId = "year",
                label = "Select years to include",
                choices = unique(fish_data$year),
                selected = unique(fish_data$year),
                multiple = TRUE),
  ),

  layout_columns(
    # A negative number in widths creates white space
    col_widths = c(-2,8,
                   6,6,
                   6,6),

    card(leafletOutput(outputId = "map")),
    card(plotOutput(outputId = "abundance_plot")),
    card(plotOutput(outputId = "richness_plot")),
    card(plotOutput(outputId = "unique_sp_plot")),
    card(plotOutput(outputId = "percent_unique_sp_plot"))
  )

)

# The server function. In a 3 file app, this would go in server.R script
server <- function(input, output, session) {

  # This reactive function will be called by other reactives and the plots, as necessary
  fish_subset <- reactive({

    fish_subset <- fish_data %>%
      filter(habitat %in% input$habitat,
             year %in% input$year)

    return(fish_subset)

  })

  # This reactive calls the function we placed in /R
  habitat_summary <- reactive({

    fish_subset <- fish_subset()
    habitat_summary <- summarize_by_habitat_year(fish_subset)

    return(habitat_summary)

  })

  output$abundance_plot <- renderPlot({

    fish_subset() %>%
      filter(habitat %in% input$habitat,
             year %in% input$year) %>%
      group_by(habitat, year, location_name) %>%
      summarize(abundance = sum(size_count, na.rm = T)) %>%
      mutate(log10_abundance = log10(abundance)) %>%
      ggplot(aes(x = log10_abundance, fill = habitat)) + geom_histogram() +
      ggtitle("Fish Abundance per Dive") +
      theme_minimal() + labs(y = "Count of Dives") +
      theme(text = element_text(size = 18))

  })

  output$richness_plot <- renderPlot({

    habitat_summary()  %>%
      ggplot(aes(x = habitat, y = n_total, fill = habitat)) +
      geom_bar(stat = "identity") +
      theme_minimal() + labs(y = "Species Richness") +
      theme(text = element_text(size = 18))

  })

  output$unique_sp_plot <- renderPlot({

    habitat_summary()  %>%
      ggplot(aes(x = habitat, y = n_unique, fill = habitat)) +
      geom_bar(stat = "identity") +
      theme_minimal() + labs(y = "Total Unique Species") +
      theme(text = element_text(size = 18))


  })

  output$percent_unique_sp_plot <- renderPlot({

    habitat_summary()  %>%
      ggplot(aes(x = habitat, y = percent_unique, fill = habitat)) +
      geom_bar(stat = "identity") +
      theme_minimal() + labs(y = "Percent Unique Species") +
      theme(text = element_text(size = 18))

  })

  output$map <- renderLeaflet({

    locations_to_map <- fish_subset() %>%
      count(location_name)

    coordinates %>%
      right_join(locations_to_map, by = "location_name") %>%
      leaflet() %>%
      addTiles() %>%
      addMarkers(lat = ~transect_decimal_latitude, lng = ~transect_decimal_longitude)

  })
}

shinyApp(ui, server)
