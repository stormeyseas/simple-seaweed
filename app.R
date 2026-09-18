library(shiny)
library(bslib)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

# Files in R/ are sourced automatically by Shiny, in alphabetical order
# Pre-calculate optimum light and velocity once, at startup.
species <- precalc_species(species_list)

# Make ggplot2 output match the app theme
thematic::thematic_shiny()

# Shared plot formatting; add to any ggplot with `+ plot_format`
plot_format <- list(
  theme_classic(),
  theme(
    legend.position = "right",
    text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14)
  )
)

# UI ------------------------------------------------------------------------

ui <- page_sidebar(
  title = "Simple seaweed growth simulator",
  theme = bs_theme(version = 5),
  sidebar = sidebar(
    title = "Environmental conditions",
    sliderInput("temperature", "Temperature (°C)", min = 5, max = 35, value = 18, step = 0.5),
    sliderInput("salinity", "Salinity (PSU)", min = 20, max = 45, value = 35, step = 0.5),
    sliderInput("nitrate", "Nitrate (µmol/L)", min = 0, max = 30, value = 10, step = 0.5),
    sliderInput("ammonium", "Ammonium (µmol/L)", min = 0, max = 30, value = 5, step = 0.5),
    helpText("Conditions are held constant over each species' growing period. Results update as you move the sliders.")
  ),
  uiOutput("value_boxes"),
  navset_card_underline(
    full_screen = TRUE,
    nav_panel("Biomass", plotOutput("biomass_plot")),
    nav_panel("Nutrients", plotOutput("nutrient_plot")),
    nav_panel(
      "Fixed settings",
      tableOutput("settings"),
      card(
        card_header("Notes"),
        shiny::tags$ul(
          shiny::tags$li("All species started with a biomass of 3 g/L."),
          shiny::tags$li("This simulation keeps farm parameters fixed: farm area = 2500m2, water depth = 30m, canopy width = 3m, top of canopy = 1m from the surface.")
        )
      )
    )
  )
)

# Server --------------------------------------------------------------------

server <- function(input, output, session) {
  # Debounce so dragging a slider doesn't trigger a model run on every pixel
  conditions <- reactive({
    list(
      temperature = input$temperature,
      salinity = input$salinity,
      nitrate = input$nitrate,
      ammonium = input$ammonium
    )
  }) |>
    debounce(300)

  # Runs at launch and whenever any slider changes
  results <- reactive({
    env <- conditions()
    lapply(names(species), function(id) {
      sp <- species[[id]]
      run_growth(
        sp,
        temperature = env$temperature,
        salinity = env$salinity,
        light = sp$light,
        velocity = sp$velocity,
        nitrate = env$nitrate,
        ammonium = env$ammonium
      ) |>
        mutate(species = sp$label, .before = 1)
    }) |>
      bind_rows() |>
      as_tibble()
  })

  output$value_boxes <- renderUI({
    req(results())
    boxes <- results() |>
      group_by(species) |>
      summarise(
        biomass_gL = mgm3_to_gL(last(B_ww.mg)),
        # Nf and Ns are mg N m-3; report as g N m-3
        N_removed = mg_to_g(last(Nf) + last(Ns) - first(Nf) - first(Ns)),
        .groups = "drop"
      )

    layout_column_wrap(
      width = "200px", fill = FALSE,
      !!!lapply(seq_len(nrow(boxes)), function(i) {
        value_box(
          title = boxes$species[i],
          value = span(
            style = "font-size: 150%;",
            sprintf("%.1f g WW/L", boxes$biomass_gL[i])
          ),
          p(sprintf("N removed: %.1f g N/m\u00b3", boxes$N_removed[i])),
          theme = if (boxes$N_removed[i] >= 0) "success" else "danger"
        )
      })
    )
  })

  output$biomass_plot <- renderPlot({
    req(results())
    ggplot(results(), aes(x = t, y = mgm3_to_gL(B_ww.mg), colour = species)) +
      geom_line(linewidth = 1) +
      scale_y_continuous(breaks = seq(0, 60, 5), limits = c(4.5, 35)) +
      scale_colour_manual(values = species_pal) +
      plot_format +
      labs(
        x = "Day", 
        y = "Wet weight biomass (g/L)", 
        colour = "Species"
      )
  })

  output$nutrient_plot <- renderPlot({
    req(results())
    results() |>
      select(species, t, up_Ni, up_Am) |>
      pivot_longer(-c(species, t), names_to = "form") |>
      mutate(form = factor(form, levels = c("up_Am", "up_Ni"), labels = c("Ammonium", "Nitrate"))) |>
      ggplot(aes(x = t, y = mg_to_g(value), colour = form)) +
      geom_line(linewidth = 1) +
      scale_y_continuous(breaks = seq(0, 2.5, 0.5), limits = c(0,2.5)) +
      scale_colour_brewer(palette = "Set1") +
      plot_format +
      facet_wrap(~species) +
      labs(x = "Day", y = "N uptake (g/d)", colour = NULL)
  })

  output$settings <- renderTable(settings_table(species))
}

shinyApp(ui, server)
