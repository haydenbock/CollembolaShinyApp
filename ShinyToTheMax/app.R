# Load required libraries
library(shiny)
library(shinythemes)
library(shinyBS)
library(shinyjs)
library(plotly)

# Simulation function
simulate_decomposition <- function(init_collembola, init_microbial_biomass, moisture, temperature_f, add_clippings = FALSE, days = 360) {
  # Convert temperature from Fahrenheit to Celsius
  temperature_c <- (temperature_f - 32) * 5 / 9
  
  # Check for extreme temperature conditions to generate flat lines
  if (temperature_f < 40 || temperature_f > 95) {
    flat_line <- rep(0, days)
    return(list(
      decomposition = flat_line,
      microbial_biomass = rep(init_microbial_biomass, days),
      collembola_density = rep(init_collembola, days)
    ))
  }
  
  # Check for system crash due to extreme Collembola densities
  if (init_collembola < 10 || init_collembola > 490) {
    return(list(
      decomposition = rep(NA, days),
      microbial_biomass = rep(NA, days),
      collembola_density = rep(NA, days)
    ))
  }
  
  # Check for extreme moisture levels that virtually stop decomposition
  if (moisture < 5 || moisture > 90) {
    return(list(
      decomposition = rep(NA, days),
      microbial_biomass = rep(NA, days),
      collembola_density = rep(NA, days)
    ))
  }
  
  # Initialize variables
  decomposition <- numeric(days)
  microbial_biomass <- numeric(days)
  collembola_density <- numeric(days)
  decomposition[1] <- 0
  microbial_biomass[1] <- init_microbial_biomass
  collembola_density[1] <- init_collembola
  
  # Parameters for Collembola population dynamics
  r_col <- 0.1  # Intrinsic growth rate
  K_col <- 500  # Carrying capacity
  
  # Environmental influences
  # Adjusted moisture effect to fine-tune decomposition
  moisture_effect <- ifelse(moisture < 30 | moisture > 70, 0.5, 1) *
    (1 - abs(moisture - 50) / 50)
  
  # Adjusted temperature effect for finer control
  if (temperature_f >= 50 && temperature_f <= 80) {
    temp_effect <- (temperature_f - 50) / (80 - 50)
  } else {
    temp_effect <- 0.5  # Reduced effect outside optimal range
  }
  
  # Adjusted base decomposition rate for lower overall decomposition
  base_decomp_rate <- 0.1  # Reduced from 0.5 to 0.1
  
  # Base standard deviations for stochasticity
  base_sigma_col <- 10
  base_sigma_microbe <- 0.5
  
  # Maximum microbial biomass for normalization
  max_microbial_biomass <- 100  # Assuming 100 is the maximum microbial biomass
  
  # Calculate environmental stochasticity factors based on temperature and moisture
  # For temperature: stochasticity increases with deviation from optimal 66°F
  temp_stochasticity <- 1 + (abs(temperature_f - 66) / 30)
  
  # For moisture: stochasticity increases when moisture is very high or very low
  moisture_stochasticity <- 1 + (abs(moisture - 50) / 50)
  
  # Total environmental stochasticity factor
  environmental_stochasticity <- temp_stochasticity * moisture_stochasticity
  
  # Parameters for stochastic events
  base_crash_prob <- 0.02  # Base probability of crash
  max_crash_prob <- 0.1    # Maximum probability of crash
  
  # Flag for clippings effect
  clippings_effect <- ifelse(add_clippings, 0.5, 1)  # Reduce decomposition by 50% if clippings are added
  
  for (day in 2:days) {
    # Update microbial biomass-dependent stochasticity factor
    microbial_stochasticity_factor <- 1 + (microbial_biomass[day - 1] / max_microbial_biomass)
    
    # Adjust sigma_col based on microbial biomass
    sigma_col <- base_sigma_col * environmental_stochasticity * microbial_stochasticity_factor
    sigma_microbe <- base_sigma_microbe * environmental_stochasticity
    
    # Determine crash probability based on Collembola density
    crash_prob <- base_crash_prob + (collembola_density[day - 1] / K_col) * (max_crash_prob - base_crash_prob)
    crash_prob <- min(crash_prob, max_crash_prob)
    
    # Adjust crash probability based on environmental stochasticity
    crash_prob <- crash_prob * environmental_stochasticity
    crash_prob <- min(crash_prob, 0.5)  # Cap crash probability at 50%
    
    # Check if a crash occurs
    if (runif(1) < crash_prob) {
      # Collembola population crashes to a random low level
      collembola_density[day] <- collembola_density[day - 1] * runif(1, 0.1, 0.3)
      # Microbial biomass increases due to reduced grazing
      microbial_biomass[day] <- microbial_biomass[day - 1] + abs(rnorm(1, mean = 2, sd = sigma_microbe))
    } else {
      # Collembola population growth with stochasticity
      col_growth <- r_col * collembola_density[day - 1] * (1 - collembola_density[day - 1] / K_col)
      stochastic_component <- rnorm(1, mean = 0, sd = sigma_col)
      collembola_density[day] <- collembola_density[day - 1] + col_growth + stochastic_component
      # Microbial biomass growth with stochasticity
      fluctuation_intensity <- max(0.1, abs(collembola_density[day] - 250) / 250)
      daily_growth <- temp_effect * (1 - decomposition[day - 1] / 100)^2 *
        (1 - fluctuation_intensity) * exp(-0.05 * decomposition[day - 1])
      daily_fluctuation <- rnorm(1, mean = 0, sd = sigma_microbe)
      microbial_biomass[day] <- microbial_biomass[day - 1] + daily_growth + daily_fluctuation
    }
    
    # Ensure populations stay within bounds
    collembola_density[day] <- max(min(collembola_density[day], K_col), 0)
    microbial_biomass[day] <- max(min(microbial_biomass[day], max_microbial_biomass), 0)
    
    # Recalculate collembola_effect based on current density
    optimal_collembola <- 200  # Adjusted optimal Collembola density
    reduction_per_100 <- 0.25  # Increased reduction factor
    deviation <- abs(collembola_density[day] - optimal_collembola) / 100
    collembola_effect <- max(0, 1 - (reduction_per_100 * deviation))
    
    # Recalculate recalcitrance effect
    recalcitrance_effect <- exp(-0.05 * decomposition[day - 1])
    
    # Calculate effect of Collembola densities on microbial biomass
    if (collembola_density[day] > 300) {
      col_over_threshold <- collembola_density[day] - 300
      col_microbe_reduction <- max(0, 1 - (col_over_threshold / 200))
    } else {
      col_microbe_reduction <- 1
    }
    
    # Calculate decomposition rate
    microbe_effect <- log(microbial_biomass[day] + 1)
    daily_decomp_rate <- max(0, collembola_effect^2 * microbe_effect * moisture_effect *
                               recalcitrance_effect * temp_effect * base_decomp_rate * clippings_effect)
    decomposition[day] <- decomposition[day - 1] + daily_decomp_rate
    decomposition[day] <- min(decomposition[day], 100)
  }
  
  return(list(
    decomposition = decomposition,
    microbial_biomass = microbial_biomass,
    collembola_density = collembola_density
  ))
}

# Define UI
ui <- fluidPage(
  theme = shinytheme("cerulean"),
  useShinyjs(),
  titlePanel("Decomposition Dynamics: The Role of Collembola and Microbes"),
  
  sidebarLayout(
    sidebarPanel(
      # Tabs for input parameters
      tabsetPanel(
        tabPanel("Organism Parameters",
                 sliderInput("collembola_density", "Initial Collembola Density (Individuals per g soil):",
                             min = 50, max = 400, value = 100, step = 10),
                 bsTooltip("collembola_density", "Collembola are microarthropods that influence decomposition."),
                 sliderInput("microbial_biomass", "Initial Microbial Biomass (relative units):",
                             min = 10, max = 100, value = 50, step = 5),
                 bsTooltip("microbial_biomass", "Microbial biomass represents the population of decomposers.")
        ),
        tabPanel("Environmental Conditions",
                 sliderInput("moisture", "Soil Moisture (0-100%):",
                             min = 0, max = 100, value = 60, step = 5),
                 bsTooltip("moisture", "Soil moisture affects microbial activity and population dynamics."),
                 sliderInput("temperature", "Temperature (°F):",
                             min = 32, max = 104, value = 66, step = 2),
                 bsTooltip("temperature", "Temperature influences biological processes and stochasticity.")
        )
      ),
      # Add Clippings Button
      actionButton("add_clippings", "Add Clippings"),
      bsTooltip("add_clippings", "Add clippings to the system, reducing decomposition rate by 50%."),
      # Info box
      br(),
      bsCollapse(
        bsCollapsePanel("Learn More About Decomposition Dynamics",
                        "This simulation models the interactions between Collembola (microarthropods) and microbial biomass in soil ecosystems. Adjust the parameters to explore how different conditions affect decomposition rates.",
                        style = "info")
      )
    ),
    
    mainPanel(
      plotlyOutput("decompositionPlot"),
      plotlyOutput("microbialCollembolaPlot"),
      uiOutput("feedback"),
      uiOutput("scoreDisplay")
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  # Reactive value to store whether clippings have been added
  clippings_added <- reactiveVal(FALSE)
  
  # Observe Add Clippings button click
  observeEvent(input$add_clippings, {
    clippings_added(TRUE)
  })
  
  # Reactive simulation results
  simulation_results <- reactive({
    simulate_decomposition(
      input$collembola_density,
      input$microbial_biomass,
      input$moisture,
      input$temperature,
      add_clippings = clippings_added()
    )
  })
  
  # Reactive expressions for outputs
  decomposition <- reactive({
    simulation_results()$decomposition
  })
  
  microbial_biomass <- reactive({
    simulation_results()$microbial_biomass
  })
  
  collembola_density <- reactive({
    simulation_results()$collembola_density
  })
  
  days <- reactive({
    1:length(decomposition())
  })
  
  # Scoring system
  score <- reactive({
    if (all(is.na(decomposition()))) {
      return(0)
    } else {
      max_decomp <- max(decomposition(), na.rm = TRUE)
      # Simple scoring formula
      score_value <- round(max_decomp)
      return(score_value)
    }
  })
  
  # Display the score and feedback
  output$scoreDisplay <- renderUI({
    if (score() == 0) {
      tags$div(
        tags$h3("Your Score: 0"),
        tags$p("Try adjusting the parameters to improve decomposition.")
      )
    } else if (score() >= 90) {
      tags$div(
        tags$h3(paste("Your Score:", score())),
        tags$p("Excellent! You've mastered the decomposition dynamics!")
      )
    } else if (score() >= 60) {
      tags$div(
        tags$h3(paste("Your Score:", score())),
        tags$p("Good job! See if you can reach a higher score.")
      )
    } else {
      tags$div(
        tags$h3(paste("Your Score:", score())),
        tags$p("Keep trying! Adjust the conditions to improve your score.")
      )
    }
  })
  
  # Feedback messages
  output$feedback <- renderUI({
    if (all(is.na(decomposition()))) {
      tags$p("The system crashed due to unstable Collembola density or extreme moisture.", style = "color: red;")
    } else if (clippings_added()) {
      tags$p("Clippings added: Decomposition rate reduced by 50%.", style = "color: purple;")
    } else if (max(decomposition(), na.rm = TRUE) >= 90) {
      tags$p("Great job! You've reached high decomposition levels!", style = "color: green;")
    } else if (max(decomposition(), na.rm = TRUE) >= 60) {
      tags$p("Good progress! Can you optimize the conditions further?", style = "color: orange;")
    } else {
      tags$p("Decomposition is low. Try adjusting Collembola density, moisture, or temperature.", style = "color: blue;")
    }
  })
  
  # Decomposition plot
  output$decompositionPlot <- renderPlotly({
    decomposition_data <- decomposition()
    if (all(is.na(decomposition_data))) {
      plot_ly() %>%
        layout(
          title = "System Crashed: Unstable Conditions",
          annotations = list(
            x = 0.5, y = 0.5, text = "System Crashed: Unstable Conditions",
            showarrow = FALSE, font = list(size = 20, color = "red")
          ),
          xaxis = list(visible = FALSE),
          yaxis = list(visible = FALSE)
        )
    } else {
      plot_ly(data.frame(Days = days(), Decomposition = decomposition_data), x = ~Days, y = ~Decomposition, type = 'scatter', mode = 'lines',
              line = list(color = 'darkgreen')) %>%
        layout(
          title = "Decomposition Over Time",
          xaxis = list(title = "Days"),
          yaxis = list(title = "Cumulative Thatch Decomposition (%)")
        )
    }
  })
  
  # Microbial biomass and Collembola density plot
  output$microbialCollembolaPlot <- renderPlotly({
    microbial_data <- microbial_biomass()
    collembola_data <- collembola_density()
    if (all(is.na(microbial_data))) {
      plot_ly() %>%
        layout(
          title = "System Crashed: Unstable Conditions",
          annotations = list(
            x = 0.5, y = 0.5, text = "System Crashed: Unstable Conditions",
            showarrow = FALSE, font = list(size = 20, color = "red")
          ),
          xaxis = list(visible = FALSE),
          yaxis = list(visible = FALSE)
        )
    } else {
      df <- data.frame(Days = days(), Microbial_Biomass = microbial_data, Collembola_Density = collembola_data)
      plot_ly(df, x = ~Days) %>%
        add_lines(y = ~Microbial_Biomass, name = 'Microbial Biomass', line = list(color = 'blue')) %>%
        add_lines(y = ~Collembola_Density, name = 'Collembola Density', line = list(color = 'red')) %>%
        layout(
          title = "Microbial Biomass and Collembola Density Dynamics",
          xaxis = list(title = "Days"),
          yaxis = list(title = "Population (relative units)"),
          legend = list(x = 0.1, y = 0.9)
        )
    }
  })
  
  # Reset clippings when parameters change
  observeEvent(list(input$collembola_density, input$microbial_biomass, input$moisture, input$temperature), {
    clippings_added(FALSE)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
