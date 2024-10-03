library(shiny)
library(ggplot2)

# UI for the Shiny app
ui <- fluidPage(
  titlePanel("Decomposition Dynamics: The Role of Collembola and Microbes"),
  
  sidebarLayout(
    sidebarPanel(
      # User inputs for microarthropods and microbes
      sliderInput("collembola_density", "Collembola Density (Individuals per g soil):",
                  min = 0, max = 500, value = 100, step = 10),
      
      sliderInput("microbial_biomass", "Microbial Biomass (relative units):",
                  min = 0, max = 100, value = 50, step = 5),
      
      # Environmental conditions
      sliderInput("moisture", "Soil Moisture (0-100%):",
                  min = 0, max = 100, value = 50, step = 5),
      
      sliderInput("temperature", "Temperature (°C):",
                  min = 0, max = 40, value = 20, step = 5),
      
      # Button to run the simulation
      actionButton("run_simulation", "Run Simulation")
    ),
    
    mainPanel(
      plotOutput("decompositionPlot"),
      textOutput("feedback")
    )
  )
)

# Server logic for the Shiny app
server <- function(input, output, session) {
  # Function to calculate decomposition over time based on input parameters
  simulate_decomposition <- function(collembola, microbial_biomass, moisture, temperature, days = 365) {
    # Check for system crash
    if (collembola < 10 || collembola > 490) {
      return(rep(NA, days))  # System crash: return NA for all days
    }
    
    # Optimal Collembola density and its effect
    optimal_collembola <- 250
    reduction_per_100 <- 0.20  # 20% reduction per ±100 deviation from optimal
    
    # Calculate effect based on distance from optimal density
    deviation <- abs(collembola - optimal_collembola) / 100
    collembola_effect <- max(0, 1 - (reduction_per_100 * deviation))  # Ensure effect is not negative
    
    # Microbial influence
    microbe_effect <- log(microbial_biomass + 1)
    
    # Environmental influences
    moisture_effect <- 1 - abs(moisture - 50) / 50
    temp_effect <- exp(-0.1 * abs(temperature - 20))
    
    # Initialize decomposition variables
    decomposition <- numeric(days)
    decomposition[1] <- 0
    
    # Daily decomposition based on combined influences
    for (day in 2:days) {
      # Calculate a daily decomposition rate affected by all factors
      daily_rate <- collembola_effect * microbe_effect * moisture_effect * temp_effect * (1 - decomposition[day - 1] / 100)
      decomposition[day] <- decomposition[day - 1] + daily_rate
      decomposition[day] <- min(decomposition[day], 100)  # Cap at 100%
    }
    
    return(decomposition)
  }
  
  # Event to run the simulation and update the plot
  observeEvent(input$run_simulation, {
    # Simulate decomposition over a year
    decomposition <- simulate_decomposition(input$collembola_density, input$microbial_biomass, 
                                            input$moisture, input$temperature)
    
    days <- 1:length(decomposition)
    
    # Output the decomposition progress over time
    output$decompositionPlot <- renderPlot({
      if (all(is.na(decomposition))) {
        ggplot() + 
          annotate("text", x = 1, y = 1, label = "System Crashed: Unstable Collembola Density", size = 6, color = "red") +
          theme_void()
      } else {
        ggplot(data.frame(days, decomposition), aes(x = days, y = decomposition)) +
          geom_line(color = "darkgreen", size = 1) +
          labs(x = "Days", y = "Cumulative Thatch Decomposition (%)", 
               title = "Decomposition Over Time") +
          theme_minimal()
      }
    })
    
    # Feedback based on the simulation outcome
    output$feedback <- renderText({
      if (all(is.na(decomposition))) {
        "The system crashed due to an unstable Collembola density (too low or too high)."
      } else if (max(decomposition) >= 90) {
        "Great job! You've reached high decomposition levels!"
      } else if (max(decomposition) >= 60) {
        "Good progress! Can you optimize the conditions further?"
      } else {
        "Decomposition is low. Try adjusting Collembola density or environmental factors."
      }
    })
  })
}

# Run the Shiny app
shinyApp(ui, server)

