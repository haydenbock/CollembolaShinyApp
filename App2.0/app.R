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
      
      sliderInput("microbial_biomass", "Initial Microbial Biomass (relative units):",
                  min = 0, max = 100, value = 50, step = 5),
      
      # Environmental conditions
      sliderInput("moisture", "Soil Moisture (0-100%):",
                  min = 0, max = 100, value = 50, step = 5),
      
      sliderInput("temperature", "Temperature (°F):",
                  min = 32, max = 104, value = 68, step = 5),
      
      # Button to run the simulation
      actionButton("run_simulation", "Run Simulation")
    ),
    
    mainPanel(
      plotOutput("decompositionPlot"),
      plotOutput("microbialPlot"),
      textOutput("feedback")
    )
  )
)

# Server logic for the Shiny app
server <- function(input, output, session) {
  # Function to calculate decomposition and microbial biomass over time based on input parameters
  simulate_decomposition <- function(collembola, initial_microbial_biomass, moisture, temperature_f, days = 360) {
    # Convert temperature from Fahrenheit to Celsius
    temperature_c <- (temperature_f - 32) * 5 / 9
    
    # Check for extreme temperature conditions to generate flat lines
    if (temperature_f < 40 || temperature_f > 95) {
      flat_line <- rep(0, days)
      return(list(decomposition = flat_line, microbial_biomass = rep(initial_microbial_biomass, days)))  # Flat lines for both
    }
    
    # Check for system crash due to extreme Collembola densities
    if (collembola < 10 || collembola > 490) {
      return(list(decomposition = rep(NA, days), microbial_biomass = rep(NA, days)))  # System crash: return NA for all days
    }
    
    # Check for extreme moisture levels that virtually stop decomposition
    if (moisture < 5 || moisture > 90) {
      return(list(decomposition = rep(NA, days), microbial_biomass = rep(NA, days)))  # No decomposition if extreme moisture
    }
    
    # Optimal Collembola density and its effect
    optimal_collembola <- 250
    reduction_per_100 <- 0.20  # 20% reduction per ±100 deviation from optimal
    
    # Calculate effect based on distance from optimal density
    deviation <- abs(collembola - optimal_collembola) / 100
    collembola_effect <- max(0, 1 - (reduction_per_100 * deviation))  # Ensure effect is not negative
    
    # Initialize decomposition and microbial biomass variables
    decomposition <- numeric(days)
    microbial_biomass <- numeric(days)
    decomposition[1] <- 0
    microbial_biomass[1] <- initial_microbial_biomass
    
    # Environmental influences
    moisture_effect <- ifelse(moisture < 20 | moisture > 80, 0.5, 1) * (1 - abs(moisture - 50) / 50)  # Reduce efficacy if moisture is extreme
    
    # Temperature effect that gradually increases decomposition from 40°F to 80°F
    if (temperature_f >= 40 && temperature_f <= 80) {
      temp_effect <- (temperature_f - 40) / (80 - 40)  # Gradually increases from 0 at 40°F to 1 at 80°F
    } else {
      temp_effect <- 1  # Normal conditions
    }
    
    # Calculate the effect of high Collembola densities on microbial populations
    if (collembola > 300) {
      col_over_threshold <- collembola - 300
      col_microbe_reduction <- max(0, 1 - (col_over_threshold / 150))  # Gradually decreases microbial biomass growth, 450 Collembola = 0 effect
    } else {
      col_microbe_reduction <- 1  # No reduction if below the threshold
    }
    
    # Scale the rate so that optimal decomposition takes around 360 days
    base_decomp_rate <- 0.3  # Base rate of decomposition per day under optimal conditions
    
    # Daily dynamics based on combined influences
    for (day in 2:days) {
      # Recalculate the effect of decreasing resources on microbial growth
      recalcitrance_effect <- exp(-0.01 * decomposition[day - 1])  # Exponential decay as thatch is decomposed
      
      # Calculate daily decomposition rate
      microbe_effect <- log(microbial_biomass[day - 1] + 1)
      daily_decomp_rate <- collembola_effect * microbe_effect * moisture_effect * recalcitrance_effect * temp_effect * base_decomp_rate
      decomposition[day] <- decomposition[day - 1] + daily_decomp_rate
      decomposition[day] <- min(decomposition[day], 100)  # Cap at 100%
      
      # Calculate daily change in microbial biomass
      fluctuation_intensity <- max(0.1, abs(collembola - optimal_collembola) / 250)  # Control fluctuation intensity
      daily_growth <- temp_effect * (1 - decomposition[day] / 100) * (1 - fluctuation_intensity) * recalcitrance_effect * col_microbe_reduction  # Resource & temp-dependent growth with Collembola effect
      daily_fluctuation <- rnorm(1, mean = 0, sd = fluctuation_intensity * 0.5)  # Random fluctuation
      
      # Update microbial biomass with a mix of growth and stochastic fluctuation
      microbial_biomass[day] <- microbial_biomass[day - 1] + daily_growth + daily_fluctuation
      microbial_biomass[day] <- max(microbial_biomass[day], 0)  # Biomass cannot be negative
    }
    
    return(list(decomposition = decomposition, microbial_biomass = microbial_biomass))
  }
  
  # Event to run the simulation and update the plot
  observeEvent(input$run_simulation, {
    # Simulate decomposition and microbial biomass over a year
    results <- simulate_decomposition(input$collembola_density, input$microbial_biomass, 
                                      input$moisture, input$temperature)
    decomposition <- results$decomposition
    microbial_biomass <- results$microbial_biomass
    days <- 1:length(decomposition)
    
    # Output the decomposition progress over time
    output$decompositionPlot <- renderPlot({
      if (all(is.na(decomposition))) {
        ggplot() + 
          annotate("text", x = 1, y = 1, label = "System Crashed: Unstable Conditions", size = 6, color = "red") +
          theme_void()
      } else {
        ggplot(data.frame(days, decomposition), aes(x = days, y = decomposition)) +
          geom_line(color = "darkgreen", size = 1) +
          labs(x = "Days", y = "Cumulative Thatch Decomposition (%)", 
               title = "Decomposition Over Time") +
          theme_minimal()
      }
    })
    
    # Output the microbial biomass dynamics over time
    output$microbialPlot <- renderPlot({
      if (all(is.na(microbial_biomass))) {
        ggplot() + 
          annotate("text", x = 1, y = 1, label = "System Crashed: Unstable Conditions", size = 6, color = "red") +
          theme_void()
      } else {
        ggplot(data.frame(days, microbial_biomass), aes(x = days, y = microbial_biomass)) +
          geom_line(color = "blue", size = 1) +
          labs(x = "Days", y = "Microbial Biomass (relative units)", 
               title = "Microbial Biomass Dynamics Over Time") +
          theme_minimal()
      }
    })
    
    # Feedback based on the simulation outcome
    output$feedback <- renderText({
      if (all(is.na(decomposition))) {
        "The system crashed due to unstable Collembola density or extreme moisture."
      } else if (max(decomposition) >= 90) {
        "Great job! You've reached high decomposition levels!"
      } else if (max(decomposition) >= 60) {
        "Good progress! Can you optimize the conditions further?"
      } else {
        "Decomposition is low. Try adjusting Collembola density, moisture, or temperature."
      }
    })
  })
}

# Run the Shiny app
shinyApp(ui, server)

