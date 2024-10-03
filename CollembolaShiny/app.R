library(shiny)

# UI for the Shiny app
ui <- fluidPage(
  titlePanel("Thatch Decomposition and Soil Dynamics Toy Model"),
  
  sidebarLayout(
    sidebarPanel(
      # User inputs to control the clipping treatment
      selectInput("clipping", "Clipping Treatment:",
                  choices = c("Coarse Frequent", "Coarse Infrequent", 
                              "Fine Frequent", "Fine Infrequent", "No Clippings")),
      
      # Slider for adjusting Collembola density
      sliderInput("collembola_density", "Collembola Density (Individuals per g soil):",
                  min = 0, max = 500, value = 50, step = 10),
      
      # Slider for simulation time (in days)
      sliderInput("time", "Simulation Time (Days):",
                  min = 0, max = 180, value = 90, step = 10),
      
      # Action button to pause and change conditions
      actionButton("pause", "Pause and Change Conditions")
    ),
    
    # Main panel to display results
    mainPanel(
      tabsetPanel(
        tabPanel("Thatch Decomposition Over Time", plotOutput("decompositionPlot")),
        tabPanel("Microbial Biomass Over Time", plotOutput("microbialPlot")),
        tabPanel("Extracellular Enzymes Over Time", plotOutput("enzymePlot")),
        tabPanel("CO2 Respiration Over Time", plotOutput("respirationPlot"))
      )
    )
  )
)

# Server logic for the Shiny app
server <- function(input, output, session) {
  # Reactive values to store the state
  rv <- reactiveValues(
    paused_day = NULL,  # Day at which the simulation was paused
    initial_decomp = NULL,  # Initial decomposition trajectory
    initial_time = NULL  # Initial time sequence
  )
  
  # Parameters for polynomial models (example coefficients)
  coefficients <- reactive({
    coeffs <- list()
    
    # Coefficients for different relationships based on your findings
    coeffs$decomposition <- list(
      "Coarse Frequent" = c(0.5, 0.02, -0.00003),
      "Coarse Infrequent" = c(0.7, 0.015, -0.00002),
      "Fine Frequent" = c(0.6, 0.018, -0.000025),
      "Fine Infrequent" = c(0.65, 0.017, -0.000024),
      "No Clippings" = c(0.3, 0.01, -0.000015)
    )
    
    coeffs$microbial_biomass <- c(100, 0.8, -0.002)
    coeffs$enzyme_production <- c(50, 0.05, -0.0001)
    coeffs$CO2_respiration <- c(200, 1.5, -0.005)
    
    return(coeffs)
  })
  
  # Polynomial model function
  polynomial_model <- function(coeffs, x) {
    coeffs[1] + coeffs[2] * x + coeffs[3] * x^2
  }
  
  # Reactive expression to calculate the initial simulation
  initial_simulation <- reactive({
    coeffs <- coefficients()
    col_density <- input$collembola_density
    decomp_coeffs <- coeffs$decomposition[[input$clipping]]
    time <- seq(0, input$time, by = 1)
    
    # Thatch decomposition over time
    decomp_rate <- polynomial_model(decomp_coeffs, col_density)
    cumulative_decomp <- decomp_rate * (1 - exp(-0.01 * time))
    
    list(time = time, decomp = cumulative_decomp)
  })
  
  # Event to pause the simulation and allow condition changes
  observeEvent(input$pause, {
    # Store the current paused day and initial simulation results
    rv$paused_day <- input$time
    rv$initial_time <- initial_simulation()$time
    rv$initial_decomp <- initial_simulation()$decomp
  })
  
  # Calculate decomposition over time with updated conditions
  output$decompositionPlot <- renderPlot({
    coeffs <- coefficients()
    col_density <- input$collembola_density
    decomp_coeffs <- coeffs$decomposition[[input$clipping]]
    time <- seq(0, input$time, by = 1)
    
    # Thatch decomposition over time
    decomp_rate <- polynomial_model(decomp_coeffs, col_density)
    cumulative_decomp <- decomp_rate * (1 - exp(-0.01 * time))
    
    # If paused, combine initial and updated trajectories
    if (!is.null(rv$paused_day) && input$time > rv$paused_day) {
      new_time <- seq(rv$paused_day, input$time, by = 1)
      updated_decomp <- decomp_rate * (1 - exp(-0.01 * (new_time - rv$paused_day)))
      final_decomp <- c(rv$initial_decomp, rv$initial_decomp[length(rv$initial_decomp)] + updated_decomp)
      final_time <- c(rv$initial_time, new_time)
      
      # Plot both initial and updated trajectories
      plot(rv$initial_time, rv$initial_decomp, type = "l", col = "darkgreen", lwd = 2,
           xlab = "Time (Days)", ylab = "Cumulative Thatch Mass Loss (%)",
           main = "Thatch Decomposition Over Time")
      lines(final_time, final_decomp, col = "blue", lwd = 2, lty = 2)
      legend("topright", legend = c("Initial Conditions", "Updated Conditions"), col = c("darkgreen", "blue"), lwd = 2, lty = c(1, 2))
    } else {
      plot(time, cumulative_decomp, type = "l", col = "darkgreen", lwd = 2,
           xlab = "Time (Days)", ylab = "Cumulative Thatch Mass Loss (%)",
           main = "Thatch Decomposition Over Time")
    }
    
    grid()
  })
  
  # Placeholder for other plots (for simplicity, only Thatch Decomposition is fully implemented)
  output$microbialPlot <- renderPlot({
    plot(1, 1, main = "Microbial Biomass Over Time (not yet implemented)")
  })
  output$enzymePlot <- renderPlot({
    plot(1, 1, main = "Extracellular Enzymes Over Time (not yet implemented)")
  })
  output$respirationPlot <- renderPlot({
    plot(1, 1, main = "CO2 Respiration Over Time (not yet implemented)")
  })
}

# Run the Shiny app
shinyApp(ui, server)
