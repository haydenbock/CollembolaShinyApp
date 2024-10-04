# Load required libraries
library(shiny)
library(shinythemes)
library(plotly)
library(deSolve)  # For solving differential equations

# Define UI
ui <- fluidPage(
  theme = shinytheme("cerulean"),
  titlePanel("Soil Decomposition Dynamics: Interactions Between Collembola and Microbes"),
  
  sidebarLayout(
    sidebarPanel(
      h4("Initial Conditions"),
      sliderInput("init_collembola", "Initial Collembola Density (Individuals per g soil):",
                  min = 0, max = 500, value = 100, step = 10),
      sliderInput("init_microbes", "Initial Microbial Biomass (Relative Units):",
                  min = 0, max = 100, value = 50, step = 5),
      
      h4("Environmental Factors"),
      sliderInput("temperature", "Soil Temperature (°C):",
                  min = 0, max = 40, value = 20, step = 1),
      sliderInput("moisture", "Soil Moisture (% Water Holding Capacity):",
                  min = 0, max = 100, value = 50, step = 5),
      
      actionButton("run", "Run Simulation")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Decomposition",
                 plotlyOutput("decompositionPlot"),
                 p("This plot shows the cumulative decomposition over time as a percentage, influenced by microbial activity and Collembola interactions.")
        ),
        tabPanel("Population Dynamics",
                 plotlyOutput("populationPlot"),
                 p("This plot illustrates the population dynamics of Collembola and microbial biomass over time.")
        ),
        tabPanel("Model Explanation",
                 h4("Overview"),
                 p("This simulation models the interactions between Collembola (soil microarthropods) and microorganisms in the decomposition of organic matter."),
                 h4("Collembola-Microbe Interaction"),
                 p("Collembola feed on microorganisms, which can reduce microbial biomass. However, their activity can also stimulate microbial growth by fragmenting organic matter and enhancing nutrient availability."),
                 h4("Environmental Effects"),
                 p("Soil temperature and moisture significantly affect biological activity. Optimal conditions promote higher rates of microbial growth and decomposition.")
        )
      )
    )
  )
)

# Define server logic
server <- function(input, output) {
  
  # Observe event for running the simulation
  observeEvent(input$run, {
    # Parameters
    params <- list(
      # Collembola parameters
      r_C = 0.1,        # Intrinsic growth rate of Collembola
      K_C = 500,        # Carrying capacity of Collembola
      a = 0.0005,       # Predation rate coefficient on microbes
      
      # Microbial parameters
      r_M = 0.2,        # Intrinsic growth rate of microbes
      K_M = 100,        # Carrying capacity of microbes
      b = 0.0005,       # Consumption rate coefficient by Collembola
      
      # Decomposition parameters
      d_0 = 0.01,       # Base decomposition rate
      
      # Environmental parameters
      T = input$temperature,     # Temperature (°C)
      T_opt = 25,                # Optimal temperature (°C)
      k_T = 0.1,                 # Temperature sensitivity
      W = input$moisture,        # Moisture (%)
      W_opt = 50,                # Optimal moisture (%)
      k_W = 0.01                 # Moisture sensitivity
    )
    
    # Initial conditions
    state <- c(
      C = input$init_collembola,    # Initial Collembola density
      M = input$init_microbes,      # Initial microbial biomass
      S = 100                       # Initial decomposable material (%)
    )
    
    # Time sequence fixed at 90 days
    times <- seq(0, 90, by = 1)
    
    # Define the model equations
    model <- function(t, state, params) {
      with(as.list(c(state, params)), {
        # Environmental effect function
        f_TW <- exp(k_T * (T - T_opt)) * exp(-k_W * (W - W_opt)^2)
        
        # Collembola population change
        dC_dt <- r_C * C * (1 - C / K_C) - a * C * M
        
        # Microbial biomass change
        dM_dt <- r_M * M * (1 - M / K_M) - b * C * M
        
        # Decomposable material change (ensuring D does not exceed 100%)
        dS_dt <- -d_0 * M * f_TW * (S / 100)
        
        # Return the rate of change
        list(c(dC_dt, dM_dt, dS_dt))
      })
    }
    
    # Solve the differential equations
    out <- ode(y = state, times = times, func = model, parms = params)
    
    # Convert output to data frame
    out_df <- as.data.frame(out)
    
    # Calculate cumulative decomposition as percentage
    out_df$D <- 100 - out_df$S
    
    # Ensure decomposition does not exceed 100%
    out_df$D[out_df$D > 100] <- 100
    
    # Decomposition Plot
    output$decompositionPlot <- renderPlotly({
      plot_ly(out_df, x = ~time, y = ~D, type = 'scatter', mode = 'lines', name = 'Decomposition', line = list(color = 'green')) %>%
        layout(
          title = "Cumulative Decomposition Over Time",
          xaxis = list(title = "Time (Days)"),
          yaxis = list(title = "Decomposition (%)", range = c(0, 100))
        )
    })
    
    # Population Dynamics Plot
    output$populationPlot <- renderPlotly({
      plot_ly(out_df, x = ~time) %>%
        add_lines(y = ~C, name = 'Collembola Density', line = list(color = 'blue')) %>%
        add_lines(y = ~M, name = 'Microbial Biomass', line = list(color = 'red')) %>%
        layout(
          title = "Population Dynamics Over Time",
          xaxis = list(title = "Time (Days)"),
          yaxis = list(title = "Population (Relative Units)"),
          legend = list(x = 0.1, y = 0.9)
        )
    })
  })
}

# Run the application
shinyApp(ui = ui, server = server)

