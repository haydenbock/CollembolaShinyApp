# Load necessary libraries
library(shiny)
library(plotly)

# Define the UI for the app
ui <- fluidPage(
  
  # App title
  titlePanel("Decomposer Dash: Interactive Thatch Decomposition Simulator"),
  
  # Sidebar layout with input controls
  sidebarPanel(
    # Input for Collembola presence
    checkboxInput("collembola", "Add Collembola", value = FALSE),
    
    # Input for clipping frequency
    selectInput("frequency", "Clipping Addition Frequency:",
                choices = c("Frequent", "Infrequent")),
    
    # Input for clipping size
    selectInput("size", "Clipping Size:",
                choices = c("No Clippings", "Coarse", "Fine")),
    
    # Day progression button
    actionButton("nextDay", "Progress Day"),
    
    # Display current day
    h3("Day Progress:"),
    textOutput("dayCount")
  ),
  
  # Main panel to show the decomposition plot
  mainPanel(
    plotlyOutput("decompPlot")
  )
)

# Define server logic for the simulation
server <- function(input, output, session) {
  
  # Reactive values to track simulation state
  rv <- reactiveValues(
    day = 1,
    cumulative_decomposition = numeric(180)
  )
  
  # Update day count
  output$dayCount <- renderText({
    paste("Day", rv$day)
  })
  
  # Reactive simulation parameters based on inputs
  sim_params <- reactive({
    collembola_effect <- if (input$collembola) 1.35 else 0.9
    
    # Determine decomposition rate based on clipping size and frequency
    decomposition_rate <- switch(input$size,
                                 "No Clippings" = 0.3,  # Greatest decomposition
                                 "Coarse" = if (input$frequency == "Infrequent") 0.2 else 0.15,  # Coarse-Infrequent > Coarse-Frequent
                                 "Fine" = if (input$frequency == "Infrequent") 0.08 else 0.03  # Fine-Infrequent > Fine-Frequent
    )
    
    list(collembola_effect = collembola_effect,
         decomposition_rate = decomposition_rate)
  })
  
  # Observe the day progression
  observeEvent(input$nextDay, {
    params <- sim_params()
    
    # Simulate decomposition for the current day
    new_decomposition <- params$decomposition_rate * params$collembola_effect
    
    # Update cumulative totals for decomposition
    rv$day <- rv$day + 1
    if (rv$day > 1) {
      rv$cumulative_decomposition[rv$day] <- rv$cumulative_decomposition[rv$day - 1] + new_decomposition
    } else {
      rv$cumulative_decomposition[rv$day] <- new_decomposition
    }
  })
  
  # Reactive line plot for decomposition over time
  output$decompPlot <- renderPlotly({
    if (rv$day > 1) {
      plot_ly(x = 1:rv$day, y = rv$cumulative_decomposition[1:rv$day], type = 'scatter', mode = 'lines',
              line = list(color = 'blue', width = 2)) %>%
        layout(title = "Thatch Decomposition Over Time",
               xaxis = list(title = "Days"),
               yaxis = list(title = "Cumulative Decomposition (%)"))
    } else {
      plot_ly(x = numeric(0), y = numeric(0), type = 'scatter', mode = 'lines',
              line = list(color = 'blue', width = 2)) %>%
        layout(title = "Thatch Decomposition Over Time",
               xaxis = list(title = "Days"),
               yaxis = list(title = "Cumulative Decomposition (%)"))
    }
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

