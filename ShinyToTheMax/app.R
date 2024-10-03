# Load required libraries
library(shiny)
library(shinythemes)
library(shinyBS)
library(shinyjs)
library(plotly)

# Simulation function
simulate_decomposition <- function(init_collembola, init_microbial_biomass,
                                   moisture, temperature_f, add_clippings = FALSE,
                                   days = 360) {
  # Convert temperature from Fahrenheit to Celsius (unused in this model)
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
  if (init_collembola < 10 || init_collembola > 500) {
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
  # Adjusted moisture effect
  moisture_effect <- ifelse(moisture < 30 | moisture > 70, 0.5, 1) *
    (1 - abs(moisture - 50) / 50)
  
  # Adjusted temperature effect, peaks at 1 when temperature is 66°F
  temp_effect <- max(1 - abs(temperature_f - 66) / 30, 0)
  
  # Adjusted base decomposition rate
  base_decomp_rate <- 0.05  # Adjusted to achieve desired decomposition percentages
  
  # Base standard deviations for stochasticity
  base_sigma_col <- 10
  base_sigma_microbe <- 0.5
  
  # Maximum microbial biomass
  max_microbial_biomass <- 100
  
  # Environmental stochasticity factors
  temp_stochasticity <- 1 + (abs(temperature_f - 66) / 30)
  moisture_stochasticity <- 1 + (abs(moisture - 50) / 50)
  environmental_stochasticity <- temp_stochasticity * moisture_stochasticity
  
  # Parameters for stochastic events
  base_crash_prob <- 0.02
  max_crash_prob <- 0.1
  
  # Flag for clippings effect
  clippings_effect <- ifelse(add_clippings, 0.5, 1)
  
  # Parameters for Collembola effect on decomposition
  optimal_collembola <- 250  # Optimal density for maximum decomposition
  width <- 50                # Controls the width of the peak in the exponential function
  
  # Parameters for microbial biomass growth
  r_max <- 0.1    # Maximum microbial growth rate per day
  r_min <- -0.15  # Minimum microbial growth rate per day
  n_microbe <- 3  # Exponent to exaggerate the effect of Collembola density on microbial growth
  
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
      
      # Microbial biomass growth rate influenced by Collembola density
      # Exaggerated inverse relationship
      microbial_growth_rate <- r_max - (r_max - r_min) * (collembola_density[day] / K_col)^n_microbe
      microbial_growth_rate <- max(microbial_growth_rate, r_min)
      
      # Calculate daily growth
      daily_growth <- microbial_growth_rate * microbial_biomass[day - 1] * temp_effect * exp(-0.05 * decomposition[day - 1])
      # Add stochastic fluctuation
      daily_fluctuation <- rnorm(1, mean = 0, sd = sigma_microbe)
      # Update microbial biomass
      microbial_biomass[day] <- microbial_biomass[day - 1] + daily_growth + daily_fluctuation
    }
    
    # Ensure populations stay within bounds
    collembola_density[day] <- max(min(collembola_density[day], K_col), 0)
    microbial_biomass[day] <- max(min(microbial_biomass[day], max_microbial_biomass), 0)
    
    # Calculate collembola_effect on decomposition
    collembola_effect <- exp(-((collembola_density[day] - optimal_collembola) / width)^2)
    # Exaggerate the effect
    collembola_effect <- collembola_effect^3  # Cubed to exaggerate
    
    # Recalculate recalcitrance effect
    recalcitrance_effect <- exp(-0.05 * decomposition[day - 1])
    
    # Calculate decomposition rate
    daily_decomp_rate <- base_decomp_rate * microbial_biomass[day] * collembola_effect *
      moisture_effect * temp_effect * recalcitrance_effect * clippings_effect
    
    decomposition[day] <- decomposition[day - 1] + daily_decomp_rate
    decomposition[day] <- min(decomposition[day], 100)
  }
  
  return(list(
    decomposition = decomposition,
    microbial_biomass = microbial_biomass,
    collembola_density = collembola_density
  ))
}

# The rest of the code (UI and server) remains the same as in the previous version.
