# Load required packages
library(burnout)  # Assuming this package contains P1_prob and final_size functions

# Function to simulate metapopulation dynamics with infection
simulate_metapopulation <- function(
    n_patches = 100,           # Number of patches
    n_years = 100,             # Number of years to simulate
    K = 1e6,                 # Carrying capacity per patch
    r = 0.5,                  # Growth rate
    c = 0.2,                  # Migration probability
    alpha = 5e-5,              # Scaling constant for infection probability
    D = 1,                  # Disease-induced mortality
    R0 = 2.5,                   # basic reproduction number
    epsilon=0.02,             # infectious period 
    initial_infected = 10,     # Number of initially infected patches
    initial_pop_ratio = 1   # Initial population as ratio of carrying capacity
) {
  
  # Initialize arrays to store population size and infection status
  # Dimensions: [patch, year, season]
  # Season: 1=beginning, 2=after growth, 3=after colonization, 4=after fizzle, 5=after burnout
  N <- array(0, dim = c(n_patches, n_years, 5))
  I <- array(0, dim = c(n_patches, n_years, 5))
  
  # Initialize an array to store total deaths for each year
  total_inf <- numeric(n_years)
  
  # Set initial conditions
  N[, 1, 1] <- initial_pop_ratio * K  # Initial population in all patches
  
  # Randomly select initial infected patches
  infected_patches <- sample(1:n_patches, initial_infected)
  I[infected_patches, 1, 1] <- 1
  
  # Simulate over years
  for (k in 1:n_years) {
    
    cat("Year:", k, "Season: 1", "Infected Patches:", sum(I[, k, 1]), "\n")
    
    # Within each year, apply the four processes
    
    # 1. Growth process
    for (i in 1:n_patches) {
      # Logistic growth
      N[i, k, 2] <- N[i, k, 1] + r * N[i, k, 1] * (1 - N[i, k, 1] / K)
      I[i, k, 2] <- I[i, k, 1]  # Infection status doesn't change during growth
    }
    cat("Year:", k, "Season: 2", "Infected Patches:", sum(I[, k, 2]), "\n")
    
    # 2. Colonization process
    # 2.1 Population movement
    total_pop <- sum(N[, k, 2])
    for (i in 1:n_patches) {
      N[i, k, 3] <- (1 - c) * N[i, k, 2] + c * total_pop / n_patches
    }
    
    # 2.2 Infection spread
    if (k == 1) {
      # First year: infection status remains unchanged during colonization
      I[, k, 3] <- I[, k, 2]
    } else {
      # Calculate infection probability based on previous year's deaths
      P_I <-1-exp(-alpha * c * total_inf[k - 1] / n_patches)  
      cat("  Year:", k, "  Season: 2", "  PI:", P_I, "\n")
      
      # Apply infection to susceptible patches
      for (i in 1:n_patches) {
        if (I[i, k, 2] == 0) {  # If patch is susceptible
          I[i, k, 3] <- ifelse(runif(1) < P_I, 1, 0)
        } else {
          I[i, k, 3] <- 1  # Already infected patches remain infected
        }
        
        #cat("Patch:", i,"  Year:", k, "  Season: 2", "  Population:", N[i, k, 3], "  Status:",I[i,k,3],"  Status change:",I[i, k, 3]-I[i,k,2], "\n")
      }
    }
    cat("Year:", k, "Season: 3", "Infected Patches:", sum(I[, k, 3]), "\n")
    
    # 3. Fizzle process
    P_f <- 1 / R0  # Fizzle probability
    cat("  Year:", k, "  Season: 3", "  Pf:", P_f, "\n")
    for (i in 1:n_patches) {
      # Infection may fizzle out
      if (I[i, k, 3] == 1) {
        I[i, k, 4] <- ifelse(runif(1) > P_f, 1, 0)
      } else {
        I[i, k, 4] <- 0  # Uninfected patches remain uninfected
      }
      
      # Population size doesn't change during fizzle
      N[i, k, 4] <- N[i, k, 3]
      
      #cat("Patch:", i,"  Year:", k, "  Season: 3", "  Population:", N[i, k, 4], "  Status:",I[i,k,4],"  Status change:",I[i, k, 4]-I[i,k,3], "\n")
    }
    cat("Year:", k, "Season: 4", "Infected Patches:", sum(I[, k, 4]), "\n")
    
    # 4. Burnout process
    
    z <- final_size(R0)  # Final size of epidemic
    
    current_year_inf <- 0  # To store total deaths in the current year
    for (i in 1:n_patches) {
      if (I[i, k, 4] == 1) {
        # Reduce population due to disease
        inf <- z * N[i, k, 4]
        N[i, k, 5] <- N[i, k, 4] - inf * D
        current_year_inf <- current_year_inf + inf
        
        # Probability of burnout
        P_b <-burnout_prob(R0,epsilon,N=N[i,k,4])
        #cat("Patch:", i,"  Year:", k, "  Season: 4", "  Population:", N[i, k, 4],"  Pb:", P_b, "\n")
        
        I[i, k, 5] <- ifelse(runif(1) > P_b, 1, 0)
      } else {
        # No change in uninfected patches
        N[i, k, 5] <- N[i, k, 4]
        I[i, k, 5] <- 0
      }
      #cat("Patch:", i,"  Year:", k, "  Season: 4","  Population:", N[i, k, 5], "  Status:",I[i,k,5],"  Status change:",I[i, k, 5]-I[i,k,4], "\n")
    }
    cat("Year:", k, "Season: 5", "Infected Patches:", sum(I[, k, 5]), "\n","\n")
    
    # Record total inf for the current year
    total_inf[k] <- current_year_inf
    
    # Set initial conditions for next year
    if (k < n_years) {
      N[, k+1, 1] <- N[, k, 5]
      I[, k+1, 1] <- I[, k, 5]
    }
  }
  
  return(list(N = N, I = I, total_inf = total_inf))
}

# Function to plot simulation results
plot_simulation_results <- function(result, n_patches = 10) {
  # Extract results
  N <- result$N
  I <- result$I
  total_inf <- result$total_inf
  n_years <- dim(N)[2]
  
  # Calculate summary statistics
  total_population <- apply(N[, , 1], 2, sum)
  infected_patches <- apply(I[, , 1], 2, sum)
  
  # Create plots
  par(mfrow = c(3, 1))
  
  # Plot total population
  plot(1:n_years, total_population, type = "l", col = "blue", 
       main = "Total Population Over Time", xlab = "Year", ylab = "Population Size")
  
  # Plot number of infected patches
  plot(1:n_years, infected_patches, type = "l", col = "red", 
       main = "Number of Infected Patches", xlab = "Year", ylab = "Count")
  
  # Plot total number of infected hosts
  plot(1:n_years, total_inf, type = "l", col = "purple",
       main = "Total inf Over Time", xlab = "Year", ylab = "inf")
}

# Run simulation
set.seed(123) 
result <- simulate_metapopulation(R0=2.5,alpha=5e-5,K=1e6,r=0.5,c=0.2,epsilon=0.02)

# Plot results
plot_simulation_results(result)