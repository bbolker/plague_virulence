#' Function to simulate metapopulation dynamics with infection
#' @param n_patches Number of patches                            
#' @param n_years Number of years to simulate                 
#' @param K Carrying capacity per patch                    
#' @param r Growth rate                                   
#' @param c0 Migration probability                         
#' @param alpha Scaling constant for infection probability   
#' @param D Disease-induced mortality                       
#' @param R0 basic reproduction number                   
#' @param epsilon infectious period (scaled to host generation time)
#' @param initial_infected Number of initially infected patches         
#' @param initial_pop_ratio Initial population as ratio of carrying capacity
#' @return A list containing three matrices, each (n_years x n_patches):
#' total_pops (total host population), infected_patches, total_inf
simulate_metapopulation <- function(
    n_patches = 100,       
    n_years = 1000,        
    K = 1e6,               
    r = 0.5,               
    c0 = 0.2,               
    alpha = 5e-5,          
    D = 1,                 
    R0 = 2.5,              
    epsilon=0.02,          
    initial_infected = 10, 
    initial_pop_ratio = 1  
) {
  if (!require("burnout", quietly = TRUE)) stop("please install burnout package")
  # Initialize arrays to store population size and infection status
  # Dimensions: [patch, year, season]
  # Season: 1=beginning, 2=after growth, 3=after colonization, 4=after fizzle, 5=after burnout
  dn <- list(patch = 1:n_patches, year = 1:n_years,
             season = c("beginning", "after_growth", "after_colonization", "after_fizzle", "after_burnout"))
  N <- array(0, dim = c(n_patches, n_years, 5), dimnames = dn)
  I <- array(0, dim = c(n_patches, n_years, 5), dimnames = dn)
  
  # Initialize an array to store total deaths for each year
  total_inf <- numeric(n_years)
  
  # Set initial conditions
  N[, 1, "beginning"] <- initial_pop_ratio * K  # Initial population in all patches
  
  # Randomly select initial infected patches
  infected_patches <- sample(1:n_patches, initial_infected)
  I[infected_patches, 1, "beginning"] <- 1
  
  # Simulate over years
  for (k in 1:n_years) {
    
    #cat("Year:", k, "Season: 1", "Infected Patches:", sum(I[, k, 1]), "\n")
    
    # Within each year, apply the four processes
    
    # 1. Growth process
    N[, k, "after_growth"] <- N[, k, "beginning"] + r * N[, k, "beginning"] * (1 - N[, k, "beginning"] / K)
    I[, k, "after_growth"] <- I[, k, "beginning"]  # Infection status doesn't change during growth
    #cat("Year:", k, "Season: 2", "Infected Patches:", sum(I[, k, 2]), "\n")
    
    # 2. Colonization process
    
    # 2.1 Population movement
    total_pop <- sum(N[, k, "after_growth"])
    N[, k, "after_colonization"] <- (1 - c0) * N[, k, "after_growth"] + c0 * total_pop / n_patches
    
    # 2.2 Infection spread
    if (k == 1) {
      I[, k, "after_colonization"] <- I[, k, "after_growth"]  # First year: infection status remains unchanged during colonization
    } else {
      P_I <-1-exp(-alpha * c0 * total_inf[k - 1] / n_patches)  # Calculate infection probability based on previous year's deaths
      
      I[, k, "after_colonization"] <- I[, k, "after_growth"]
      susceptible_mask <- I[, k, "after_growth"] == 0
      I[susceptible_mask, k, "after_colonization"] <- rbinom(sum(susceptible_mask), 1, P_I)
      
    }
    #cat("Year:", k, "Season: 3", "Infected Patches:", sum(I[, k, 3]), "\n")
    
    # 3. Fizzle process
    P_f <- 1 / R0  # Fizzle probability
    
    # Population doesn't change during fizzle
    N[, k, "after_fizzle"] <- N[, k, "after_colonization"]
    
    infected_mask <- I[, k, "after_colonization"] == 1
    I[, k, "after_fizzle"] <- 0  # Initialize all as 0
    I[infected_mask, k, "after_fizzle"] <- rbinom(sum(infected_mask), 1, 1 - P_f)
    
    #cat("Year:", k, "Season: 4", "Infected Patches:", sum(I[, k, 4]), "\n")
    
    # 4. Burnout process
    
    z <- final_size(R0)  # Final size of epidemic
    
    N[, k, "after_burnout"] <- N[, k, "after_fizzle"]
    I[, k, "after_burnout"] <- 0
    
    infected_mask <- I[, k, "after_fizzle"] == 1
    
    # Process infected patches only
    if (sum(infected_mask) > 0) {
      # Calculate infections for infected patches only
      inf_infected <- z * N[infected_mask, k, "after_fizzle"]
      
      # Reduce population due to disease
      N[infected_mask, k, "after_burnout"] <- N[infected_mask, k, "after_fizzle"] - inf_infected * D
      
      # Record total infections for the current year
      total_inf[k] <- sum(inf_infected)
      
      # Calculate burnout probabilities 
      P_b_vec <- burnout_prob(R0, epsilon, N = N[infected_mask, k, "after_fizzle"])
      
      # Generate burnout outcomes
      I[infected_mask, k, "after_burnout"] <- rbinom(sum(infected_mask), 1, 1 - P_b_vec)
    } else {
      # No infected patches
      total_inf[k] <- 0
    }
    
    #cat("Year:", k, "Season: 5", "Infected Patches:", sum(I[, k, 5]), "\n","\n")
    
    # Set initial conditions for next year
    if (k < n_years) {
      N[, k+1, 1] <- N[, k, "after_burnout"]
      I[, k+1, 1] <- I[, k, "after_burnout"]
    }
  }
  
  return(list(N = N, I = I, total_inf = total_inf))
}

## default parameters
params0 <- list(n_patches = 100,
                n_years = 1000,
                K = 1e6,
                r = 0.5,
                c0 = 0.2,
                alpha = 5e-5,
                D = 1,
                R0 = 2.5,
                epsilon=0.02)

#' multiple replicate simulations, same parameters
#' @param nsim number of replicate simulations
#' @param params parameter list
#' @param verbose print progress info? (doesn't work for parallel sims)
#' @param ncores number of parallel cores
#' @param seed random-number seed
#' @return 
mult_sim_mp <- function(nsim, params, verbose = FALSE,
                        ncores = getOption("sim.ncores", 4),
                        seed = NULL) {
  require("parallel", quietly = TRUE)
  if (ncores>1) {
    cl <- makeCluster(ncores)
    clusterExport(cl, c("simulate_metapopulation"))
    clusterExport(cl, c("params"), envir = environment())
    clusterEvalQ(cl, "library(burnout)")
    if (!is.null(seed)) clusterSetRNGStream(cl, seed)
    results <- parLapply(cl, 1:nsim,
                         function(i) {
                           if (verbose) cat(".")
                           do.call(simulate_metapopulation, params)
                         })
    stopCluster(cl)
    rm(cl)
  } else {
    results <- replicate(nsim, {
      if (verbose) cat(".")
      do.call(simulate_metapopulation, params)
    }, simplify = FALSE)
  }
  total_pops <- sapply(results, function(res) {
    apply(res$N[, , 1], 2, sum)
  })
  infected_patches <- sapply(results, function(res) {
    apply(res$I[, , 1], 2, sum)
  })
  total_inf <- sapply(results, function(res) {
    res$total_inf
  })
  tibble::lst(total_pops, infected_patches, total_inf)
}

## TO DO:
##  progress bars?
##  different summaries?

plotfun1 <-  function(res, quasi_eq = FALSE) {
  
  main_labs <- c("Total Population over Time",
                 "Number of Infected Patches",
                 "Total infections over Time")
  
  y_labs <- c("Total Population", "Infected Patches", "Total Infections")
  
  col_vec <- c("blue", "red", "purple")
  
  vars <- c("total_pops", "infected_patches", "total_inf")
  
  par(mfrow = c(3, 1), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  
  pfun <- function(var, ...,  meancol = "blue") {
    x <- res[[var]]
    if (quasi_eq) {
      x[res[["infected_patches"]] == 0] <- NA
    }
    matplot(x, type = "l", lty = 1, col = "grey",
            ..., xlab = "Year")
    lines(rowMeans(x, na.rm = TRUE), col = meancol, lwd = 3)
  }
  
  mapply(pfun, var = vars, main = main_labs, ylab = y_labs, meancol = col_vec)
  ## pfun(x$total_pops, main = "Total Population over Time", ylab = "Total Population",
  ##      meancol = "blue")
  
  ## pfun(x$infected_patches, main = "Number of Infected Patches", ylab = "Infected Patches",
  ##      meancol = "red")
  
  ## pfun(x$total_inf, main = "Total Infections over Time", ylab = "Total Infections",
  ##      meancol = "purple")
  
  par(mfrow = c(1, 1))
}


##' average last N steps of state variables
##' @param x 'raw' sim output
sumfun1 <- function(x, nsteps = 100) {
  nyr <- nrow(x[[1]])
  yr_vec <- seq(nyr-nsteps, nyr)
  ## average across pops and years
  mfun <- function(y) mean(y[yr_vec,], na.rm = TRUE)
  r1 <- sapply(x, mfun)
  ## quasi-equilibrium version
  ## replace values from extinct metapopulations (inf_patches == 0) with NA,
  ##  compute mean with na.rm = TRUE (as above)
  x_qe <- lapply(x, function(z) {z[x$infected_patches==0] <- NA; z})
  r2 <- sapply(x_qe, mfun)
  names(r2) <- paste0(names(r2), "_qe")
  c(r1, r2)
}
