if (!file.exists("polyfit.rda")) {
  stop("please `make polyfit.rda`")
}
load("polyfit.rda")

#' Two-Strain Metapopulation Simulation with Vectorized Logic
#' @param n_patches Number of patches in the metapopulation
#' @param n_years Number of years to simulate
#' @param K Carrying capacity per patch
#' @param r Growth rate of the host population
#' @param c0 Migration probability (proportion of hosts moving)
#' @param nu Scaling constant (dead hosts -> total environmental fleas)
#' @param rho Carrying capacity of fleas per migrating host (The Bottleneck)
#' @param alpha Transmission efficiency per flea (infectivity)
#' @param D Disease-induced mortality (usually 1 for plague)
#' @param R01 Basic reproduction number for Strain 1 (Resident)
#' @param R02 Basic reproduction number for Strain 2 (Invader)
#' @param invade_year Year when Strain 2 is introduced
#' @param initial_inf_ratio_1 Proportion of patches initially infected with Strain 1
#' @param initial_inf_ratio_2 Proportion of patches to seed Strain 2 during invasion
#' @param initial_pop_ratio Initial population as ratio of K
#' @param early_stop If TRUE, stops simulation when both strains are extinct
#' @param seed random-number seed
#' @param coinf_approx approximation to use for coinfection dynamics
#' @return A list containing N (Host dynamics), I (Indicators), S (Mortality), and Total Infection counts
simulate_metapopulation_2strain <- function(
    n_patches = 100,
    n_years = 500,
    K = 1e6,
    r = 0.5,
    c0 = 0.2,
    nu = 5,        
    rho = 3,       
    alpha = 5e-6,   
    D = 1,
    R01 = 2.0,     
    R02 = 1.8,    
    invade_year = 100,      
    initial_inf_ratio_1 = 0.1, 
    initial_inf_ratio_2 = 0.05, 
    initial_pop_ratio = 1,
    early_stop = FALSE,
    seed = NULL,
    coinf_approx = c("yy", "polyfit")
) {
  if (!require("burnout", quietly = TRUE)) stop("Please install the 'burnout' package")
  if (!is.null(seed)) set.seed(seed)
  coinf_approx <- match.arg(coinf_approx)

  z_vec <- c(burnout::final_size(R01),
             burnout::final_size(R02),
             ## only needed for "YY" approximation
             burnout::final_size((R01+R02)/2))


  # Vectorized helper for the Final Size Equation (1 - z = exp(-R0 * z))
  z_func_vec <- function(R0_vec) {
    res <- numeric(length(R0_vec))
    idx <- R0_vec > 1
    if(any(idx)) {
      # Apply final_size to each valid R0 entry
      res[idx] <- sapply(R0_vec[idx], burnout::final_size)
    }
    return(res)
  }

  ## FIXME: compute final sizes *once* for strain-1-only (R01),
  ## strain-2-only (R02), coinfected patches (either YY approximation
  ## or approximation derived from polynomial regression)

  ## FIXME: named dimnames for N, I, S, etc.
  
  ## Initialize recording structures

  ## Total host population across 4 stages
  dn0 <- list(patch = seq.int(n_patches), year = seq.int(n_years))
  N <- array(0, dim = c(n_patches, n_years, 4),
             dimnames = c(dn0, list(stage = c("begin", "after_growth", "after_colonization", "end"))))

  smat <- function() matrix(0, nrow = n_patches, ncol = n_years, dimnames = dn0)
  I1 <- smat()  # Indicator for Strain 1
  I2 <- smat()  # Indicator for Strain 2
  S1 <- smat()  # Deaths from Strain 1
  S2 <- smat() # Deaths from Strain 2
  total_inf1 <- rep(0, n_years) # Global infection count Strain 1
  total_inf2 <- rep(0, n_years) # Global infection count Strain 2

  Z_total <- rep(NA_real_, n_patches)

  # Initial Setup: Seed Resident Strain (Strain 1)
  N[, 1, "begin"] <- initial_pop_ratio * K
  infected_init_1 <- sample(1:n_patches,
                            size = max(1, initial_inf_ratio_1 * n_patches))
  z_init <- z_vec[1]
  S1[infected_init_1, 1] <- z_init * D * N[infected_init_1, 1, "begin"]
  I1[infected_init_1, 1] <- 1 
  total_inf1[1] <- sum(S1[, 1])
  
  ## pre-invasion colonization of strain 2
  b <- rep(0, n_patches)
  
  for (k in 2:n_years) {
    
    # --- Stage 1 -> 2: Host Growth ---
    # Retrieve previous year-end population (Stage 4)
    N_prev_end <- if(k == 2) N[, 1, 1] - S1[, 1] else N[, k-1, 4]
    
    # Logistic growth formula
    N_after_growth <- N_prev_end + r * N_prev_end * (1 - N_prev_end / K)
    N[, k, 2] <- N_after_growth
    
    # --- Stage 2 -> 3: Host Migration (Colonization) ---
    pop_total_network <- sum(N_after_growth)
    N_after_col <- (1 - c0) * N_after_growth + c0 * (pop_total_network / n_patches)
    N[, k, 3] <- N_after_col
    
    # --- Transmission Dynamics: Force of Infection ---
    S_last_total <- S1[, k-1] + S2[, k-1]
    
    # Proportions of strains in the flea pool (q-vectors)
    q1 <- ifelse(S_last_total > 0, S1[, k-1] / S_last_total, 0)
    q2 <- ifelse(S_last_total > 0, S2[, k-1] / S_last_total, 0)
    
    # Calculate Transport Bottleneck (Min function for fleas vs migration capacity)
    E_total <- pmin(nu * S_last_total, rho * c0 * N_after_growth)
    R_total <- nu * S_last_total - E_total # Retained fleas
    
    # Force of Infection 
    lambda1 <- alpha * (q1 * R_total + sum(q1 * E_total) / n_patches) * max(0, 1 - 1/R01)
    lambda2 <- alpha * (q2 * R_total + sum(q2 * E_total) / n_patches) * max(0, 1 - 1/R02)
    
    # Poisson sampling for initial infection seeds (a and b)
    a <- rpois(n_patches, lambda1)
    if (k >= invade_year) b <- rpois(n_patches, lambda2)
    
    # --- Invasion Event: Introduction of Strain 2 ---
    if (k == invade_year) {
      invade_idx <- sample(1:n_patches, max(1, initial_inf_ratio_2 * n_patches))
      b[invade_idx] <- b[invade_idx] + 1 # Force establishment in selected patches
    }
    
    # --- Stage 3 -> 4: Epidemic Outbreak & Resource Partitioning ---
    ind1 <- a > 0 # Indicator for successful establishment
    ind2 <- b > 0
    ind_any <- (a + b) > 0
    denom <- a + b
    
    I1[, k] <- as.numeric(ind1)
    I2[, k] <- as.numeric(ind2)
    
    if (coinf_approx == "polyfit") {
      coinf <- which(ind1 & ind2)
      if (length(coinf) > 0) {
        coinf_res <- sapply(coinf,
                            \(i) pred_outcomes_poly(R01, R02, a[i]/N_after_col[i], b[i]/N_after_col[i]))
      }
    }
    browser()

    # Calculate total death toll Z using vectorized final size
    Z_total[ind1 & !ind2] <- z_vec[1] * D * N_after_col
    Z_total[ind2 & !ind1] <- z_vec[2] * D * N_after_col
    if (coinf_approx == "YY") {
      Z_total[ind1 &  ind2] <- z_vec[3] * D * N_after_col
    } else {
      Z_total[ind1 &  ind2] <- coinf_res[,1] * D * N_after_col
    }
      
    # Partition Z based on initial seed ratios (Galton-Watson allocation)
    # denom is (a + b); we only divide when ind_any is TRUE to avoid 0/0

    if (coinf_approx == "YY") {
      S1[, k] <- ifelse(ind_any, Z_total * a / denom, 0)
      S2[, k] <- ifelse(ind_any, Z_total * b / denom, 0)
    } else {
      stop("not implemented")
    }
    
    # Update end-of-year surviving population
    N[, k, 4] <- N_after_col - (S1[, k] + S2[, k])
    
    # Record global time-series statistics
    total_inf1[k] <- sum(S1[, k])
    total_inf2[k] <- sum(S2[, k])
    
    if (k < n_years) N[, k+1, 1] <- N[, k, 4]
    
    # Survival check: exit if both strains disappear after invasion starts
    if (early_stop && sum(ind_any) == 0 && k > invade_year) break
  }
  
  return(list(
    N = N, 
    I1 = I1, I2 = I2, 
    S1 = S1, S2 = S2, 
    total_inf1 = total_inf1, 
    total_inf2 = total_inf2
  ))
}

#' Wrapper for parallel execution of two-strain model
#' @param nsim Number of simulations
#' @param params List of parameters
#' @param verbose Print progress
#' @param ncores Number of cores
#' @param seed Master random seed
mult_sim_2strain <- function(nsim, params, verbose = FALSE,
                             ncores = getOption("sim.ncores", 4),
                             seed = NULL) {
  require("parallel", quietly = TRUE)
  
  wrapper_fun <- function(i) {
    if (verbose) cat(".")
    do.call(simulate_metapopulation_2strain, params)
  }
  
  if (ncores > 1) {
    cl <- makeCluster(ncores)
    lp <- .libPaths()
    clusterExport(cl, "lp", envir = environment())
    clusterEvalQ(cl, .libPaths(lp))
    clusterExport(cl, c("simulate_metapopulation_2strain", "params"), envir = environment())
    clusterEvalQ(cl, library(burnout))
    if (!is.null(seed)) clusterSetRNGStream(cl, seed)
    
    results <- parLapply(cl, 1:nsim, wrapper_fun)
    stopCluster(cl)
  } else {
    if (!is.null(seed)) set.seed(seed)
    results <- replicate(nsim, wrapper_fun(1), simplify = FALSE)
  }
  
  # Extract components into lists/matrices
  list(
    patches1   = sapply(results, function(res) colSums(res$I1)),
    patches2   = sapply(results, function(res) colSums(res$I2)),
    total_inf1 = sapply(results, function(res) res$total_inf1),
    total_inf2 = sapply(results, function(res) res$total_inf2),
    total_pops = sapply(results, function(res) colSums(res$N[, , 1])) # Get total host pop
  )
}

#' Summary function for two-strain outputs (Standardized Names)
sumfun_2strain <- function(x, nsteps = 100, invade_year = 100) {
  
  nyr <- nrow(x$patches1)  
  nsim <- ncol(x$patches1)
  
  # --- Strain 1 (Resident) Extinction Stats ---
  extinct_idx1 <- apply(x$patches1, 2, function(col) {
    idx <- which(col == 0)
    if (length(idx) > 0) return(idx[1]) else return(NA_integer_)
  })
  extinction_rate_1 <- sum(!is.na(extinct_idx1)) / nsim
  mean_extinct_time_1 <- mean(extinct_idx1, na.rm = TRUE)
  
  # --- Strain 2 (Invader) Extinction Stats ---
  extinct_idx2 <- apply(x$patches2, 2, function(col) {
    idx <- which(col == 0 & seq_along(col) > invade_year)
    if (length(idx) > 0) return(idx[1]) else return(NA_integer_)
  })
  extinction_rate_2 <- sum(!is.na(extinct_idx2)) / nsim
  mean_extinct_time_2 <- mean(extinct_idx2, na.rm = TRUE)
  
  # --- Quasi-Equilibrium Means ---
  yr_vec <- seq(max(1, nyr - nsteps + 1), nyr) 
  
  calc_qe <- function(value_mat, patch_mat) {
    value_mat[patch_mat == 0] <- NA
    mean(value_mat[yr_vec, ], na.rm = TRUE)
  }
  
  infected_patches_1 <- calc_qe(x$patches1, x$patches1)
  infected_patches_2 <- calc_qe(x$patches2, x$patches2)
  total_inf_1 <- calc_qe(x$total_inf1, x$patches1)
  total_inf_2 <- calc_qe(x$total_inf2, x$patches2)
  
  # Duplicate host population for both strains so ggplot can draw them together
  total_pops_1 <- mean(x$total_pops[yr_vec, ], na.rm = TRUE)
  total_pops_2 <- total_pops_1
  
  # Return fully standardized names
  c(
    extinction_rate_1 = extinction_rate_1,
    extinction_rate_2 = extinction_rate_2,
    mean_extinct_time_1 = mean_extinct_time_1,
    mean_extinct_time_2 = mean_extinct_time_2,
    total_pops_1 = total_pops_1,
    total_pops_2 = total_pops_2,
    infected_patches_1 = infected_patches_1,
    infected_patches_2 = infected_patches_2,
    total_inf_1 = total_inf_1,
    total_inf_2 = total_inf_2
  )
}

## Default parameters for the two-strain model
params0_2strain <- list(
  n_patches = 100,
  n_years = 500,
  K = 1e6,
  r = 0.5,
  c0 = 0.2,
  nu = 5,
  rho = 3,
  alpha = 5e-6,
  D = 1,
  R01 = 2.0,
  R02 = 1.9,
  invade_year = 100,
  initial_inf_ratio_1 = 0.1,
  initial_inf_ratio_2 = 0.02,
  initial_pop_ratio = 1,
  early_stop = FALSE
)
