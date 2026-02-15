#' Function to simulate metapopulation dynamics with infection(poisson number version)
#' @param n_patches Number of patches
#' @param n_years Number of years to simulate
#' @param K Carrying capacity per patch
#' @param r Growth rate
#' @param c0 Migration probability (host movement), corresponds to c in the .tex
#' @param eta Scaling constant linking host movement to propagule export (tex: \eta)
#' @param kappa Infection-introduction intensity per propagule 
#' @param D Disease-induced mortality (tex: D; often 1)
#' @param R0 basic reproduction number
#' @param initial_infected Number of initially infected patches
#' @param initial_pop_ratio Initial population as ratio of carrying capacity
#' @return A list containing arrays:
#' N (n_patches x n_years x 5), I (n_patches x n_years x 5), total_inf (n_years)
simulate_metapopulation <- function(
    n_patches = 100,
    n_years = 20,
    K = 1e6,
    r = 0.5,
    c0 = 0.2,
    eta = 0.5,
    kappa = 1e-5,
    D = 1,
    R0 = 2.5,
    initial_inf_ratio = 0.1,
    initial_pop_ratio = 1,
    early_stop = TRUE               
) {
  if (!require("burnout", quietly = TRUE)) stop("please install burnout package")
  
  dn <- list(patch = 1:n_patches, year = 1:n_years,
             season = c("beginning", "after_growth", "after_colonization", "after_fizzle", "after_burnout"))
  N <- array(0, dim = c(n_patches, n_years, 5), dimnames = dn)
  I <- array(0, dim = c(n_patches, n_years, 5), dimnames = dn)
  
  total_inf <- rep(NA_real_, n_years)   
  S <- array(NA_real_, dim = c(n_patches, n_years))  
  extinct_year <- NA_integer_        
  
  N[, 1, "beginning"] <- initial_pop_ratio * K
  
  infected_patches <- sample(1:n_patches, initial_inf_ratio * n_patches)
  I[infected_patches, 1, "beginning"] <- 1
  
  P_f <- 1 / R0
  z <- final_size(R0)
  g <- pmin(1, eta * c0)
  
  for (k in 1:n_years) {
    
    # 1. Growth
    N[, k, "after_growth"] <- N[, k, "beginning"] + r * N[, k, "beginning"] * (1 - N[, k, "beginning"] / K)
    I[, k, "after_growth"] <- I[, k, "beginning"]
    
    # 2. Colonization
    total_pop <- sum(N[, k, "after_growth"])
    N[, k, "after_colonization"] <- (1 - c0) * N[, k, "after_growth"] + c0 * total_pop / n_patches
    
    I[, k, "after_colonization"] <- I[, k, "after_growth"]
    
    if (k > 1) {
      lambda_G <- ((1 - g) * S[, k-1] + (g / n_patches) * total_inf[k-1])
      
      lambda_X <- kappa * lambda_G
      X_new <- rpois(n_patches, lambda = pmax(lambda_X, 0))
      I[, k, "after_colonization"] <- I[, k, "after_growth"] + X_new
    }
    
    # 3. Fizzle
    N[, k, "after_fizzle"] <- N[, k, "after_colonization"]
    
    m <- I[, k, "after_colonization"]
    P_est <- 1 - (P_f ^ m)
    I[, k, "after_fizzle"] <- rbinom(n_patches, 1, pmax(pmin(P_est, 1), 0))
    
    # 4. Major epidemic 
    N[, k, "after_burnout"] <- N[, k, "after_fizzle"]
    I[, k, "after_burnout"] <- 0
    
    infected_patches <- I[, k, "after_fizzle"]
    S[, k] <- z * D * N[, k, "after_fizzle"] * infected_patches
    N[, k, "after_burnout"] <- N[, k, "after_fizzle"] - S[, k]
    total_inf[k] <- sum(S[, k])
    
    # early stop if extinct
    if (early_stop && sum(infected_patches) == 0) {
      extinct_year <- k
      
      if (k < n_years) {
        N[, (k+1):n_years, ] <- NA_real_
        I[, (k+1):n_years, ] <- NA_real_
        S[, (k+1):n_years] <- NA_real_
        total_inf[(k+1):n_years] <- NA_real_
      }
      break
    }
    
    # Next year init
    if (k < n_years) {
      N[, k+1, "beginning"] <- N[, k, "after_burnout"]
      I[, k+1, "beginning"] <- I[, k, "after_burnout"]
    }
  }
  
  return(list(N = N, I = I, S = S, total_inf = total_inf, extinct_year = extinct_year))
}

## default parameters 
params0 <- list(
  n_patches = 100,
  n_years = 1000,
  K = 1e6,
  r = 0.5,
  c0 = 0.2,
  eta = 1,
  kappa = 1e-5,
  D = 1,
  R0 = 2.5,
  initial_inf_ratio = 0.1,
  initial_pop_ratio = 1
)

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
    lp <- .libPaths()
    clusterExport(cl, "lp")
    clusterEvalQ(cl, .libPaths(lp))
    
    clusterExport(cl, c("simulate_metapopulation"))
    clusterExport(cl, c("params"), envir = environment())
    
    clusterEvalQ(cl, quote(library(burnout)))
    
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
    apply(res$I[, , 3], 2, sum)   ##infection status after colonization
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