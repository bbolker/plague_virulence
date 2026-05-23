library(odin)
library(dde) ## odin insists on this
library(tidyr)
library(dplyr)
library(ggplot2)

K <- 1e4
r <- 0.125
n_patch  <- 100
n_strains <- 2
nt       <- 1000
alpha <- 1e-3


## patch-level parameters (vectors, length n_patch)
K_vec <- rep(K, n_patch)
r_vec <- rep(r, n_patch)

set.seed(101)
## strain x patch parameters: matrices [n_patch, n_strains], rows = patches, cols = strains
beta_vec  <- c(1.5, 2.5)
I_ini_mat <- matrix(rpois(n_strains*n_patch, lambda = 10),
                    ncol = n_strains)

I_ini_mat[seq(n_patch/2),] <- 0  ## test of patch coupling

S_ini_vec <- K_vec - rowSums(I_ini_mat)

odin_file <- here::here("metapop/odin", "odin_twostrain_npatch.R")
gen <- odin::odin(odin_file)

args <- tibble::lst(beta    = beta_vec,
                    r       = r_vec,
                    K       = K_vec,
                    I_ini   = I_ini_mat,
                    S_ini   = S_ini_vec,
                    alpha,
                    n_patch)

mod <- do.call(gen$new, args)

system.time(res <- mod$run(0:nt))

conv_fun <- function(x) {
  ret <- as.data.frame.table(x[,!colnames(res) %in% "step"]) |>
    as_tibble() |>
    transmute(step = as.numeric(Var1), ## grab factor levels
              state = case_when(
                grepl("^S", Var2) ~ gsub("\\[", ".", gsub("\\]", "", Var2)),
                grepl("^I", Var2) ~ gsub("\\[([0-9]+),([0-9]+)\\]", "\\2.\\1", Var2)),
              value = Freq) |>
    mutate(patch = as.numeric(sub(".*\\.", "", state)),
           state = sub("\\..*", "", state)) |>
    arrange(step, patch)
  return(ret)
}

cc <- conv_fun(res)


ggplot(cc, aes(step, value, colour = state)) +
  geom_line(aes(group = interaction(state, patch)))  +
  scale_y_log10()

ggsave("odin_twostrain_run_patch.png")
