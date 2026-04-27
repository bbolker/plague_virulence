options(repos = (CRAN = "https://cloud.r-project.org"))
pkgs <- c("gsl", "tidyverse", "ggrastr", "biscale", "cowplot",
          "deSolve", "future", "progressr", "furrr")
i1 <- installed.packages()
pkgs <- setdiff(pkgs, rownames(i1))
if (length(pkgs)>0) {
  install.packages(pkgs)
}
remotes::install_github("davidearn/burnout")
