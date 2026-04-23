pkgs <- c("gsl", "tidyverse", "ggrastr", "biscale", "cowplot",
          "deSolve", "future", "progressr", "furrr")
i1 <- installed.packages()
pkgs <- setdiff(pkgs, rownames(i1))
if (length(pkgs)>0) {
  install.packages(pkgs)
}
