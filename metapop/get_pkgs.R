options(repos = (CRAN = "https://cloud.r-project.org"))
pkgs <- c(
  "biscale"
 , "colorspace"
 , "cowplot"
 , "deSolve"
 , "furrr"
 , "future"
 , "ggrastr"
 , "gsl"
 , "optparse"
 , "patchwork"
 , "progressr"
 , "tidyverse"
 , "tinysnapshot"
 , "tinytest"
)

i1 <- installed.packages()
pkgs <- setdiff(pkgs, rownames(i1))
if (length(pkgs)>0) {
  install.packages(pkgs)
}
remotes::install_github("davidearn/burnout")
remotes::install_github("canmod/macpan2")
