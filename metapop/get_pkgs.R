options(repos = (CRAN = "https://cloud.r-project.org"))
pkgs <- c(
  "biscale"
 , "colorspace"
 , "cowplot"
 , "dde" ## cryptic odin requirement
 , "deSolve"
 , "furrr"
 , "future"
 , "ggrastr"
 , "gsl"
 , "odin"
 , "optparse"
, "patchwork"
, "pkgbuild" ## cryptic odin requirement
   , "pkgload" ## cryptic odin requirement
 , "progressr"
 , "remotes"
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
