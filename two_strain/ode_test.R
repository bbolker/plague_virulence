library(deSolve)
library(dplyr)
library(purrr)
library(future)
library(progressr)
library(furrr)
plan(multicore, workers = 8)

## FIXME: use odin, integrate on log scale, etc. ...
## this is probably good enough
sir2 <- function(t, y, parms, ...) {
  with(as.list(c(parms, y)),
       {
         incidence <- S*c(R01*I1, R02*I2)
         list(grad = c(S = -sum(incidence),
                       I1tot = incidence[1],
                       I2tot = incidence[2],
                       I1 = incidence[1] - I1,
                       I2 = incidence[2] - I2))
       })
}

intfun <- function(allparms, last_only = TRUE) {
  with(as.list(allparms), {
    res <- ode(y = c(S = 1 - I10 - I20,
                     I1tot = 0,
                     I2tot = 0,
                     I1 = I10,
                     I2 = I20),
               parms = c(R01 = R01, R02 = R02),
               times = seq(0, 50, by = 0.1),
               func = sir2)
    if (last_only) return(res[nrow(res),])
    res
  })
}

dd <- intfun(list(R01 = 2, R02 = 1.5, I10 = 0.01, I20 = 0.02), last_only = FALSE)
par(las = 1)
matplot(dd[,1], dd[,-1], type = "l", log = "y", xlab = "", ylim = c(1e-6, 10))

parvals <- expand.grid(R01 = seq(1.1, 3, by = 0.1),
                       R02 = seq(1.1, 3, by = 0.1),
                       I10 = 10^seq(-3, -1, by = 0.25),
                       I20 = 10^seq(-3, -1, by = 0.25)) |>
  filter(R01 > R02)

## 1.25 seconds for 100 rows
system.time(
  apply(parvals[1:100,], 1, intfun)
)

## 3.2 minutes?
(nrow(parvals)/100*1.25)/60

nn <- nrow(parvals)
res <- with_progress({
  p <- progressor(steps = nn)
  future_map(seq(nn),
             function(i) {
               p()
               intfun(parvals[i,])
             })
}) |>
  do.call(what = "rbind")


res2 <- cbind(parvals, res)
saveRDS(res2, file = "outputs_2strain/finalsize.rds")
## now compare with Yuyang's approximations ((R01+R02)/2 final size, ??)
