
library(burnout)
## Figure 5 computations

n <- 101
epsvec <- seq(1e-4, 0.05, length.out = n) ## warnings if eps=0
R0vec <- seq(1, 2.5, length.out = n)
B<-grad_R0 <- persist <- matrix(NA, length(epsvec), length(R0vec))
for (i in seq_along(epsvec)) {
   for (j in seq_along(R0vec)) {
     persist[i,j] <- P1_prob(R0 = R0vec[j], eps = epsvec[i],N=1e5)/not_fizzle_prob(R0=R0vec[j]) ## subject to N
   }
}
dR0 <- diff(R0vec)[1]
for (j in 2:length(R0vec)) {
   grad_R0[,j] <- (persist[,j]-persist[,j-1])/dR0
}

par(las=1, bty = "l")
image(epsvec, R0vec, grad_R0>0, col = c(adjustcolor("blue", alpha =
0.3), "#FFFFFF"))
contour(epsvec, R0vec, persist,
         levels = seq(0, 2, by = 0.1),
         add = TRUE)
contour(epsvec, R0vec, grad_R0, level = 0, add = TRUE, col = 2, lwd = 2)

