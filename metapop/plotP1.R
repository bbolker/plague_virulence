library(burnout)
## Figure 5 computations

n <- 101
epsvec <- seq(1e-4, 0.02, length.out = n) ## warnings if eps=0
R0vec <- seq(1, 2.5, length.out = n)
grad_R0 <- persist <- matrix(NA, length(epsvec), length(R0vec))
for (i in seq_along(epsvec)) {
  for (j in seq_along(R0vec)) {
    persist[i,j] <- P1_prob(R0 = R0vec[j], eps = epsvec[i]) ## uses N=10^6, I(0)=k=1
  }
}
dR0 <- diff(R0vec)[1]
for (j in 2:length(R0vec)) {
  grad_R0[,j] <- (persist[,j]-persist[,j-1])/dR0
}


par(las=1, bty = "l")
image(epsvec, R0vec, grad_R0>0, col = c(adjustcolor("blue", alpha = 0.3), "#FFFFFF"))
contour(epsvec, R0vec, log10(persist),
        levels = c(-12, -8, -4, -2, -1, log10(0.2), log10(0.3)),
        add = TRUE)
contour(epsvec, R0vec, grad_R0, level = 0, add = TRUE, col = 2, lwd = 2)