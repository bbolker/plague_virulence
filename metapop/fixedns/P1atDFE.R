library(burnout)

getDFE<-function(c0=0.01, ## was: 0.5
                 S=1e6, ## smaller population makes the decreasing range larger
                 R0=3,
                 D=1,
                 r=0.5
) {
  z=final_size(R0)
  N_DFE=c0*S/(1-(1+r)*(1-z*D)*(1-c0))
  return(N_DFE)}

P1_at_DFE<-function(R0=2,epsilon=0.02){
  N_DFE <- getDFE(R0=R0)
  P1 <- P1_prob(R0,epsilon=epsilon,k=1,N=N_DFE)
}

n <- 101
epsvec <- seq(1e-4, 0.03, length.out = n) ## warnings if eps=0
R0vec <- seq(1, 2.5, length.out = n)
grad_R0 <- persist <- matrix(NA, length(epsvec), length(R0vec))
for (i in seq_along(epsvec)) {
  for (j in seq_along(R0vec)) {
    persist[i,j] <- P1_at_DFE(R0 = R0vec[j], eps = epsvec[i]) 
    ## cat(" eps=", epsvec[i], " R0=", R0vec[j], " P1=", persist[i,j], "\n")
  }
}
dR0 <- diff(R0vec)[1]
for (j in 2:length(R0vec)) {
  grad_R0[,j] <- (persist[,j]-persist[,j-1])/dR0
}



par(las=1, bty = "l")
  image(epsvec, R0vec, grad_R0>0, col = c(adjustcolor("blue", alpha = 0.3), "#FFFFFF"),
      main = "P1 at DFE")
contour(epsvec, R0vec, log10(persist),
        levels = c(-12, -8, -4, -2, -1, log10(0.2), log10(0.3)),
        add = TRUE)
  contour(epsvec, R0vec, grad_R0, level = 0, add = TRUE, col = 2, lwd = 2)

dd <- data.frame(eps = rep(epsvec, length(R0vec)),
                 R0 = rep(R0vec, each = length(epsvec)), grad_R0 = c(grad_R0),
                 persist = c(persist))

## hack only the part we want of the contour
cc <- contourLines(epsvec, R0vec, grad_R0, level = 0)
cc <- do.call(data.frame, cc[[1]][c("x", "y")])[-(1:138),]


library(ggplot2); theme_set(theme_bw(base_size=16))
ggplot(dd, aes(eps, R0)) +
  geom_raster(aes(fill = persist)) +
  scale_fill_viridis_c(name = "equilibrium\nmetapop\noccupancy") +
  annotate(x = 0.008, y = 1.8,
           label = expression(Delta*P/Delta*R[0] < 0),
           geom = "label", size = 10,
           fill = "gray") +
  geom_contour(aes(z = persist), colour = "gray") +
    geom_path(data = cc, aes(x, y), colour = "red", linewidth = 2) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = expression("disease generation time "*(epsilon)),
       y = expression(R[0]))

ggsave("analytic_MP_occ.png")
