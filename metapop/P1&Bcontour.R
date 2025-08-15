library(burnout)

getDFE<-function(c0=0.5,
                 S=1e6,
                 R0=3,
                 D=1,
                 r=0.5
) {
  
  z=final_size(R0)
  N_DFE=c0*S/(1-(1+r)*(1-z*D)*(1-c0))
  return(N_DFE)}


# P1_at_DFE<-function(R0=2,epsilon=0.01){
#   N_DFE<-getDFE(R0=R0)
#   P1<-P1_prob(R0,epsilon=epsilon,k=1,N=N_DFE)
# }

thresholdB<-function(R0=2,epsilon=0.01){
  N_DFE<-getDFE(R0=R0)
  P1<-P1_prob(R0,epsilon=epsilon,k=1,N=N_DFE)
  B=1/P1
}



plot_P1 <- function(P_func=P1_at_DFE, x_start = 1, x_end = 5, step = 0.01) {
  x_vals <- seq(x_start, x_end, by = step)
  y_vals <- P_func(x_vals)
  
  plot(x_vals, y_vals, type = "l", col = "blue", lwd = 2,
       xlab = "R0", ylab = "Persistence prob", main = "change of P1 with R0")
}

#plot_P1()

plot_P1_2d <- function(R0_range = c(1, 5), epsilon_range = c(1e-4, 0.02),
                       R0_steps = 100, epsilon_steps = 100) {
  
  R0_vals <- seq(R0_range[1], R0_range[2], length.out = R0_steps)
  eps_vals <- seq(epsilon_range[1], epsilon_range[2], length.out = epsilon_steps)
  
  P1_mat <- outer(R0_vals, eps_vals, Vectorize(function(R0, eps) {
    P1_at_DFE(R0 = R0, epsilon = eps)
  }))
  
  filled.contour(R0_vals, eps_vals, P1_mat,
                 xlab = "R0",
                 ylab = "epsilon",
                 main = "P1 with R0 and epsilon",
                 color.palette = terrain.colors)
}

#plot_P1_2d()

plot_P1_2d <- function(R0_range = c(1, 3), epsilon_range = c(1e-4, 0.03),
                       R0_steps = 100, epsilon_steps = 100) {
  
  R0_vals <- seq(R0_range[1], R0_range[2], length.out = R0_steps)
  eps_vals <- seq(epsilon_range[1], epsilon_range[2], length.out = epsilon_steps)
  
  # compute P1 matrix
  P1_mat <- outer(R0_vals, eps_vals, Vectorize(function(R0, eps) {
    P1_at_DFE(R0 = R0, epsilon = eps)
  }))
  
  # plot with filled.contour
  logP1_mat <- log10(P1_mat)
  levels <- seq(-12, -1, by = 0.5)
  
  filled.contour(eps_vals, R0_vals, t(logP1_mat),
                 levels = levels,
                 xlab = "epsilon",
                 ylab = "R0",
                 main = "P1 with epsilon and R0",
                 color.palette = terrain.colors)
}

#plot_P1_2d()

plot_B <- function(R0_range = c(1, 3), epsilon_range = c(1e-4, 0.03),
                   R0_steps = 100, epsilon_steps = 100) {
  
  R0_vals <- seq(R0_range[1], R0_range[2], length.out = R0_steps)
  eps_vals <- seq(epsilon_range[1], epsilon_range[2], length.out = epsilon_steps)
  
  # compute P1 matrix
  P1_mat <- outer(R0_vals, eps_vals, Vectorize(function(R0, eps) {
    thresholdB(R0 = R0, epsilon = eps)
  }))
  
  # plot with filled.contour
  logP1_mat <- log10(P1_mat)
  levels <- seq(-1, 12, by = 0.5)
  
  filled.contour(eps_vals, R0_vals, t(logP1_mat),
                 levels = levels,
                 xlab = "epsilon",
                 ylab = "R0",
                 main = "threshold B with epsilon and R0",
                 color.palette = terrain.colors)
}

plot_B()