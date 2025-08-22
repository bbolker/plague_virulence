library(burnout)

logN <- seq(1, 7, length.out = 100)

P1 <- P1_prob(R0=2, epsilon=0.02, k=1, N=10^logN)

plot(logN, P1, type = "l", col = "blue", lwd = 2)
