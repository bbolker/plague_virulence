library(burnout)

logN <- seq(1, 7, length.out = 100)

Pb <- burnout_prob(R0=2.5, epsilon=0.02, N=10^logN)

plot(logN, Pb, type = "l", col = "blue", lwd = 2)