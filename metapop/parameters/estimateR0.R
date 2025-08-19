n=1e6
K_f    <- rgamma(n, shape=5, scale=1.4)
betar <- rgamma(n, shape=2, scale=0.25)
life_f <- rgamma(n, shape=4, scale=1)

res <- K_f * betar * life_f
res <- res[res<100]

h<-hist(res,breaks=200,xlab="R0",freq=FALSE)
max_idx <- which.max(h$counts)
max_x <- h$mids[max_idx]
#cat(max_idx,max_x, "\n")
