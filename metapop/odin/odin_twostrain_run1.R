library(odin)
library(dde) ## odin insists on this
odin_fn <- here::here("metapop/odin", "odin_twostrain1.R")
twostrain_generator <- odin::odin(odin_fn)

K <- 1e4
I_ini <- c(5, 5)
x <- twostrain_generator$new(beta = c(2.1, 1.1),
                             r = 0.1,
                             K = K,
                             I_ini = I_ini,
                             S_ini = K - sum(I_ini),
                             stochastic = 1)

pfun  <- function(res, ymin  = NULL, ymax = NULL) {
    ymin <- ymin %||% min(res[,-1])
    ymax <- ymax %||% max(res[,-1])
    matplot(res[,1], res[,-1],
            type = "l", log = "y", lwd = 2.5, col = c(1,2,3,4), lty = 1,
            ylim = c(ymin, ymax))
}

set.seed(101)
res <- x$run(0:500)
pfun(res, ymin = 1e-3)

ext <- which(rowSums(is.na(res))==0 & res[,"I[1]"]+res[,"I[2]"] == 0)[1]
res <- res[1:(ext-1),]
print(nrow(res))


pdf("odin_twostrain_run.Rout.pdf")
legend("right", legend = colnames(res)[-1], lty = 1, col = 1:4)
dev.off()
