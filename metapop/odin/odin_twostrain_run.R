library(odin)
library(dde) ## not sure why we need this?

## test expm1
## -expm1(-x) = 1-exp(-x)
## stopifnot(all.equal(-expm1(-0.2), 1-exp(-0.2)))

twostrain_generator <-odin::odin("./odin_twostrain0.R")
x <- twostrain_generator$new(beta = c(0.2, 0.19),
                        gamma = c(0.1, 0.1),
                        I_ini = c(5, 2))

set.seed(101)
res <- x$run(0:1000)
ext <- which(rowSums(is.na(res))==0 & res[,"I[1]"]+res[,"I[2]"] == 0)[1]
res <- res[1:(ext-1),]

pdf("odin_twostrain_run.Rout.pdf")
matplot(res[,1], res[,-1], type = "l", log = "y")
legend("right", legend = colnames(res)[-1], lty = 1, col = 1:4)
dev.off()
