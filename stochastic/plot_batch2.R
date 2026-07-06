library(ggplot2); theme_set(theme_bw())
library(tidyverse)

dd <- readRDS("outputs/sim_batch3.rds")
print(nrow(dd))

pivfun <- function(x) {
  res <- tidyr::pivot_longer(x, contains("_"))
  if ("alphavec" %in% names(res)) {
    res <- res |> mutate(alphavec_f =
                           ordered(alphavec, levels = unique(alphavec), labels = signif(unique(alphavec), 2)))
  }
  res
}

dd_long <- pivfun(dd) |> rename(R0 = "R0vec", K = "Kvec")

out_dir <- "sim_batch3_plot"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# produce plot for each distinct K
Ks <- sort(unique(dd_long$K))
for (k in Ks) {
  ddk <- dd_long %>% filter(K == k)

  gg1 <- ggplot(ddk, aes(R0, value, colour = alphavec_f)) +
    geom_point() +
    ## GAM gives nicer plots than default loess smoothing, but fails in a few cases
    ## (because quasi-eq doesn't leave enough points etc.)
    geom_smooth(method = "gam", formula = y ~ s(x)) +
    facet_wrap(~ name, scale = "free", ncol = 2)
  
  p <- gg1 + ggtitle(paste("K =", format(k, scientific = FALSE)))
  
  fname <- file.path(out_dir, paste0("plot_batch3_K_", format(k, scientific = FALSE), ".pdf"))
  ggsave(fname, p, width = 10, height = 5)
  message("Saved: ", fname)
  
}