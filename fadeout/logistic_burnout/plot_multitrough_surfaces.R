source(file.path("fadeout","logistic_burnout","logistic_burnout_functions.R"))
if(!requireNamespace("ggplot2",quietly=TRUE))stop("ggplot2 required")
od <- file.path("fadeout","logistic_burnout","outputs")
fd <- file.path("fadeout","logistic_burnout","figures")
dir.create(od,recursive=TRUE,showWarnings=FALSE)
dir.create(fd,recursive=TRUE,showWarnings=FALSE)
K <- 10000; I0 <- 1; max_troughs <- 5
grid <- expand.grid(R0=seq(1.05,5,length.out=41),
                    r=seq(.01,.5,length.out=41))
cache_file <- file.path(od,"logistic_multitrough_R0_r_grid.csv")
if(file.exists(cache_file) && !"--recompute" %in% commandArgs(TRUE) &&
   nrow(read.csv(cache_file))==nrow(grid)*max_troughs) {
  d <- read.csv(cache_file)
} else {
res <- vector("list",nrow(grid))
for(i in seq_len(nrow(grid))) {
  if(i%%100==0||i==1)message(i,"/",nrow(grid))
  p <- grid[i,]
  z <- tryCatch(logistic_multitrough_probabilities(
    p$R0,p$r,K,I0,max_troughs,dt=.05,rtol=1e-8,atol=1e-10,
    rel.tol=1e-7),error=function(e)e)
  blank <- data.frame(R0=p$R0,r=p$r,K=K,I0=I0,
    trough_index=seq_len(max_troughs),x_in=NA_real_,q_lineage=NA_real_,
    Q_burnout=NA_real_,P_persistence=NA_real_,
    cumulative_persistence=NA_real_,burnout_at_this_trough=NA_real_,
    m_raw=NA_real_,m_used=NA_real_,entry_status="not_available",
    integration_converged=FALSE,troughs_found=0L,
    termination_status=if(inherits(z,"error"))paste0("error: ",conditionMessage(z))
      else z$termination_status)
  if(!inherits(z,"error")&&nrow(z$trough_table)) {
    tt <- z$trough_table; jj <- tt$trough_index
    blank[jj,c("x_in","q_lineage","Q_burnout","P_persistence",
      "cumulative_persistence","burnout_at_this_trough","m_raw","m_used",
      "entry_status","integration_converged")] <-
      tt[,c("x_in","q_lineage","Q_burnout","P_persistence",
        "cumulative_persistence","burnout_at_this_trough","m_raw","m_used",
        "entry_status","integration_converged")]
    blank$troughs_found <- z$troughs_found
  }
  res[[i]] <- blank
}
d <- do.call(rbind,res)
write.csv(d,cache_file,row.names=FALSE)
}
status <- aggregate(R0~troughs_found+termination_status,d[d$trough_index==1,],
                    length); names(status)[3]<-"n_grid_cells"
write.csv(status,file.path(od,"logistic_multitrough_status_summary.csv"),
          row.names=FALSE)
d$trough <- factor(paste0("P",d$trough_index),levels=paste0("P",1:5))
cap <- paste0("K = 10,000; I(0) = 1; y_BL = r(R0-1)/R0^2. ",
 "P_j is conditional on reaching trough j; early stochastic fizzle is excluded.\n",
 "Grey cells mean that the j-th deterministic boundary entry was not found or defined.")
base_theme <- ggplot2::theme_bw()+ggplot2::theme(
  panel.grid=ggplot2::element_blank(),plot.caption=ggplot2::element_text(hjust=0,size=8))
pj <- ggplot2::ggplot(d,ggplot2::aes(R0,r,fill=P_persistence))+
  ggplot2::geom_tile()+ggplot2::facet_wrap(~trough,ncol=3)+
  ggplot2::scale_fill_viridis_c(limits=c(0,1),na.value="grey70",name="P_j")+
  ggplot2::labs(x="R0",y="r",title="Conditional persistence at successive troughs",
                caption=cap)+base_theme
ggplot2::ggsave(file.path(fd,"logistic_Pj_R0_r_facets.png"),pj,width=11,height=8,dpi=180)
ggplot2::ggsave(file.path(fd,"logistic_Pj_R0_r_facets.pdf"),pj,width=11,height=8)
for(j in 1:3) {
  x <- d[d$trough_index==j,]
  p <- ggplot2::ggplot(x,ggplot2::aes(R0,r,fill=P_persistence))+
    ggplot2::geom_tile()+ggplot2::scale_fill_viridis_c(
      limits=c(0,1),na.value="grey70",name=paste0("P",j))+
    ggplot2::labs(x="R0",y="r",title=paste("Conditional persistence at trough",j),
                  caption=cap)+base_theme
  ggplot2::ggsave(file.path(fd,paste0("logistic_P",j,"_R0_r_heatmap.png")),p,width=9,height=7,dpi=180)
  ggplot2::ggsave(file.path(fd,paste0("logistic_P",j,"_R0_r_heatmap.pdf")),p,width=9,height=7)
}
d$log10Q <- ifelse(is.na(d$Q_burnout),NA,pmax(-12,log10(d$Q_burnout)))
qplot <- ggplot2::ggplot(d,ggplot2::aes(R0,r,fill=log10Q))+
  ggplot2::geom_tile()+ggplot2::facet_wrap(~trough,ncol=3)+
  ggplot2::scale_fill_viridis_c(limits=c(-12,0),na.value="grey70",
                               name="log10(Q_j), clipped")+
  ggplot2::labs(x="R0",y="r",title="Conditional burnout at successive troughs",
                caption=paste0(cap,"\nValues below 1e-12 are clipped only for display."))+
  base_theme
ggplot2::ggsave(file.path(fd,"logistic_log10_Qj_R0_r_facets.png"),qplot,width=11,height=8,dpi=180)
ggplot2::ggsave(file.path(fd,"logistic_log10_Qj_R0_r_facets.pdf"),qplot,width=11,height=8)
selected <- do.call(rbind,lapply(c(1.5,2,2.5,3,4,5),function(R) {
  z<-logistic_multitrough_probabilities(R,.1,K,I0,5,dt=.02)
  if(!nrow(z$trough_table))return(NULL)
  transform(z$trough_table,R0_label=factor(R))
}))
lp <- ggplot2::ggplot(selected,ggplot2::aes(trough_index,P_persistence,
  colour=R0_label,group=R0_label))+ggplot2::geom_line()+
  ggplot2::geom_point(size=2)+ggplot2::scale_x_continuous(breaks=1:5)+
  ggplot2::scale_y_continuous(limits=c(0,1))+
  ggplot2::labs(x="Trough index j",y="Conditional persistence P_j",
    colour="R0",title="Conditional persistence across successive troughs",
    caption=paste0("Fixed r = 0.1. ",cap))+base_theme
ggplot2::ggsave(file.path(fd,"logistic_Pj_by_trough_selected_R0.png"),lp,width=9,height=7,dpi=180)
ggplot2::ggsave(file.path(fd,"logistic_Pj_by_trough_selected_R0.pdf"),lp,width=9,height=7)
counts <- aggregate(P_persistence~trough_index,d,function(x)sum(is.finite(x)))
print(counts); print(status)
