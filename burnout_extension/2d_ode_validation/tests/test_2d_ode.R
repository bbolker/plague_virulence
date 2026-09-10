source('2d_ode_validation/R/trough_kendall_psi.R')
z<-trough_kendall_2d(1.5,.05,0,.5,1000)
stopifnot(z$status=='OK',z$t_minus<z$t_t,z$t_t<z$t_plus,z$alpha_t>0,
  z$L_minus>0,z$L_plus>0,z$P_trough>0,z$P_trough<1,
  z$P_kendall>0,z$P_kendall<1)
g_residual<-1.5*z$x_t*(z$x_t+z$y_t)^(-.5)-1
stopifnot(abs(g_residual)<1e-7,
  abs(z$y_t-1.417e-2)<4e-4,
  abs(z$alpha_t-.01452)<5e-4,
  abs(z$P_kendall-.4639)<.006,
  abs(z$P_trough-.4941)<.006)
cat('2D ODE trough/Kendall regression tests passed\n')
