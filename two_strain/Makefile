
default: ode_test_proc.Rout
R=R CMD BATCH --vanilla

outputs/%.rds:
	mkdir -p outputs
	$(R) $<
	mv $(notdir $(basename $@))*.rds outputs

outputs_2strain/plot_2strain_points.pdf:
	$(R) plot_2strain.R

outputs_2strain/finalsize.rds: ode_test.R
	$(R) $<

ode_test_proc.Rout: ode_test_proc.R
	$(R) $<

polyfit.rda: ode_test_proc.R
	$(R) $<

two_strain_test.pdf: two_strain_test.R simulation_funs.R
	$(R) $<

two_strain.pdf: two_strain.tex
	texi2dvi -p $<

simulation_funs.Rout: simulation_funs.R
	$(R) $<

sim_pip.Rout: sim_pip.R simulation_funs.R
	./sim_pip.R --cores 26 --nsim 200 >sim_pip.Rout

outputs_pip/sim_pip_alpha_rho.rds: sim_pip.Rout
	./assemble_pip.R >assemble_pip.Rout

outputs_pip/pip_hybrid.pdf: outputs_pip/sim_pip_alpha_rho.rds plot_pip.R
	./plot_pip.R >plot_pip.Rout

pkgs: get_pkgs.R
	$(R) $<

clean:
	rm -f *~ \#* *.aux *.log Rplots.pdf

