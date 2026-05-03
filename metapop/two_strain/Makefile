
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

outputs_pip/sim_pip_alpha_rho.rds: sim_pip.R simulation_funs.R
	./sim_pip.R --cores 26 --nsim 2000 >sim_pip.Rout

test:
	Rscript --vanilla -e "tinytest::run_test_file('sim_tests.R')"

pkgs: get_pkgs.R
	$(R) $<

clean:
	rm -f *~ \#* *.aux *.log Rplots.pdf

