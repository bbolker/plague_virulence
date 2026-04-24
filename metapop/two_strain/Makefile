
default: ode_test_proc.Rout
R=R CMD BATCH --vanilla

outputs/%.rds:
	mkdir -p outputs
	$(R) $<
	mv $(notdir $(basename $@))*.rds outputs

outputs_2strain/finalsize.rds: ode_test.R
	$(R) $<

ode_test_proc.Rout: ode_test_proc.R
	$(R) $<

two_strain.pdf: two_strain.tex
	texi2dvi -p $<

sim_fun.Rout: sim_fun.R
	$(R) $<

outputs_pip/sim_pip_alpha_rho.rds: sim_pip.R sim_fun.Rout
	$(R) $<

pkgs: get_pkgs.R
	$(R) $<

clean:
	rm *~ \#*
