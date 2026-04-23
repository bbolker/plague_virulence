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

clean:
	rm *~ \#*
