outputs/%.rds:
	mkdir -p outputs
	R CMD BATCH --vanilla $<
	mv $(notdir $(basename $@))*.rds outputs

outputs/sim_batch0.rds: sim_batch0.R simulation_funs.R

outputs/sim_batch1.rds: sim_batch1.R simulation_funs.R

plot_batch1.pdf: outputs/sim_batch1.rds outputs/sim_batch0.rds simulation_funs.R plot_batch1.R
	R CMD BATCH --vanilla plot_batch1.R

clean:
	rm *~ \#*
