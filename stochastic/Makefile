outputs/%.rds:
	mkdir -p outputs
	R CMD BATCH --vanilla $<
	mv $(notdir $(basename $@))*.rds outputs

outputs/sim_batch0.rds: sim_batch0.R simulation_funs.R

outputs/sim_batch1.rds: sim_batch1.R simulation_funs.R

plot_batch1.pdf: outputs/sim_batch1.rds outputs/sim_batch0.rds simulation_funs.R plot_batch1.R
	R CMD BATCH --vanilla plot_batch1.R

plot_batch0.pdf: outputs/sim_batch0.rds simulation_funs.R plot_batch0.R
	R CMD BATCH --vanilla plot_batch0.R

plot_batch2.pdf: outputs/sim_batch1.rds outputs/sim_batch0.rds ../poisson/plot_alt.R
	../poisson/plot_alt.R --which 1 -i outputs/sim_batch0.rds -o $@

clean:
	rm *~ \#*
