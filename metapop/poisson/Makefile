R=R CMD BATCH --vanilla

outputs/sim.rds: ./sim.R simulation_funs.R
	$<

outputs/sim_fixed.rds: ./sim.R simulation_funs.R
	./sim.R -o $@ --fixedpop -r -n 2000 -c 26 >sim_fixed.Rout
