R=R CMD BATCH --vanilla

outputs/sim.rds: ./sim.R simulation_funs.R
	$<

outputs/sim_fixed.rds: ./sim.R simulation_funs.R
	./sim.R -o $@ --fixedpop -r -n 2000 -c 26 >sim_fixed.Rout

outputs/poisson_plot_alt.pdf: outputs/sim.rds plot_alt.R
	./plot_alt.R -i outputs/sim.rds -o $@

outputs/poisson_plot_alt_fixed.pdf: outputs/sim_fixed.rds plot_alt.R
	./plot_alt.R -i outputs/sim_fixed.rds -o $@

## outputs/sim_fixed.ppa.pdf:
outputs/%.ppa.pdf: plot_alt.R outputs/%.rds
	./$< -i outputs/$*.rds -o $@

## for CSEE talk
outputs/poisson_ext.pdf: plot_alt.R
	./plot_alt.R -o outputs/poisson_ext.pdf --which 2 -w 8 --height 8

outputs/poisson_qepop.pdf: plot_alt.R
	./plot_alt.R -o outputs/poisson_qepop.pdf --which 3 -w 8 --height 8 --smooth line

%.png: %.pdf
	convert -density 300 $< $@


332:
	git commit -am "Meeting with Ben"
	git pull
	git push
