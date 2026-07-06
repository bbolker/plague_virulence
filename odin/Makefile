## This is plague_virulence/odin

current: target
-include target.mk
Ignore = target.mk

vim_session:
	bash -ic "vmt"

R=R CMD BATCH --vanilla
OUTDIR=sharcnet/outputs

######################################################################

euler_twostrain_singlepatchintro_examples_plot.Rout: euler_twostrain_singlepatchintro_examples_plot.R outputs/euler_twostrain_singlepatchintro_examples.rds
	$(pipeR)

outputs/euler_twostrain_singlepatchintro_examples.rds: euler_twostrain_singlepatchintro_examples.R
	$(pipeR)

compare_platforms.Rout: compare_platforms.R
	$(pipeR)

ode_trough_linear.pdf: ode_trough_linear.Rout ;

ode_trough.rds:	oud_trough_run.Rout 
	$(pipeR)
ode_trough_linear.Rout: ode_trough_plot.R outputs/ode_trough.rds
	$(pipeR)

tests.Rout: tinytest/tests.R
	$(pipeR)

discrete_onepatch_twostrain_extinct.rds: discrete_onepatch_twostrain_extinct_run.R
	$(pipeR)

euler_onepatch_onestrain_extinct.rds: euler_onepatch_onestrain_extinct_run.R
	$(R) $<

euler_onepatch_onestrain_extinct_%.png: euler_onepatch_onestrain_extinct_plot.R $(OUTDIR)/euler_onepatch_onestrain_extinct_%.rds
	Rscript $< --combo $*

euler_onepatch_onestrain_extinct_comb.png: euler_onepatch_onestrain_comb_plot.R $(OUTDIR)/euler_onepatch_onestrain_extinct_logistic_continuous.rds
	Rscript $< --combo $*

euler_onestrain_plot.png: euler_onestrain_plot.R $(OUTDIR)/euler_onestrain.rds
	Rscript $< --output $@

euler_onepatch_onestrain_example.pdf: euler_onepatch_onestrain_example.R
	$(R) $<

discrete_onepatch_twostrain_extinct.png: discrete_onepatch_twostrain_extinct_plot.R discrete_onepatch_twostrain_extinct.rds
	$(pipeR)

## not working yet
tests:
	R --vanilla -e "tinytest::run_test_dir('tinytest')"

######################################################################

### Makestuff

Sources += Makefile

Ignore += makestuff
msrepo = https://github.com/dushoff

## ln -s ../makestuff . ## Do this first if you want a linked makestuff
Makefile: makestuff/00.stamp
makestuff/%.stamp: | makestuff
	- $(RM) makestuff/*.stamp
	cd makestuff && $(MAKE) pull
	touch $@
makestuff:
	git clone --depth 1 $(msrepo)/makestuff

-include makestuff/os.mk

-include makestuff/pipeR.mk

-include makestuff/git.mk
-include makestuff/visual.mk
