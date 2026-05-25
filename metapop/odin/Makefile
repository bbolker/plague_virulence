## This is plague_virulence/metapop/odin

current: target
-include target.mk
Ignore = target.mk

vim_session:
	bash -ic "vmt"

######################################################################

compare_platforms.Rout: compare_platforms.R
	$(pipeR)

tests.Rout: tinytest/tests.R
	$(pipeR)

discrete_onepatch_twostrain_extinct.png: odin_twostrain_run1.R
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
