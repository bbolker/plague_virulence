## This is a PERSONAL MAKEFILE for plague_virulence
## It is version-controlled as jd.MK and linked and used as extras.mk
## `make jd.extras` to use

current: target
-include target.mk

vim_session:
	bash -ic "vi extras.mk target.mk Makefile"

Sources += $(wildcard *.MK) Makefile
Ignore += extras.mk 

######################################################################

Ignore += *.pdf
Ignore += *.docx
Ignore += *.html
Ignore += *_files/

######################################################################

## Inherited from wrap

## metapop/README.md

Sources += $(wildcard metapop/*.R metapop/*.md)

burnout5.Rout: metapop/burnout5.R
	$(pipeR)

## metapop/stochastic/

## Hastily because Ben is here
## When opportunity arises, make a Makefile in metapop/stochastic
## metapop/stochastic/Rplots.pdf: 
## See alsometapop/stochastic/sim_batch1.R
metapop/stochastic/sim_batch0.rds: metapop/stochastic/sim_batch0.R
	cd $(dir $<) && R --vanilla < $(notdir $<)

######################################################################

plagueMetapop.pkg: %.pkg: %/NAMESPACE
	R CMD INSTALL $*
	$(touch)

######################################################################

autopipeR = defined

## fadeout/single_strain_metapop.Rout: fadeout/single_strain_metapop.R
## fadeout/seasonal_single_strain_metapop.Rout: fadeout/seasonal_single_strain_metapop.R

######################################################################

### Makestuff

Ignore += makestuff

msrepo = https://github.com/dushoff
Makefile: makestuff/01.stamp
makestuff/%.stamp: | makestuff
	- $(RM) makestuff/*.stamp
	cd makestuff && $(MAKE) pull
	touch $@
makestuff:
	git clone --depth 1 $(msrepo)/makestuff

-include makestuff/os.mk

-include makestuff/texj.mk
-include makestuff/pipeR.mk

-include makestuff/git.mk
-include makestuff/visual.mk



