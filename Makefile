ALL: plague_conf.pdf

R=R CMD BATCH --vanilla

README.md: README.qmd
	quarto render $< 

plagueMetapop:
	Rscript -e "install.packages('plagueMetapop', repos = NULL)"

get_pkgs:
	$(R) get_pkgs.R

fastslow.pdf: fastslow.R
	R CMD BATCH --vanilla fastslow.R

FACTORIAL_FIGS = fadeout/output/occupancy_factorial_100y/figures
LATE_SUMMARY = $(FACTORIAL_FIGS)/occupancy_factorial_late_summary

## one run of the script writes all the factorial figures (PDF and PNG)
$(LATE_SUMMARY).pdf $(LATE_SUMMARY).png &: fadeout/seasonal/occupancy_factorial_100y_plotx.R fadeout/seasonal/seasonal_model_metapop.R
	Rscript $<

plague_conf.pdf: $(LATE_SUMMARY).pdf

## main.html: main.qmd virulence.bib
## main.pdf: main.qmd virulence.bib

%.pdf: %.qmd virulence.bib
	quarto render $< --to pdf

%.html: %.qmd virulence.bib
	quarto render $<

%.docx: %.qmd virulence.bib
	quarto render $< --to docx

.PRECIOUS: %.pdf
%.open: %.pdf
	open "$<"

# Optional features (add your own MK file, or use someone else's)
-include extras.mk
## jd.extras: jd.MK
%.extras: %.MK
	/bin/ln -fs $< extras.mk


