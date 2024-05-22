fastslow.pdf: fastslow.R
	R CMD BATCH --vanilla fastslow.R

%.html: %.md virulence.bib
	quarto render $<

%.html: %.qmd virulence.bib
	quarto render $<
