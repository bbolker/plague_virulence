fastslow.pdf: fastslow.R
	R CMD BATCH --vanilla fastslow.R

abstract.html: abstract.md virulence.bib
	quarto render $<
