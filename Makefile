
fastslow.pdf: fastslow.R
	R CMD BATCH --vanilla fastslow.R

%.html: %.md virulence.bib
	quarto render $<

## abstract2.html: abstract2.qmd virulence.bib
%.html: %.qmd virulence.bib
	quarto render $<

%.docx: %.qmd virulence.bib
	quarto render $< --to docx

