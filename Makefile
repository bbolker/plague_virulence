
fastslow.pdf: fastslow.R
	R CMD BATCH --vanilla fastslow.R

## main.html: main.qmd virulence.bib
## main.pdf: main.qmd virulence.bib

%.pdf: %.qmd virulence.bib
	quarto render $< --to pdf

%.html: %.qmd virulence.bib
	quarto render $<

%.docx: %.qmd virulence.bib
	quarto render $< --to docx
