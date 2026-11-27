R=R CMD BATCH --vanilla
DATADIR=validation/data
SIMFILES=$(DATADIR)/stochastic_scan_results.csv $(DATADIR)/stochastic_scan_analytic.csv $(DATADIR)/fig11_scan_analytic.csv
ALL: paper

.PHONY: paper-figures paper landscape-figures

paper-figures: $(SIMFILES) $(DATADIR)/fig11_scan_results.csv
	Rscript validation/make_section7_figures.R

paper: burnout_theta_paper.pdf

## Figures actually included by the paper (commented-out \includegraphics lines
## are skipped). Figures with rules below are rebuilt when stale; the rest are
## used as they are.
PAPER_FIGS := $(shell sed -n 's/^[^%]*\\includegraphics\(\[[^]]*\]\)\{0,1\}{\([^}]*\)}.*/\2/p' burnout_theta_paper.tex)

## BibTeX manages the independent reference database; unsrt uses citation order.
burnout_theta_paper.pdf: burnout_theta_paper.tex references.bib $(PAPER_FIGS)
	pdflatex -interaction=nonstopmode -halt-on-error $<
	bibtex burnout_theta_paper
	pdflatex -interaction=nonstopmode -halt-on-error $<
	pdflatex -interaction=nonstopmode -halt-on-error $<

######################################################################
## Analytic persistence landscapes (no stochastic simulation needed).
## Each script takes well under a minute from empty caches.

PAPERFIG=validation/figures/paper
LANDSCAPE_LIB=validation/R/theory.R validation/R/kendall.R validation/R/extrema.R
COND_CACHE=$(DATADIR)/paper_conditional_landscape_theta5_R32.csv
EXT_CACHE=$(DATADIR)/unconditional_landscape_extrema_theta5_K2_R32.csv
STRUCT_CACHES=$(DATADIR)/paper_unconditional_landscape.csv \
              $(DATADIR)/paper_extrema_branches.csv $(DATADIR)/paper_theta_critical.csv

## The scripts reuse these caches whenever they exist (the extremum cache has
## no version check), so drop any cache older than the analytic library code.
define drop_stale_caches
	@for c in $(1); do for l in $(LANDSCAPE_LIB); do \
	  if [ -e $$c ] && [ $$l -nt $$c ]; then echo "removing stale $$c"; rm -f $$c; fi; \
	done; done
endef

landscape-figures: $(PAPERFIG)/conditional_persistence_landscape.pdf \
                   $(PAPERFIG)/unconditional_persistence_structure.pdf \
                   $(PAPERFIG)/unconditional_persistence_landscape.pdf

## One run of the section 7 script writes both PDFs (grouped target, make >= 4.3).
$(PAPERFIG)/conditional_persistence_landscape.pdf $(PAPERFIG)/unconditional_persistence_structure.pdf &: \
    validation/make_section7_figures.R $(LANDSCAPE_LIB)
	mkdir -p $(DATADIR)
	$(call drop_stale_caches,$(COND_CACHE) $(STRUCT_CACHES))
	Rscript validation/make_section7_figures.R --analytic-only

## Both scripts write COND_CACHE; the order-only prerequisite keeps them from
## running concurrently under make -j (the second one then reuses the cache).
$(PAPERFIG)/unconditional_persistence_landscape.pdf: validation/make_unconditional_landscape_figure.R $(LANDSCAPE_LIB) \
    | $(PAPERFIG)/conditional_persistence_landscape.pdf
	mkdir -p $(DATADIR)
	$(call drop_stale_caches,$(COND_CACHE) $(EXT_CACHE))
	Rscript validation/make_unconditional_landscape_figure.R --tikz

######################################################################

validation/figures/fig13_comb.pdf: validation/plot_fig13_comb.R $(SIMFILES)
	$(R) $<
validation/figures/fig13_BI_conditional_stochastic_validation.pdf: validation/plot_BI_validation.R $(SIMFILES)
	$(R) $<

$(DATADIR)/stochastic_scan_results.csv: validation/run_stochastic_scan.R
	mkdir -p $(DATADIR)
	$(R) $<

######################################################################

## Making these csvs is expensive
$(DATADIR)/stochastic_scan_analytic.csv: validation/plot_stochastic_scan.R
	$(R) $<

$(DATADIR)/fig11_scan_results.csv: validation/run_fig11_scan.R
	$(R) $<

$(DATADIR)/fig11_scan_analytic.csv: validation/plot_fig11_scan.R validation/R/scan.R $(DATADIR)/fig11_scan_results.csv
	$(R) $<

# looking for fig13_BI_conditional_stochastic_validation.pdf, eventually

## Final figures with LaTeX typography. Requires tikzDevice and XeLaTeX;
## TIKZ_XELATEX and TIKZ_TECTONIC can override the compiler paths.
.PHONY: paper-tikz-figures
paper-tikz-figures:
	Rscript validation/plot_phaseplane_parsons_style.R
	Rscript validation/make_paper_figures.R --fig1-only
	Rscript validation/make_section7_figures.R --conditional-landscape-only --tikz
	Rscript validation/make_unconditional_validation_figure.R
	Rscript validation/make_unconditional_landscape_figure.R --tikz
	Rscript validation/plot_entry_curve_comparison.R --plot-only --paper-export
