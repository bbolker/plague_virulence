R=R CMD BATCH --vanilla
DATADIR=validation/data
SIMFILES=$(DATADIR)/stochastic_scan_results.csv $(DATADIR)/stochastic_scan_analytic.csv $(DATADIR)/fig11_scan_analytic.csv

validation/figures/fig13_comb.pdf: validation/plot_fig13_comb.R $(SIMFILES)
	$(R) $<
validation/figures/fig13_BI_conditional_stochastic_validation.pdf: validation/plot_BI_validation.R $(SIMFILES)
	$(R) $<

$(DATADIR)/stochastic_scan_results.csv: validation/run_stochastic_scan.R
	mkdir -p $(DATADIR)
	$(R) $<

$(DATADIR)/stochastic_scan_analytic.csv: validation/plot_stochastic_scan.R
	$(R) $<

$(DATADIR)/fig11_scan_results.csv: validation/run_fig11_scan.R
	$(R) $<

$(DATADIR)/fig11_scan_analytic.csv: validation/plot_fig11_scan.R validation/R/scan.R $(DATADIR)/fig11_scan_results.csv
	$(R) $<

# looking for fig13_BI_conditional_stochastic_validation.pdf, eventually
