R=R CMD BATCH --vanilla

validation/data/stochastic_scan_results.csv: validation/run_stochastic_scan.R
	mkdir -p validation/data
	$(R) $<
