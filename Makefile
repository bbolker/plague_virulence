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

## one run of the script writes the summary data and all the factorial
## figures (PDF and PNG)
$(LATE_SUMMARY).pdf $(LATE_SUMMARY).png &: fadeout/seasonal/occupancy_factorial_100y.R fadeout/seasonal/seasonal_model_metapop.R
	Rscript $<

plague_conf.pdf: $(LATE_SUMMARY).pdf

BATCH_DATA = fadeout/output/occupancy_factorial_100y_batch/data
BATCH_FIGS = fadeout/output/occupancy_factorial_100y_batch/figures
BATCH_HEATMAPS = $(foreach h,occupancy persistence burnout fadeout_hazard amp_burnout amp_fadeout_hazard,\
	$(foreach ext,pdf png,$(BATCH_FIGS)/heatmap_$(h).$(ext)))

## the batch simulations (9800 runs, ~2 h on 15 workers / 8 physical cores;
## set N_WORKERS to override) write all the summary data
$(BATCH_DATA)/combo_summary.csv $(BATCH_DATA)/run_summary.csv $(BATCH_DATA)/settings.csv &: fadeout/seasonal/occupancy_factorial_100y_batch.R fadeout/seasonal/seasonal_model_metapop.R fadeout/seasonal/seasonal_fadeout_functions.R
	Rscript $<

## the heatmaps are drawn from the summary data in a few seconds
$(BATCH_HEATMAPS) &: fadeout/seasonal/occupancy_factorial_100y_batch_plots.R $(BATCH_DATA)/combo_summary.csv $(BATCH_DATA)/run_summary.csv $(BATCH_DATA)/settings.csv
	Rscript $<

batch_heatmaps: $(BATCH_HEATMAPS)

PIP_DATA = fadeout/output/pairwise_invasion_100y/data
PIP_FIGS = fadeout/output/pairwise_invasion_100y/figures
PIP_PLOTS = $(foreach h,pip_outcomes invader_survival occupancy_series,\
	$(foreach ext,pdf png,$(PIP_FIGS)/$(h).$(ext)))

## pairwise invasion runs at K = 5000, amplitude 0.2, alpha = 5e-5
## (3600 runs, ~55 min on 15 workers; set N_WORKERS to override)
$(PIP_DATA)/combo_summary.csv $(PIP_DATA)/run_summary.csv $(PIP_DATA)/occupancy_thinned.rds $(PIP_DATA)/settings.csv &: fadeout/seasonal/pairwise_invasion_100y.R fadeout/seasonal/seasonal_model_metapop.R
	Rscript $<

$(PIP_PLOTS) &: fadeout/seasonal/pairwise_invasion_100y_plots.R $(PIP_DATA)/combo_summary.csv $(PIP_DATA)/run_summary.csv $(PIP_DATA)/occupancy_thinned.rds $(PIP_DATA)/settings.csv
	Rscript $<

pip_plots: $(PIP_PLOTS)

## main.html: main.qmd virulence.bib
## main.pdf: main.qmd virulence.bib

%.pdf: %.qmd virulence.bib
	quarto render $< --to pdf

%.html: %.qmd virulence.bib
	quarto render $<

%.docx: %.qmd virulence.bib
	quarto render $< --to docx

## GitHub-flavored markdown, for documents meant to be read on GitHub
%.md: %.qmd virulence.bib
	quarto render $< --to gfm

.PRECIOUS: %.pdf
%.open: %.pdf
	open "$<"

# Optional features (add your own MK file, or use someone else's)
-include extras.mk
## jd.extras: jd.MK
%.extras: %.MK
	/bin/ln -fs $< extras.mk


