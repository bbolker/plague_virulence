# plague virulence metapopulation model


<!-- don't edit the .md file directly! make it from README.qmd -->

## overview

This repo started with a collection of ideas about how to *verbally*
explain the evolution of attenuated virulence noted by Sidhu et al.
(2025), to a series of models trying to instantiate and explore those
verbal models; conceptually, this builds on Parsons et al. (2024),
although we have moved some distance from those analytical models. We
started with analytical (patch-occupancy, Levins-style) models that used
final-size and burnout equations to model demographic changes and
probabilities of extinction with a time step of a single (synchronized!)
epidemic “season”; a “2-D” version of this considered only a 2D state
space for the metapopulation (proportion of patches infected, average
patch population size); the 3-D version tracked proportion infected as
well as average pop size for both susceptible and infectious patches)
model. We could do more analysis on the 2D (`/fixedns` subdir) than the
3D (`/3d` subdir) model, but these models convinced us that there were
at least some plausible parameter ranges (`/parameters`) where a reduced
within-patch $R_0$ could lead to higher patch occupancy (these were all
single-strain models).

From newest to oldest:

- `odin`: mostly continuous-time simulations with
  [odin](https://mrc-ide.github.io/odin/) (and also comparable
  implementations with [macpan2](https://canmod.github.io/macpan2/). The
  `sharcnet` subdirectory has machinery for submitting large jobs as
  SLURM job arrays (and reassembling the results when finished).
- `plagueMetapop`: mini-R package implementing the functionality used in
  the `odin` runs
- `talks`: BMB talks for Canadian Society of Ecology and Evolution, May
  2026, and Statistical Society of Canada, June 2026
- `parameters`: discussion of estimating/guessing orders of magnitude
  for parameters
- `poisson`: a single-strain, stochastic model driven by rat
  colonization limitation (colonization uses Poisson deviates)
- `two_strain`: a two-strain version of the Poisson model (PIP
  computations)
- `stochastic`: stochastic model with burnout
- `two_strain_burnout_probability`: calculations of burnout probability
  in two-strain (coinfected) patches
- `3d`: Levins model, state variables {average host pop size in S
  patches; average host pop size in I patches; fraction of patches
  infected}
- `within_season_transmission`: ?
- `report_for_summer`: YZ’s report
- `fixedns`: : 2D model Levins model, state variables: {average host pop
  size; fraction of patches infected}
- `averaging_all_patches`: ?
- `outputs`
- `abstract2.qmd`, `main.qmd`, `notes/notes_*` are from our early
  discussions of the verbal model

## references

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-Pars+2024" class="csl-entry">

Parsons, Todd L., Benjamin M. Bolker, Jonathan Dushoff, and David J. D.
Earn. 2024. “The Probability of Epidemic Burnout in the Stochastic SIR
Model with Vital Dynamics.” *Proceedings of the National Academy of
Sciences* 121 (5): e2313708120.
<https://doi.org/10.1073/pnas.2313708120>.

</div>

<div id="ref-sidhuAttenuationVirulenceYersinia2025" class="csl-entry">

Sidhu, Ravneet Kaur, Guillem Mas Fiol, Pierre Lê-Bury, Christian E.
Demeure, Emelyne Bougit, Rémi Beau, Charlotte Balière, et al. 2025.
“Attenuation of Virulence in *Yersinia Pestis* Across Three Plague
Pandemics.” *Science* 388 (6750): eadt3880.
<https://doi.org/10.1126/science.adt3880>.

</div>

</div>
