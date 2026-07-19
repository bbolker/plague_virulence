# Established Patch Occupancy

The raw metapopulation occupancy statistic classifies a patch as occupied
whenever `I > 0`. This includes newly introduced infections that disappear
quickly and therefore may not represent successful establishment after the
first local burnout.

For a forward persistence window `tau`, patch `i` is instead classified as
established at time `t` when

\[
I_i(s) > 0 \quad \text{for every recorded } s \in [t,t+\tau].
\]

The analysis uses `tau = 50` disease generations. This conservative window was
chosen to reduce classification of long transient infections as established.
Earlier sensitivity checks at `tau = 5, 10, 20, 50` changed most trajectories
only slightly, so retaining several window values added complexity without much
scientific information. The threshold remains a pragmatic definition rather
than a uniquely identifiable biological constant.

This definition uses future observations, so it is an offline diagnostic rather
than a state variable that could drive the simulation in real time. When the
full future window is unavailable near the end of a trajectory, established
occupancy is recorded as missing (`NA`) rather than inferred. Metapopulation
fractions exclude these missing values.

The odin transmission and demographic equations are unchanged. The new runner
uses the same scenarios, seeds, explicit virgin-soil initial conditions, odin
generator, and stochastic simulation logic as the existing occupancy analysis;
only the summary statistic and output location differ. Existing outputs under
`fadeout/output/patch_occupancy/` and
`fadeout/output/patch_occupancy_approximation/` remain untouched.

Run from the repository root:

```bash
Rscript fadeout/run_established_occupancy_exploration.R
Rscript fadeout/compare_established_occupancy_approximation.R
```

All new results are stored under
`fadeout/output/patch_occupancy_established/`. The first script saves compact
raw and established occupancy summaries, diagnostics, and raw-versus-established
figures. The second script performs the primary comparison of established
stochastic occupancy with the existing logistic analytical approximation and
writes to the nested `analytical_comparison/` directory. It generates one
`tau=50` trajectory overlay for every varied parameter; raw occupancy is not
used as the stochastic curve in these primary plots.

The primary outputs are the established-versus-analytical PDFs in
`analytical_comparison/figures/`. Figures are saved as PDF only; the plotted
curve values and diagnostics are also saved as CSV.

Here `tau` is measured in disease-generation time units, not stored rows. With
`dt = 0.1`, `tau = 50` uses 500 forward simulation steps. For comparison, the
separate seasonal single-strain raster uses days and
`gamma = 0.2/day`, so one mean infectious period is about five days: windows of
50 disease generations correspond approximately to 250 seasonal days. The
seasonal raster contains some long terminated infection bouts, motivating a
conservative window, but it does not by itself calibrate a threshold for the
non-seasonal model.

The statistic has important limitations. Its value depends on `tau`; repeated
immigration can keep `I` positive through a window even if infection would not
be self-sustaining locally; and rescue, rapid extinction followed by
reintroduction, or the discrete observation interval can change classification.
It therefore measures uninterrupted observed infection over a chosen window,
not a causal probability of permanent local establishment.
