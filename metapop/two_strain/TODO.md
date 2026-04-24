## 24 April

* understand full pipeline
* re-run with more sims (why aren't curves smooth)?
* modify PIPs with more realistic {I10, I20, R01, R02} → {final size, fraction by strain} mapping
   * is it worth truncating coefficients?
   * better heuristics? (JD, DJDE?)
   * run on SHARCnet?? (via job arrays: factorial design is size 3600. MaxArraySize/MaxSubmit are both 10,000 on nibi (chunks of size 100 probably make sense. Design is 20 (R01) x 20 (R02) x 3 (alpha) x 3 (rho). nsim = 200? 2000 ?
* assume initial numbers are small ? (yes: computing final size based on classic formula)
* understanding invasion: benefit of lower R0 will depend on
   * probability of patch coinfection (unchanged to first order, but determines balance of factors below)
   * change in burnout probability (decrease [if near burnout max ~ R0=2.0], good)
   * change in fizzle probability (increase, bad)
   * change in final size in singly infected patches (decrease, bad)
   * change in share of co-infected patches (decrease, bad)
 Can we find a sensible way to combine these?
* check PIP summaries again before running ...
