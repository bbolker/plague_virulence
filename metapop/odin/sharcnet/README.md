# sharcnet — Alliance Canada (Nibi) job array scripts

Submit all jobs from this directory (`metapop/odin/sharcnet/`).  
Before the first submission: `mkdir -p logs outputs`

Need to install all the required packages (this takes a while):

```
module load r/4.5.0
module load gsl ## not sure about this, maybe required for some packages?
Rscript --vanilla get_pkgs.R  ## in metapop/ dir
```

The `plagueMetapop` package must be installed on the cluster:
```
R CMD INSTALL ../../plagueMetapop
```

---

## Jobs

### euler_onepatch_onestrain_extinct

One-patch, one-strain extinction grid: varies R0 × K.  
Four combinations of demography (logistic / linear restoring force) × dynamics (continuous Euler / Reed-Frost).

| file | purpose |
|------|---------|
| `submit_euler_extinct_logistic_continuous.sh` | logistic demography, continuous Euler (520 tasks, 5 min) |
| `submit_euler_extinct_logistic_reedfrost.sh` | logistic demography, Reed-Frost (520 tasks, 5 min) |
| `submit_euler_extinct_linear_continuous.sh` | linear demography, continuous Euler (520 tasks, 5 min) |
| `submit_euler_extinct_linear_reedfrost.sh` | linear demography, Reed-Frost (520 tasks, 5 min) |
| `submit_euler_extinct_mini.sh` | mini test (32 tasks, coarser grid) |
| `euler_onepatch_onestrain_extinct_run_array.R` | shared array script (`--lineargrowth`, `--reedfrost` flags) |
| `euler_onepatch_onestrain_extinct_run_array_mini.R` | mini array script (separate parameters) |
| `euler_onepatch_onestrain_extinct_combine.R` | combine outputs for one combination (takes combination name as argument) |
| `euler_onepatch_onestrain_extinct_mini_combine.R` | combine mini outputs |

Grid (full): `R0 = seq(1.1, 5, by=0.1)` × `K = 10^seq(3, 6, by=0.25)` — 40 × 13 = 520 tasks  
Grid (mini): `R0 = seq(1.1, 5, by=0.5)` × `K = 10^seq(3, 6)` — 32 tasks  
Parameters: `nsim=1000`, `n_patch=1`; `dt=0.1` (continuous), `dt=1` (Reed-Frost)

```bash
sbatch submit_euler_extinct_logistic_continuous.sh   # logistic demography, continuous Euler
sbatch submit_euler_extinct_logistic_reedfrost.sh    # logistic demography, Reed-Frost
sbatch submit_euler_extinct_linear_continuous.sh     # linear demography, continuous Euler
sbatch submit_euler_extinct_linear_reedfrost.sh      # linear demography, Reed-Frost

# after each run completes, combine with the matching combination name:
Rscript euler_onepatch_onestrain_extinct_combine.R logistic_continuous
Rscript euler_onepatch_onestrain_extinct_combine.R logistic_reedfrost
Rscript euler_onepatch_onestrain_extinct_combine.R linear_continuous
Rscript euler_onepatch_onestrain_extinct_combine.R linear_reedfrost
```

Outputs: `outputs/euler_onepatch_onestrain_extinct_{logistic|linear}_{continuous|reedfrost}_task_NNNN.rds`  
Combined: `outputs/euler_onepatch_onestrain_extinct_{logistic|linear}_{continuous|reedfrost}.rds`

---

### euler_onestrain

Multi-patch, one-strain grid: varies R0 × K × alpha.  
A single `_run_array.R` handles all modes; the submit script passes the appropriate flag.

| file | purpose |
|------|---------|
| `submit_euler_onestrain.sh` | full run (1400 tasks, 32G, 1 CPU) |
| `submit_euler_onestrain_mini.sh` | mini test (24 tasks, 2G, 4 CPUs) |
| `submit_euler_onestrain_batch2.sh` | batch 2: higher alpha range (840 tasks, 32G, 1 CPU) |
| `euler_onestrain_run_array.R` | shared array script (`--mini` / `--batch2` flags) |
| `euler_onestrain_combine.R` | combine outputs (`--mini` / `--batch2` flags) |

Grid (full):   `R0 = seq(1.1, 5, by=0.1)` × `K = 10^seq(3, 6, by=0.5)` × `alpha = 10^seq(-5.5, -3.5, by=0.5)` — 1400 tasks  
Grid (batch2): `R0 = seq(1.1, 5, by=0.1)` × `K = 10^seq(3, 6, by=0.5)` × `alpha = 10^seq(-3, -2, by=0.5)` — 840 tasks  
Grid (mini):   `R0 = seq(1.1, 3, by=0.5)` × `K = 10^seq(4, 6)` × `alpha = 10^seq(-5, -4)` — 24 tasks  
Parameters (full/batch2): `nsim=200`, `dt=0.1`, `n_patch=200`  
Parameters (mini): `nsim=50`, `dt=0.2`, `n_patch=50`

```bash
sbatch submit_euler_onestrain_mini.sh   # test first
Rscript euler_onestrain_combine.R --mini

sbatch submit_euler_onestrain.sh        # full run
Rscript euler_onestrain_combine.R

sbatch submit_euler_onestrain_batch2.sh # higher-alpha extension
Rscript euler_onestrain_combine.R --batch2
```

Outputs: `outputs/euler_onestrain[_mini|_batch2]_task_NNNN.rds`  
Combined: `euler_onestrain[_mini|_batch2].rds`

---

### euler_twostrain

Multi-patch, two-strain PIP grid: varies R01 × R02 × K × alpha.  
Strain 1 runs alone for `strain2_delay = round(100/dt)` steps before strain 2 is seeded.
Full and mini runs stop when either strain goes extinct (`stop_either_extinct()`);
mini2 stops only when both strains are simultaneously extinct (`stop_both_extinct`).

The full run uses **META-Farm** (dynamic load balancing) rather than a SLURM array;
the mini test still uses a plain SLURM array.

| file | purpose |
|------|---------|
| `meta_euler_twostrain/job_script.sh` | META metajob SLURM script (32G, 1 CPU, 8 h) |
| `meta_euler_twostrain/single_case.sh` | per-case worker: sets `SLURM_ARRAY_TASK_ID=$ID`, runs R script |
| `meta_euler_twostrain/make_table.sh` | generates `table.dat` (lines 1–14400) |
| `submit_euler_twostrain_mini.sh` | mini SLURM array (96 tasks, 4G, 4 CPUs) |
| `submit_euler_twostrain_mini2.sh` | mini2 SLURM array (100 tasks, 32G, 1 CPU) |
| `euler_twostrain_run_array.R` | shared R script (`--mini` / `--mini2` flags) |
| `euler_twostrain_combine.R` | combine outputs (`--mini` / `--mini2` flags, short: `-m` / `-2`) |

Grid (full):  `R01 = R02 = seq(1.1, 5, by=0.1)` × `K = 10^seq(3, 5)` × `alpha = 10^seq(-5, -3)` — 14400 tasks  
Grid (mini):  `R01 = R02 = seq(1.1, 3, by=0.5)` × `K = 10^seq(3, 5)` × `alpha = 10^seq(-5, -4)` — 96 tasks  
Grid (mini2): `R01 = R02 = seq(1.1, 4, length.out=10)` × `K = 1e4` × `alpha = 10^-5.5` — 100 tasks; uses `stop_both_extinct`  
Parameters (full/mini2): `nsim=200`, `dt=0.1`, `n_patch=200`  
Parameters (mini): `nsim=50`, `dt=0.2`, `n_patch=50`

**Wall-time derivation (from mini sacct, job 15051381):** Mini tasks took ~6.5 s elapsed with 4 CPUs → 0.52 CPU-s/sim. Full run: 4× more patches, 2× more steps, 200 sims on 1 CPU → ~4.2 CPU-s/sim × 200 = ~840 s ≈ 14 min/task. With 14.4 tasks/metajob → ~3.4 h worst case; 8 h provides ~2.4× buffer.

**Mini (SLURM array):**
```bash
sbatch submit_euler_twostrain_mini.sh
Rscript euler_twostrain_combine.R --mini
```

**Mini2 (SLURM array — fixed K=1e4, alpha=10^-5.5, full-run parameters):**
```bash
sbatch submit_euler_twostrain_mini2.sh
Rscript euler_twostrain_combine.R --mini2
```

**Full (META-Farm) — initialise once, then submit from `meta_euler_twostrain/`:**
```bash
module load meta-farm
farm_init.run meta_euler_twostrain   # run once from sharcnet/; fails if dir already exists
cd meta_euler_twostrain
bash make_table.sh                   # generates table.dat (14400 lines)
mkdir -p logs
submit.run 1000                      # 1000 metajobs × 8 h wall time
```

Note current meta run is 15056098_[1-1000]

To resubmit after failure without re-initialising, just re-run `submit.run` from inside `meta_euler_twostrain/`. To start completely fresh (new farm), run `clean.run` first — but note that `clean.run` deletes `make_table.sh` (a git-tracked file), so restore it afterwards:
```bash
clean.run
git checkout -- make_table.sh
bash make_table.sh
mkdir -p logs
submit.run 1000
```

Monitor progress (run from inside `meta_euler_twostrain/`):
```bash
Status.run    # note capital S — shows COMPLETED/RUNNING/FAILED counts
```

After completion (run from `sharcnet/`):
```bash
Rscript euler_twostrain_combine.R
```

Outputs: `outputs/euler_twostrain[_mini]_task_NNNNNN.rds`  
Combined: `euler_twostrain[_mini].rds`

---

## Utilities

```bash
bash check_jobs.sh      # state summary for most recent job array
bash check_jobs.sh 3    # state summaries for 3 most recent job arrays
sq                      # running/pending jobs
sacct -j <jobid> --format=JobID,JobName%20,CPUTime,Elapsed,State
sshare -l -u bolker -A def-bolker_cpu
```

