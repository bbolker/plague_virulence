# sharcnet — Alliance Canada (Nibi) job array scripts

Submit all jobs from this directory (`metapop/odin/sharcnet/`).  
Before the first submission: `mkdir -p logs outputs`

The `plagueMetapop` package must be installed on the cluster:
```
R CMD INSTALL ../../plagueMetapop
```

---

## Jobs

### euler_onepatch_onestrain_extinct

One-patch, one-strain extinction grid: varies R0 × K.

| file | purpose |
|------|---------|
| `submit_euler_extinct.sh` | full run (520 tasks, R0 × K grid) |
| `submit_euler_extinct_mini.sh` | mini test (32 tasks, coarser grid) |
| `euler_onepatch_onestrain_extinct_run_array.R` | shared array script (full and mini) |
| `euler_onepatch_onestrain_extinct_run_array_mini.R` | mini array script (separate parameters) |
| `euler_onepatch_onestrain_extinct_combine.R` | combine full outputs |
| `euler_onepatch_onestrain_extinct_mini_combine.R` | combine mini outputs |

Grid (full): `R0 = seq(1.1, 5, by=0.1)` × `K = 10^seq(3, 6, by=0.25)` — 40 × 13 = 520 tasks  
Grid (mini): `R0 = seq(1.1, 5, by=0.5)` × `K = 10^seq(3, 6)` — 32 tasks  
Parameters: `nsim=1000`, `dt=0.1`, `n_patch=1`

```bash
sbatch submit_euler_extinct.sh          # or submit_euler_extinct_mini.sh
Rscript euler_onepatch_onestrain_extinct_combine.R      # after completion
```

Outputs: `outputs/euler_onepatch_onestrain_extinct_task_NNNN.rds`  
Combined: `outputs/euler_onepatch_onestrain_extinct[_mini].rds`

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
Strain 1 runs alone for `strain2_delay = round(100/dt)` steps before strain 2 is seeded;
run stops when either strain goes extinct (`stop_either_extinct()`).

| file | purpose |
|------|---------|
| `submit_euler_twostrain.sh` | full run (14400 tasks, 32G, 1 CPU) |
| `submit_euler_twostrain_mini.sh` | mini test (96 tasks, 4G, 4 CPUs) |
| `euler_twostrain_run_array.R` | shared array script (`--mini` flag) |
| `euler_twostrain_combine.R` | combine outputs (`--mini` flag) |

Grid (full): `R01 = R02 = seq(1.1, 5, by=0.1)` × `K = 10^seq(3, 5)` × `alpha = 10^seq(-5, -3)` — 14400 tasks  
Grid (mini): `R01 = R02 = seq(1.1, 3, by=0.5)` × `K = 10^seq(3, 5)` × `alpha = 10^seq(-5, -4)` — 96 tasks  
Parameters (full): `nsim=200`, `dt=0.1`, `n_patch=200`  
Parameters (mini): `nsim=50`, `dt=0.2`, `n_patch=50`

```bash
sbatch submit_euler_twostrain_mini.sh   # test first
Rscript euler_twostrain_combine.R --mini

sbatch submit_euler_twostrain.sh        # full run
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
```
