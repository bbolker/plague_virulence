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
Parameters: `nsim=500`, `dt=0.1`, `n_patch=1`

```bash
sbatch submit_euler_extinct.sh          # or submit_euler_extinct_mini.sh
Rscript euler_onepatch_onestrain_extinct_combine.R      # after completion
```

Outputs: `outputs/euler_onepatch_onestrain_extinct_task_NNNN.rds`  
Combined: `outputs/euler_onepatch_onestrain_extinct[_mini].rds`

---

### euler_onestrain

Multi-patch, one-strain grid: varies R0 × K × alpha.  
Uses a single `_run_array.R`; the submit script passes `--mini` to switch modes.

| file | purpose |
|------|---------|
| `submit_euler_onestrain.sh` | full run (1400 tasks) |
| `submit_euler_onestrain_mini.sh` | mini test (24 tasks) |
| `euler_onestrain_run_array.R` | shared array script (full and mini via `--mini`) |
| `euler_onestrain_combine.R` | combine outputs (`--mini` for mini) |

Grid (full): `R0 = seq(1.1, 5, by=0.1)` × `K = 10^seq(3, 6, by=0.5)` × `alpha = 10^seq(-5.5, -3.5, by=0.5)` — 1400 tasks  
Grid (mini): `R0 = seq(1.1, 3, by=0.5)` × `K = 10^seq(4, 6)` × `alpha = 10^seq(-5, -4)` — 24 tasks  
Parameters (full): `nsim=200`, `dt=0.1`, `n_patch=200`  
Parameters (mini): `nsim=50`, `dt=0.2`, `n_patch=50`

```bash
sbatch submit_euler_onestrain_mini.sh   # test first
Rscript euler_onestrain_combine.R --mini

sbatch submit_euler_onestrain.sh        # full run
Rscript euler_onestrain_combine.R
```

Outputs: `outputs/euler_onestrain[_mini]_task_NNNN.rds`  
Combined: `euler_onestrain[_mini].rds`

---

## Checking job status

```bash
sq                              # running/pending jobs
sacct -j <jobid> --format=JobID,JobName%20,CPUTime,Elapsed,State
```
