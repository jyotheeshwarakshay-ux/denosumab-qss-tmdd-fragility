# Estimation Fragility in Multi-Start Profile-Likelihood Assessment of Kss in a QSS-TMDD Model of Denosumab

Code and data to reproduce the analysis in "Estimation Fragility in
Multi-Start Profile-Likelihood Assessment of Kss in a QSS-TMDD Model of
Denosumab." Simulation, multi-start SAEM profiling, ground-truth recovery,
and reproducibility diagnostics for a mass-balance-corrected QSS-TMDD model
of denosumab, in R (nlmixr2 / rxode2).

---

## What this reproduces

This analysis examines the practical identifiability of the quasi-steady-state
constant (Kss) in a QSS-TMDD model of denosumab, using profile likelihood with
multi-start SAEM estimation on a mass-balance-corrected model. The principal
finding is that single-fit SAEM profile-likelihood estimation is fragile in
this setting: repeated fits reach different objective-function-value (OFV)
solution regions depending on the random seed, and the resulting profile shape
is not reproducible across two independent simulated N=30 datasets. This
disagreement was not attributable to the reproducibility dataset having
received more optimization effort, since equalizing both datasets' fit counts
(giving the pilot dataset the same additional starts) did not resolve it.
Increasing the SAEM run length also did not resolve the underlying
solution-region variability, and moved some fits away from the low-OFV region
rather than toward it, arguing against insufficient convergence as an
explanation. A further forensic experiment showed that even an identical,
explicit random seed, repeated across independent freshly started processes
under a verified-identical computational environment, produced different
discrete outcomes; this non-determinism persisted when R's own random-number
generator was also explicitly seeded, and the specific computational mechanism
was not localized.

Beyond this fragility, a direct ground-truth check found a more fundamental
limitation: on the reproducibility dataset, the objective-function value of
the true, known parameter vector (evaluated directly, without estimation) was
substantially worse than that of the lowest-OFV solution the multi-start
estimation actually found — meaning the estimation's own selection criterion
would not have picked the true parameters, and this held even though the
evaluation itself was verified deterministic. This indicates that, in this
configuration, objective-function value does not track proximity to the true
parameters, a concern more fundamental than fragility because it would persist
even for a perfectly reproducible estimator. This finding was demonstrated on
the reproducibility dataset; an analogous check on the independent pilot
dataset was directionally consistent but not a clean replication (one
subject's evaluation failed numerically and was excluded), so it is reported
as suggestive rather than as independent confirmation. The recovery failure
extended beyond Kss to most of the other estimated structural parameters.

This behavior did not straightforwardly generalize to a larger (N=80) dataset
— profile fits there converged to a substantially higher-OFV regime, and
warm-starting from the N=30 low-OFV solution region failed to recover it —
though because the N=80 dataset was an independent realization rather than a
resampling of the N=30 data, this analysis cannot distinguish sample-size
effects from dataset-realization variability. Kss is only weakly constrained
by pharmacokinetic-only data, consistent with the shrinkage-based limitation
reported by Choi et al. Robust profiling therefore requires multi-start
estimation with explicit solution-region and convergence assessment, together
with direct reproducibility checks — and even reproducibility checks alone
would not have been sufficient to establish that a result is near the true
parameters.

---

## Requirements

- R 4.6.1
- rxode2 5.1.5
- nlmixr2 7.0.1
- dplyr 1.2.1
- ggplot2 4.0.3 (used only by the corrected model script and the figure script)

---

## How to reproduce (run order)

The pipeline has two independently-runnable stages (data generation) followed
by analysis stages, several of which are themselves multi-driver sequences
that must run in order because they append to a shared results file. Running
steps out of order, or skipping a driver in a sequence, produces an incomplete
or wrong result — this is not a stylistic preference, it is how the crash-safe
append-and-skip pattern used throughout this pipeline works. Steps 7 and 8
(ground-truth recovery and OFV-at-truth) are independent of steps 2–6 and of
each other; they only need the raw datasets from step 1.

### 1. Generate the three datasets (any order, independent)
```
Rscript generate_corrected_full_data.R    # -> simulated_denosumab_data_v2_corrected_full.csv   (N=80)
Rscript generate_corrected_pilot_data.R   # -> simulated_denosumab_data_v2_corrected_pilot.csv   (N=30, draw 1)
Rscript generate_corrected_repro_data.R   # -> simulated_denosumab_data_v2_corrected_repro.csv   (N=30, draw 2)
```

### 2. N=80 profile
```
./run_multistart_corrected_full.sh
```
Drives `multistart_profile_Kss_corrected_full.R <grid_point 1-7> <seed>` over a
7-point grid × 5 seeds (35 fits) → `multistart_results_corrected_full.csv`.

### 3. Pilot profile — TWO drivers, in sequence
```
./run_multistart_corrected_pilot.sh          # seeds 42, 123, 999  (21 fits)
./run_multistart_corrected_pilot_topup.sh    # seeds 7, 2024       (14 fits)
```
Both drive `multistart_profile_Kss_corrected_pilot.R <grid_point> <seed>` and
append-and-skip to the same `multistart_results_corrected_pilot.csv`. **Both
must run** — the first alone produces only 21 of the 35 reported fits.

### 3b. Pilot equalization (depends on step 1; extends step 3)
```
./run_weighted_topup_pilot.sh
```
Drives `weighted_topup_pilot.R <grid_point> <seed>` for the same 4 additional
seeds (1, 2025, 314, 777) at the same under-resolved grid points as the
reproducibility dataset's own top-up (step 4) — 12 fits, appended to a
**new** file, `weighted_topup_pilot.csv` (does not touch
`multistart_results_corrected_pilot.csv`). This equalizes both N=30 datasets
at 47 fits each, so the profile-shape comparison is not confounded by
differing optimization depth. **Must run before step 6 (Figure 1)**, which
now reads this file too.

### 4. Reproducibility profile — THREE drivers, in strict sequence, shared CSV
```
./run_gate_test_N30_Kss1263_repro.sh        # grid point 4 only (Kss=1.263, the profile centre), 5 seeds -> 5 fits
./run_gate_test_N30_repro_remaining.sh      # the remaining grid points 1,2,3,5,6,7 x same 5 seeds -> 30 fits
./run_weighted_topup_corrected_repro.sh     # 4 extra seeds at low-resolution grid points -> +12 fits
```
5 + 30 + 12 = **47 fits**, matching the paper's reported total exactly. All
three append-and-skip to the **same** `gate_test_N30_Kss1263_repro.csv`.
**Note the argument-signature difference, and why it exists**: the first
script's worker (`gate_test_N30_Kss1263_repro.R`) has `kss_value <- 1.263`
hardcoded — it only ever fits the profile's centre grid point, which is why it
takes a single argument, `<seed>`. It was originally a preliminary "gate"
check at that one point (hence the name) before the full 7-point profile was
run, and its output was later folded into the reported 47-fit profile rather
than discarded. The other two workers (`gate_test_N30_repro.R`,
`weighted_topup_corrected_repro.R`) profile arbitrary grid points and so each
take two arguments, `<grid_point> <seed>`. Do not assume a uniform signature
across the three drivers. Final row count: 47 fits + header = 48 lines.

### 5. N=80 warm-start (depends on step 4)
```
./run_warmstart_N80_Kss1263.sh
```
Drives `warmstart_N80_Kss1263.R <source_seed>` for seeds 42, 999, 2024. Each
run reads its starting parameter vector directly from
`gate_test_N30_Kss1263_repro.csv` (the low-OFV fits at those seeds) — **this
step must run after step 4 completes**, not before or during. Output:
`warmstart_N80_Kss1263.csv`.

### 6. Figure 1 (depends on steps 3, 3b, and 4)
```
Rscript plot_fig1_corrected_two_draws.R
```
Reads `gate_test_N30_Kss1263_repro.csv` (step 4),
`multistart_results_corrected_pilot.csv` (step 3), and
`weighted_topup_pilot.csv` (step 3b) — the pilot curve is now computed from
the equalized 47-fit profile (steps 3+3b combined), matching the
reproducibility dataset's own 47-fit structure. Output:
`fig1_corrected_two_draws.png`.

### 7. Ground-truth recovery (depends only on step 1's reproducibility dataset)
```
./run_ground_truth_recovery_repro.sh
```
Drives `ground_truth_recovery_repro.R <seed>` for the same 5 seeds as the
main profile analysis, with Kss fixed at its true simulation value (1.56,
not a profile grid point) on the reproducibility dataset. Output:
`ground_truth_recovery_repro.csv`.

### 8. OFV-at-truth evaluation (depends only on step 1's datasets)
```
./run_ofv_at_truth_repro.sh          # 4 independent repeats -> ofv_at_truth_repro.csv
./run_ofv_at_truth_pilot.sh          # 3 independent repeats -> ofv_at_truth_pilot.csv
./run_ofv_at_truth_pilot_29subj.sh   # 3 independent repeats -> ofv_at_truth_pilot_29subj.csv
```
Each evaluates the objective function at the fully true parameter vector (all
theta, omega, and residual-error terms fixed at their exact simulation
values) via FOCEi's empirical-Bayes-only mode, with no estimation — the
mechanism needed because nlmixr2's SAEM crashes when zero population
parameters are free. Each driver repeats the evaluation as independent,
freshly started R processes to verify determinism. The full-pilot evaluation
(`run_ofv_at_truth_pilot.sh`) returns a corrupted value (17656.89) because one
subject's individual optimization diverges to physically implausible values
(eta estimates tens of standard deviations from zero) where the ODE cannot be
integrated, producing undefined individual predictions rather than a genuine
likelihood contribution; `run_ofv_at_truth_pilot_29subj.sh` excludes that
subject and reproduces the paper's reported 29-subject value (4006.04)
exactly. Both pilot scripts are identical except for this one exclusion —
same true parameter vector, same FOCEi method, same default tolerances.
Independent of steps 2–7.

**Dependency summary:**
```
generate_*_data.R (any order)
        |
        +--> run_multistart_corrected_full.sh -------------------> N=80 profile
        |
        +--> run_multistart_corrected_pilot.sh
        |    -> run_multistart_corrected_pilot_topup.sh
        |       -> run_weighted_topup_pilot.sh (step 3b) -----------> pilot profile (equalized) --+
        |                                                                                          |
        +--> run_gate_test_N30_Kss1263_repro.sh                                                    |
             -> run_gate_test_N30_repro_remaining.sh                                                |
             -> run_weighted_topup_corrected_repro.sh --------------> repro profile ----------------+--> Figure 1
                        |
                        +----> run_warmstart_N80_Kss1263.sh (after repro profile)

run_ground_truth_recovery_repro.sh -------------------------------> ground-truth recovery (independent)
run_ofv_at_truth_repro.sh / _pilot.sh / _pilot_29subj.sh ----------> OFV-at-truth (independent)
```

---

## Repository contents

**Model**
- `denosumab_QSS_TMDD_v2.R` — the mass-balance-corrected QSS-TMDD model
  (rxode2), used by every downstream script.

**Data generation**
- `generate_corrected_full_data.R` — simulates the N=80 dataset.
- `generate_corrected_pilot_data.R` — simulates the first independent N=30
  dataset (pilot).
- `generate_corrected_repro_data.R` — simulates the second independent N=30
  dataset (reproducibility).

**N=80 profile**
- `multistart_profile_Kss_corrected_full.R` — one-fit worker (grid point,
  seed).
- `run_multistart_corrected_full.sh` — driver, 7 grid points × 5 seeds.

**Pilot profile**
- `multistart_profile_Kss_corrected_pilot.R` — one-fit worker.
- `run_multistart_corrected_pilot.sh`, `run_multistart_corrected_pilot_topup.sh`
  — two drivers, run in sequence (see "How to reproduce").

**Pilot equalization**
- `weighted_topup_pilot.R` — one-fit worker, same design as the
  reproducibility dataset's own top-up.
- `run_weighted_topup_pilot.sh` — driver, 4 seeds × 4 grid points (12 fits
  total, see "How to reproduce").

**Reproducibility profile**
- `gate_test_N30_Kss1263_repro.R` — one-fit worker, profile centre only
  (Kss=1.263).
- `gate_test_N30_repro.R` — one-fit worker, arbitrary grid point.
- `weighted_topup_corrected_repro.R` — one-fit worker, arbitrary grid point.
- `run_gate_test_N30_Kss1263_repro.sh`, `run_gate_test_N30_repro_remaining.sh`,
  `run_weighted_topup_corrected_repro.sh` — three drivers, run in sequence
  (see "How to reproduce").

**N=80 warm-start**
- `warmstart_N80_Kss1263.R` — one-fit worker; reads its starting vector from
  the reproducibility profile's CSV.
- `run_warmstart_N80_Kss1263.sh` — driver, 3 source seeds.

**Ground-truth recovery**
- `ground_truth_recovery_repro.R` — one-fit worker, Kss fixed at its true
  value (1.56).
- `run_ground_truth_recovery_repro.sh` — driver, 5 seeds.

**OFV-at-truth evaluation**
- `ofv_at_truth_repro.R`, `ofv_at_truth_pilot.R`, `ofv_at_truth_pilot_29subj.R`
  — one-evaluation workers (FOCEi, empirical-Bayes only, no estimation,
  default tolerances). The pilot and 29-subject-pilot scripts are identical
  except for subject-9 exclusion.
- `run_ofv_at_truth_repro.sh`, `run_ofv_at_truth_pilot.sh`,
  `run_ofv_at_truth_pilot_29subj.sh` — drivers, 4/3/3 independent repeats
  respectively (determinism verification).

**Figure**
- `plot_fig1_corrected_two_draws.R` — regenerates Figure 1 from the pilot
  (equalized, steps 3+3b) and reproducibility (step 4) profile CSVs.

**Data** (shipped alongside the scripts that regenerate them — belt and
suspenders: use the provided data directly, or regenerate it yourself)
- `simulated_denosumab_data_v2_corrected_full.csv` — N=80 input dataset.
- `simulated_denosumab_data_v2_corrected_pilot.csv` — N=30 pilot input dataset.
- `simulated_denosumab_data_v2_corrected_repro.csv` — N=30 reproducibility
  input dataset.
- `multistart_results_corrected_full.csv` — N=80 profile results.
- `multistart_results_corrected_pilot.csv` — pilot profile results (35 fits;
  combine with `weighted_topup_pilot.csv` for the equalized 47-fit profile).
- `weighted_topup_pilot.csv` — pilot equalization top-up results (12 fits).
- `gate_test_N30_Kss1263_repro.csv` — reproducibility profile results
  (47 fits).
- `warmstart_N80_Kss1263.csv` — N=80 warm-start results.
- `ground_truth_recovery_repro.csv` — ground-truth recovery results (5 fits).
- `ofv_at_truth_repro.csv` — OFV-at-truth, reproducibility dataset (4 repeats,
  deterministic).
- `ofv_at_truth_pilot.csv` — OFV-at-truth, pilot dataset, all 30 subjects
  (3 repeats, deterministic; **corrupted** — see "How to reproduce", step 8).
- `ofv_at_truth_pilot_29subj.csv` — OFV-at-truth, pilot dataset, subject 9
  excluded (3 repeats, deterministic; reproduces the paper's reported value).
- `fig1_corrected_two_draws.png` — Figure 1.

**Model equations and parameters**
- `SUPPLEMENT.md` — the paper's Supplementary Material: full model
  equations and the ground-truth parameter table (sourced from Choi et al.
  2025, Table 3).

**Docs**
- `README.md` — this file.
- `LICENSE` — MIT.
- `CITATION.cff` — machine-readable citation metadata.

---

## What each stage produces (maps to the paper's reported results)

| Output file | Paper location |
|---|---|
| `multistart_results_corrected_full.csv` | Results, N=80 analysis (best OFV per grid point ~6739–8928; full set of fits extends to 11118) |
| `multistart_results_corrected_pilot.csv` + `weighted_topup_pilot.csv` (combined, 47 fits) | Results (pilot ΔOFV max 26.04, interior minimum at Kss=1.263, on equalized footing); Figure 1 |
| `gate_test_N30_Kss1263_repro.csv` | Results (repro ΔOFV: 0.00/13.38/5.82/4.73/1.39/6.65/4.36); Figure 1 |
| `warmstart_N80_Kss1263.csv` | Results (migration to OFV 7373.774/8031.354/8322.857) |
| `ground_truth_recovery_repro.csv` | Results (recovery-failure %-deviations: ka +79.5%, Vc +70.3%, Vp −26.2%, kint −33.2%, ksyn −47.1%, R0 +22.1%; CL −4.0%) |
| `ofv_at_truth_repro.csv` | Results (OFV at true parameter vector = 4746.271, falling in the high-OFV region, not the low) |
| `ofv_at_truth_pilot.csv` | Results (initial pilot evaluation = 17656.89; identified as corrupted — subject 9's individual optimization diverged) |
| `ofv_at_truth_pilot_29subj.csv` | Results (29-subject re-evaluation = 4006.04, reported as directionally consistent, not independent confirmation) |
| `fig1_corrected_two_draws.png` | Figure 1 |
| `SUPPLEMENT.md` | Supplementary Material (equations, parameter table) |

---

## Citation

If you use this code or data, please cite:

> Ravikumar, J.A. (2026+). Estimation Fragility in Multi-Start
> Profile-Likelihood Assessment of Kss in a QSS-TMDD Model of Denosumab.
> *Journal of Pharmacokinetics and Pharmacodynamics*. [DOI to be added on
> publication.]

Archived code/data DOI (Zenodo): *[to be added once minted]*

See `CITATION.cff` in this repository for structured citation metadata
(ORCID and DOI fields to be added once available).

---

## License

MIT — see `LICENSE`.
