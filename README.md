# README draft for the curated public reproducibility repo
(Working draft — assembled here section by section, verified against source
before each is locked. Not the repo's actual root README yet; that gets
created when the clean folder is assembled per the manifest in FLAG_STATUS.md.)

---

## Estimation Fragility in Multi-Start Profile-Likelihood Assessment of Kss in a QSS-TMDD Model of Denosumab

Code and data to reproduce the analysis in "Estimation Fragility in
Multi-Start Profile-Likelihood Assessment of Kss in a QSS-TMDD Model of
Denosumab." Simulation, multi-start SAEM profiling, and warm-start
diagnostics for a mass-balance-corrected QSS-TMDD model of denosumab, in R
(nlmixr2 / rxode2).

**Status:** title/one-liner matches the paper's own title exactly (see
`paper_draft.md`); no separate claim made here.

---

## What this reproduces

This analysis examines the practical identifiability of the quasi-steady-state
constant (Kss) in a QSS-TMDD model of denosumab, using profile likelihood with
multi-start SAEM estimation on a mass-balance-corrected model. The principal
finding is that single-fit SAEM profile-likelihood estimation is fragile in
this setting: fits reach different parameter basins depending on the random
seed, and the resulting profile shape is not reproducible across independent
simulated N=30 datasets. This behavior did not generalize to a larger (N=80)
dataset — profile fits there converged to a substantially higher
objective-function-value regime, and warm-starting from the N=30 low-OFV basin
failed to recover it, indicating the identifiability behavior observed at N=30
is sample-size-dependent. Kss is only weakly constrained by pharmacokinetic-only
data, consistent with the shrinkage-based limitation reported by Choi et al.
Robust profiling therefore requires multi-start estimation with explicit
assessment of basin and convergence consistency.

**Status:** LOCKED (checked against all three of the paper's reported
result-legs and both framing bans — no "confirms," no "one-sided").

---

## Requirements

- R 4.6.1
- rxode2 5.1.5
- nlmixr2 7.0.1
- dplyr 1.2.1
- ggplot2 4.0.3 (used only by the corrected model script and the figure script)

**Status:** verified directly against the installed environment and each
script's `library()`/`suppressPackageStartupMessages()` calls, not assumed.

---

## How to reproduce (run order)

The pipeline has two independently-runnable stages (data generation) followed
by four analysis stages, two of which are themselves multi-driver sequences
that must run in order because they append to a shared results file. Running
steps out of order, or skipping a driver in a sequence, produces an incomplete
or wrong result — this is not a stylistic preference, it is how the crash-safe
append-and-skip pattern used throughout this pipeline works.

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

### 4. Reproducibility profile — THREE drivers, in strict sequence, shared CSV
```
./run_gate_test_N30_Kss1263_repro.sh        # grid point 4 only (Kss=1.263, the profile centre), 5 seeds -> 5 fits
./run_gate_test_N30_repro_remaining.sh      # the remaining grid points 1,2,3,5,6,7 x same 5 seeds -> 30 fits
./run_weighted_topup_corrected_repro.sh     # 4 extra seeds at low-resolution grid points -> +12 fits
```
5 + 30 + 12 = **47 fits**, matching Methods 2.4's reported total exactly. All
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

### 6. Figure 1 (depends on steps 3 and 4)
```
Rscript plot_fig1_corrected_two_draws.R
```
Reads `gate_test_N30_Kss1263_repro.csv` (step 4) and
`multistart_results_corrected_pilot.csv` (step 3); computes each dataset's
best-of-N profile with ΔOFV relative to its own minimum, and overlays them.
Output: `fig1_corrected_two_draws.png`.

**Dependency summary:**
```
generate_*_data.R (any order)
        |
        +--> run_multistart_corrected_full.sh -------------------> N=80 profile
        |
        +--> run_multistart_corrected_pilot.sh
        |    -> run_multistart_corrected_pilot_topup.sh -----------> pilot profile --+
        |                                                                            |
        +--> run_gate_test_N30_Kss1263_repro.sh                                      |
             -> run_gate_test_N30_repro_remaining.sh                                 |
             -> run_weighted_topup_corrected_repro.sh --------------> repro profile -+--> Figure 1
                        |                                                            
                        +----> run_warmstart_N80_Kss1263.sh (after repro profile)
```

**Status:** every invocation, argument signature, and file dependency in this
section was traced directly from the actual scripts (commandArgs() calls,
results_csv paths, driver .sh loops) and cross-checked (row counts, unique
seed sets) against the numbers reported in the paper — not written from
memory or from script names.

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

**Figure**
- `plot_fig1_corrected_two_draws.R` — regenerates Figure 1 from the pilot
  and reproducibility profile CSVs.

**Data** (shipped alongside the scripts that regenerate them — belt and
suspenders: use the provided data directly, or regenerate it yourself)
- `simulated_denosumab_data_v2_corrected_full.csv` — N=80 input dataset.
- `simulated_denosumab_data_v2_corrected_pilot.csv` — N=30 pilot input dataset.
- `simulated_denosumab_data_v2_corrected_repro.csv` — N=30 reproducibility
  input dataset.
- `multistart_results_corrected_full.csv` — N=80 profile results.
- `multistart_results_corrected_pilot.csv` — pilot profile results.
- `gate_test_N30_Kss1263_repro.csv` — reproducibility profile results
  (47 fits).
- `warmstart_N80_Kss1263.csv` — N=80 warm-start results.
- `fig1_corrected_two_draws.png` — Figure 1.

**Model equations and parameters**
- `SUPPLEMENT.md` — the paper's Supplementary Material: full model
  equations (verified against `denosumab_QSS_TMDD_v2.R`) and the
  ground-truth parameter table (verified against Choi et al. 2025, Table 3).
  *(To create: this is the verified content of `supplement_draft.md` S1+S2,
  copied in — not the working file itself, which stays private.)*

**Docs**
- `README.md` — this file.
- `LICENSE` — MIT.
- `CITATION.cff` — machine-readable citation metadata. *(To create.)*

**Status:** every IN file above is exactly the locked manifest recorded in
this repo's `FLAG_STATUS.md` — nothing added or dropped here.

---

## What each stage produces (maps to the paper's reported results)

| Output file | Paper location |
|---|---|
| `multistart_results_corrected_full.csv` | Results ¶3, N=80 analysis (best OFV per grid point ~6739–8928; full set of fits extends to 11118) |
| `multistart_results_corrected_pilot.csv` | Results ¶1 (pilot ΔOFV max 26.04, interior minimum at Kss=1.263); Figure 1 |
| `gate_test_N30_Kss1263_repro.csv` | Results ¶1, ¶3 (repro ΔOFV: 0.00/13.38/5.82/4.73/1.39/6.65/4.36); Figure 1 |
| `warmstart_N80_Kss1263.csv` | Results ¶3 (migration to OFV 7373.774/8031.354/8322.857) |
| `fig1_corrected_two_draws.png` | Figure 1 |
| `SUPPLEMENT.md` | Supplementary Material (equations, parameter table) |

**Status:** derived directly from the run-order section above and the
values already verified against `paper_draft.md` in earlier sessions
(recorded in `FLAG_STATUS.md`) — not new claims, a consolidation of what's
already been checked.

---

## Citation

If you use this code or data, please cite:

> Ravikumar, J.A. (2026+). Estimation Fragility in Multi-Start
> Profile-Likelihood Assessment of Kss in a QSS-TMDD Model of Denosumab.
> *Journal of Pharmacokinetics and Pharmacodynamics*. [DOI to be added on
> publication.]

Archived code/data DOI (Zenodo): *[to be added once minted]*

A `CITATION.cff` file will be added with the same information once both
DOIs exist.

**Status:** author/title exactly match the confirmed manuscript title and
author block; journal name matches the confirmed target; both DOI fields
are honestly marked as not-yet-existing rather than invented.

---

## License

MIT — see `LICENSE`.

---

## Still to do (not README content — repo-assembly steps)
- Create `SUPPLEMENT.md` from the verified `supplement_draft.md` S1+S2
  content.
- Create `CITATION.cff`.
- Rewrite/replace the currently-stale `README.md` (see FLAG_STATUS.md —
  it says the analysis "confirms" the finding, which contradicts the
  paper's locked "consistent with" framing) with this file's content once
  finalized.
- Assemble the clean folder against the locked manifest (Mac, you) → push
  new public repo (GitHub, you) → Zenodo release → DOI → put the real URL
  into `paper_draft.md`'s Declarations data-availability line (here,
  together).
