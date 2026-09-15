# Supplementary Material

## S1. Model Equations

Verified against `denosumab_QSS_TMDD_v2.R` (commit `8fb40db`, the mass-balance-corrected
model; file confirmed unchanged since that commit) — every term below checked line-for-line
against the source rxode2 code before inclusion here.

### S1.1 Model structure and compartments

The model structure — a depot absorption compartment, central and peripheral drug
compartments, and a target (RANKL) compartment, linked through a quasi-steady-state
approximation for free drug — follows Choi et al. (2025). Four state variables are tracked
over time:

- `depot`: amount of drug in the absorption depot (nmol)
- `Ctot`: total (free + receptor-bound) drug concentration in the central compartment (nmol/L)
- `Cp`: drug concentration in the peripheral compartment (nmol/L)
- `Rtot`: total (free + drug-bound) receptor (RANKL) concentration (nmol/L)

Two further quantities are computed algebraically at each time point from the
quasi-steady-state approximation, rather than tracked as separate differential states:

- `Cfree`: free (unbound) drug concentration in the central compartment (nmol/L)
- `RC`: drug–receptor complex concentration in the central compartment (nmol/L)

### S1.2 Quasi-steady-state relations

$$\Delta = (C_{tot} - R_{tot} - K_{ss})^2 + 4 K_{ss} C_{tot}$$
$$\Delta^+ = \max(\Delta,\ 0)$$
$$C = \tfrac{1}{2}\left[(C_{tot} - R_{tot} - K_{ss}) + \sqrt{\Delta^+}\right]$$
$$C_{free} = \max(C,\ 0)$$
$$RC = C_{tot} - C_{free}$$

$\Delta$ and $C$ are each floored at zero before use, guarding against small negative values
arising from numerical round-off near the boundary of the feasible region.

### S1.3 Differential equations

$$\frac{d(depot)}{dt} = -k_a \cdot depot$$

$$\frac{dC_{tot}}{dt} = \frac{k_a \cdot depot}{V_c} - \frac{CL}{V_c} C_{free} - \frac{Q}{V_c} C_{free} + \frac{Q}{V_c} C_p - k_{int} \cdot RC$$

$$\frac{dC_p}{dt} = \frac{Q}{V_p} C_{free} - \frac{Q}{V_p} C_p$$

$$\frac{dR_{tot}}{dt} = k_{syn} - k_{deg}(R_{tot} - RC) - k_{int} \cdot RC, \qquad k_{deg} = \frac{k_{syn}}{R_0}$$

Inter-compartmental drug transfer is driven by a shared amount-flux — $Q \cdot C_{free}$ from
central to peripheral, and $Q \cdot C_p$ from peripheral to central — with each compartment's
concentration-rate equation dividing that flux by its own volume ($V_c$ for the
central-compartment terms, $V_p$ for the peripheral-compartment terms), so that drug amount
is conserved between the two compartments.

### S1.4 Units

| Quantity | Units |
|---|---|
| depot | nmol |
| Ctot, Cp, Cfree, RC | nmol/L |
| Rtot, R0, Kss | nmol/L |
| ka, kint, kdeg | 1/h |
| Vc, Vp | L |
| CL, Q | L/h |
| ksyn | nmol·L⁻¹·h⁻¹ |

## S2. Parameter Table

Ground-truth parameter values used to simulate all datasets in this study were drawn from
Choi et al. (2025) Table 3, using their postmenopausal osteoporosis (PMO) population
estimates consistently for parameters with a reported healthy-volunteer/PMO split (ka, R0,
Q), together with their Caucasian apparent clearance (CL/F) and apparent volume estimates
(VC/F, VP/F). Values are transcribed exactly as used in `generate_corrected_full_data.R`
(commit `8fb40db`), cross-checked against Table 3 directly.

### Table S1. Population parameters

| Parameter | Value | Unit | Source (Choi et al. 2025, Table 3) |
|---|---|---|---|
| ka | 0.0078 | 1/h | ka_PMO |
| Vc | 1.58 | L | VC/F |
| Vp | 6.06 | L | VP/F |
| CL | 0.006 | L/h | CL/F, Caucasian |
| Q | 0.20 | L/h | Q/F_PMO |
| kint | 0.022 | 1/h | Table 3 (no subgroup split) |
| Kss | 1.56 | nmol/L | Table 3 (no subgroup split) |
| ksyn | 0.01 | nmol·L⁻¹·h⁻¹ | Table 3 (no subgroup split) |
| R0 | 15.23 | nmol/L | R0_PMO |
| kdeg | 6.566×10⁻⁴ | 1/h | Derived: kdeg = ksyn / R0 (not a direct Table 3 value) |

**Residual error model** (Choi et al. 2025, Table 3): additive σ = 0.72 nmol/L,
proportional σ = 0.07 (7%).

### Table S2. Inter-individual variability (IIV)

Inter-individual variability was simulated using Choi et al. (2025) Table 3's reported
%CV for each parameter, converted to log-normal variance via ω² = ln(1 + (CV/100)²).
The %CV values are Choi's reported quantities; the ω² column is our derived transformation
of them, not a value reported directly in Table 3.

| Parameter | IIV %CV (Choi Table 3) | ω² = ln(1+(CV/100)²) |
|---|---|---|
| ka | 56.57% | 0.2776 |
| Vc | 61.69% | 0.3225 |
| Vp | 15.55% | 0.0239 |
| CL | 26.39% | 0.0673 |
| Q | 295.99% | 2.2784 |
| kint | 7.88% | 0.0062 |
| Kss | 58.12% | 0.2910 |
| ksyn | 22.46% | 0.0492 |
| R0 | 158.3% | 1.2544 |
