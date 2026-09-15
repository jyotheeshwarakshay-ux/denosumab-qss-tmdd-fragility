# =============================================================================
# PROJECT 3 (v2): QSS-Approximation TMDD Model — Denosumab
# CORRECTED after systematic debugging of v1
#
# Ground truth: Choi S, Park S, Jung J, Baek S, Lim H-S (2025). Front Pharmacol
#               16:1631034. doi:10.3389/fphar.2025.1631034 (CC-BY). Table 3.
# Model theory: Dua P, Hawkins E, van der Graaf PH (2015). CPT PSP 4(6):324-337.
#               Gibiansky L et al. (2008). J Pharmacokinet Pharmacodyn 35:573-591.
# Author:       Jyotheeshwar Akshay Ravikumar
#
# -----------------------------------------------------------------------------
# WHAT WENT WRONG IN v1, AND HOW IT WAS FOUND
# -----------------------------------------------------------------------------
# v1 fit converged but underpredicted systematically: mean(DV-IPRED) = +5.01,
# >75% of residuals positive, max IPRED 70.6 against data max 120.8, and
# add.err inflated to 2.68 against a true 0.72.
#
# Six hypotheses were tested and eliminated by direct measurement:
#   1. Negative-value clipping bias   -> too small (-0.02, needed +5)
#   2. Sparse sampling design          -> data DID distinguish parameter sets
#                                          (RMSE 6.15 nmol/L between true and
#                                           fitted curves at sampled times)
#   3. Compartment (CMT) mis-mapping   -> IPRED was exactly equal to Ctot
#   4. Variance (eta) collapse         -> REAL, but fixing it left bias at +5.67
#   5. Simulation/fit model mismatch   -> models identical to ~1e-12
#   6. Data-handling in pk_data        -> pk_data and eventTable pathways agreed
#                                          (max Ctot 37.23 vs 37.24)
#
# ACTUAL CAUSE: sampling in deep washout.
#   The 3600 h and 4380 h post-dose offsets sit at ~1.4 and ~0.10 nmol/L on the
#   typical curve. With an additive residual error of 0.72 nmol/L (the value
#   Table 3 reports), the noise at the 4380 h trough is ~7x the signal. SAEM
#   weights observations by inverse variance, so hundreds of these
#   uninformative points dominated the likelihood and dragged the fitted curve
#   downward. The inflated add.err was SAEM absorbing that misspecification.
#
#   The source paper avoided this without needing to: their assay LLOQ was
#   20 ng/mL and 26.49% of post-dose samples fell below it and were set to
#   MDV=1 (excluded from the likelihood). I simulated the assay but not its
#   detection limit, generating observations their dataset never contained.
#
# THE TWO CORRECTIONS IN v2
#   (a) SAMPLING: drop the 3600 h and 4380 h washout offsets; add 72 h, 336 h,
#       and 1440 h in the informative absorption/distribution region. Minimum
#       typical concentration rises from 0.095 to 8.55 nmol/L; noise-to-signal
#       at the worst point falls from 7.0 to 0.08.
#   (b) INITIAL VARIANCES: initialise every eta AT its true omega^2 rather than
#       below it. In v1, eta_Vc was initialised at 0.20 against a true 0.3225
#       and collapsed 20-fold to 0.0156; eta_R0 was initialised at 0.50 against
#       a true 1.0195 and collapsed 40-fold. Confirmed fixed by diagnostic test
#       (eta_Vc ratio 1.15x).
#   (c) Q FIXED: in v1, Q had 499% RSE and BSV 902% CV -- it was not being
#       estimated in any meaningful sense. It is fixed to its true value here,
#       the same standard remedy applied to Ka in Project 1 v2.
#
# -----------------------------------------------------------------------------
# METHODOLOGY NOTE (must appear in the README verbatim)
# -----------------------------------------------------------------------------
#   Patient-level denosumab concentration data are not publicly available. What
#   is public is the FITTED PARAMETER TABLE (Choi et al. Table 3), estimated by
#   those authors from 6,583 real serum concentrations in 615 real subjects.
#   This script simulates a dataset FROM those published parameters, then
#   re-fits to recover them. This is simulation-based parameter recovery, NOT
#   an analysis of real patient data. State this plainly in the README.
#
# UNITS: concentrations nmol/L, volumes L, time hours. Dose converted from mg
# to nmol using denosumab MW ~147 kDa.
# =============================================================================


# =============================================================================
# STAGE 1: PACKAGES
# =============================================================================
library(nlmixr2)
library(rxode2)
library(ggplot2)
library(dplyr)
library(xpose.nlmixr2)

set.seed(42)


# =============================================================================
# STAGE 2: GROUND-TRUTH PARAMETERS (Choi et al. 2025, Table 3 — PMO population)
# =============================================================================
true_params <- list(
  ka    = 0.0078,   # 1/h   absorption rate constant (PMO patients)
  Vc    = 1.58,     # L     central volume (apparent, Vc/F)
  Vp    = 6.06,     # L     peripheral volume (apparent, Vp/F)
  CL    = 0.006,    # L/h   clearance (Caucasian, apparent CL/F)
  Q     = 0.20,     # L/h   inter-compartmental clearance (PMO)
  kint  = 0.022,    # 1/h   drug-target complex internalization rate
  Kss   = 1.56,     # nmol/L quasi-steady-state constant
  ksyn  = 0.01,     # 1/h   target (RANKL) synthesis rate
  R0    = 15.23     # nmol/L baseline target concentration (PMO)
)

# kdeg is NOT reported separately in Table 3. Derived from the steady-state
# assumption stated in the paper's Methods (their Eq. 7): ksyn = kdeg * R0.
# This derivation is documented here and MUST be documented in the README.
true_params$kdeg <- true_params$ksyn / true_params$R0
cat("Derived kdeg =", true_params$kdeg, "1/h  (= ksyn / R0)\n")

# IIV (%CV from Table 3) -> log-normal variance: omega^2 = log(1 + CV^2)
cv_to_omega2 <- function(cv_pct) log(1 + (cv_pct / 100)^2)

true_omega2 <- list(
  ka   = cv_to_omega2(56.57),
  Vc   = cv_to_omega2(61.69),
  Vp   = cv_to_omega2(15.55),
  kint = cv_to_omega2(7.88),
  Kss  = cv_to_omega2(58.12),
  ksyn = cv_to_omega2(22.46),
  R0   = cv_to_omega2(158.3),
  CL   = cv_to_omega2(26.39),
  Q    = cv_to_omega2(295.99)
)

true_sigma_add  <- 0.72   # nmol/L, Table 3
true_sigma_prop <- 0.07   # 7%,     Table 3

MW_kDa    <- 147
dose_mg   <- 60
dose_nmol <- (dose_mg * 1e6) / (MW_kDa * 1e3)
cat("Dose in nmol:", round(dose_nmol, 1), "\n")


# =============================================================================
# STAGE 3: SIMULATE THE DATASET  ** SAMPLING DESIGN CORRECTED **
# =============================================================================
n_subjects <- 80
month_h    <- 30.44 * 24

# ---- CORRECTION (a): informative sampling only -------------------------------
# v1 used: c(12, 24, 168, 720, 2160, 3600, 4380)
#   -> 3600 h gives ~1.41 nmol/L and 4380 h gives ~0.095 nmol/L on the typical
#      curve. Against add.err = 0.72 nmol/L that is noise-dominated, and those
#      points drove the systematic underprediction diagnosed above.
# v2 uses the schedule below: washout offsets removed, absorption/distribution
#   coverage improved. Minimum typical concentration 8.55 nmol/L.
sample_offsets_h <- c(12, 24, 72, 168, 336, 720, 1440, 2160)
dose_times_h     <- c(0, 6, 12) * month_h

cat("\nSampling offsets (h post-dose):", paste(sample_offsets_h, collapse = ", "), "\n")
cat("  (v1 used 12,24,168,720,2160,3600,4380 -- last two were washout)\n\n")

qss_model_rxode <- rxode2({
  ka    <- exp(lka + eta_ka)
  Vc    <- exp(lVc + eta_Vc)
  Vp    <- exp(lVp + eta_Vp)
  CL    <- exp(lCL + eta_CL)
  Q     <- exp(lQ  + eta_Q)
  kint  <- exp(lkint + eta_kint)
  Kss   <- exp(lKss + eta_Kss)
  ksyn  <- exp(lksyn + eta_ksyn)
  R0    <- exp(lR0 + eta_R0)
  kdeg  <- ksyn / R0

  disc  <- (Ctot - Rtot - Kss)^2 + 4 * Kss * Ctot
  discP <- max(disc, 0)
  C     <- 0.5 * ((Ctot - Rtot - Kss) + sqrt(discP))
  Cfree <- max(C, 0)
  RC    <- Ctot - Cfree

  d/dt(depot) <- -ka * depot
  d/dt(Ctot)  <- (ka * depot) / Vc - (CL / Vc) * Cfree - (Q / Vc) * Cfree +
                 (Q / Vc) * Cp - kint * RC
  d/dt(Cp)    <- (Q / Vp) * Cfree - (Q / Vp) * Cp
  d/dt(Rtot)  <- ksyn - kdeg * (Rtot - RC) - kint * RC
})

sim_params <- data.frame(
  id       = 1:n_subjects,
  eta_ka   = rnorm(n_subjects, 0, sqrt(true_omega2$ka)),
  eta_Vc   = rnorm(n_subjects, 0, sqrt(true_omega2$Vc)),
  eta_Vp   = rnorm(n_subjects, 0, sqrt(true_omega2$Vp)),
  eta_CL   = rnorm(n_subjects, 0, sqrt(true_omega2$CL)),
  eta_Q    = rnorm(n_subjects, 0, sqrt(true_omega2$Q)),
  eta_kint = rnorm(n_subjects, 0, sqrt(true_omega2$kint)),
  eta_Kss  = rnorm(n_subjects, 0, sqrt(true_omega2$Kss)),
  eta_ksyn = rnorm(n_subjects, 0, sqrt(true_omega2$ksyn)),
  eta_R0   = rnorm(n_subjects, 0, sqrt(true_omega2$R0))
)

theta <- log(c(
  true_params$ka, true_params$Vc, true_params$Vp, true_params$CL,
  true_params$Q, true_params$kint, true_params$Kss,
  true_params$ksyn, true_params$R0
))
names(theta) <- c("lka","lVc","lVp","lCL","lQ","lkint","lKss","lksyn","lR0")

sim_results <- list()

for (i in 1:n_subjects) {
  p    <- sim_params[i, ]
  R0_i <- exp(theta["lR0"] + p$eta_R0)

  # NOTE: eventTable() called WITHOUT units arguments. Specifying
  # amount.units/time.units requires the 'units' package, which is not
  # installed on Posit Cloud; this caused a silent loop failure in v1.
  ev <- eventTable()
  for (dt in dose_times_h) {
    ev$add.dosing(dose = dose_nmol, start.time = dt, dosing.to = "depot")
  }
  samp_times <- sort(unique(as.vector(outer(dose_times_h, sample_offsets_h, "+"))))
  samp_times <- samp_times[samp_times <= max(dose_times_h) + max(sample_offsets_h)]
  ev$add.sampling(samp_times)

  inits <- c(depot = 0, Ctot = 0, Cp = 0, Rtot = as.numeric(R0_i))

  subj_params <- c(theta,
                   eta_ka = p$eta_ka, eta_Vc = p$eta_Vc, eta_Vp = p$eta_Vp,
                   eta_CL = p$eta_CL, eta_Q = p$eta_Q, eta_kint = p$eta_kint,
                   eta_Kss = p$eta_Kss, eta_ksyn = p$eta_ksyn, eta_R0 = p$eta_R0)

  sol           <- rxSolve(qss_model_rxode, subj_params, ev, inits = inits)
  sol_df        <- as.data.frame(sol)
  sol_df$ID     <- i
  sim_results[[i]] <- sol_df
}

sim_all <- bind_rows(sim_results)

# Residual error. No pmax() clipping this time: with washout sampling removed,
# concentrations never approach zero, so negative simulated values do not arise
# in practice. Any that did would signal a problem worth seeing, not hiding.
sim_all <- sim_all %>%
  mutate(
    IPRED    = Ctot,
    err_add  = rnorm(n(), 0, true_sigma_add),
    err_prop = rnorm(n(), 0, true_sigma_prop),
    DV       = IPRED * (1 + err_prop) + err_add
  )

cat("Simulated dataset:", nrow(sim_all), "records across", n_subjects, "subjects\n")
cat("DV summary:\n"); print(summary(sim_all$DV))
cat(sprintf("\n  negative DV values: %d  (should be 0)\n", sum(sim_all$DV < 0)))
cat(sprintf("  DV below 1 nmol/L : %d  (v1 had hundreds)\n\n", sum(sim_all$DV < 1)))


# =============================================================================
# STAGE 4: FORMAT AS NLMIXR2 (NONMEM-STYLE) DATASET
# =============================================================================
obs_rows <- sim_all %>%
  filter(time %in% unlist(lapply(dose_times_h, function(dt) dt + sample_offsets_h))) %>%
  transmute(ID, TIME = time, AMT = 0, DV = DV, EVID = 0, CMT = 2)

dose_rows <- expand.grid(ID = 1:n_subjects, TIME = dose_times_h) %>%
  transmute(ID, TIME, AMT = dose_nmol, DV = 0, EVID = 1, CMT = 1)

pk_data <- bind_rows(obs_rows, dose_rows) %>%
  arrange(ID, TIME, desc(EVID))

dir.create("outputs", showWarnings = FALSE)
write.csv(pk_data, "outputs/simulated_denosumab_data_v2.csv", row.names = FALSE)
cat("Dataset saved. Rows:", nrow(pk_data),
    " (obs:", sum(pk_data$EVID == 0), " dose:", sum(pk_data$EVID == 1), ")\n\n")


# =============================================================================
# STAGE 5: FIT  ** INITIAL VARIANCES CORRECTED, Q FIXED **
# =============================================================================
qss_fit_model_v2 <- function() {
  ini({
    lka   <- log(0.01)
    lVc   <- log(2)
    lVp   <- log(5)
    lCL   <- log(0.008)
    lQ    <- fixed(log(0.20))   # CORRECTION (c): 499% RSE in v1, not estimable
    lkint <- log(0.02)
    lKss  <- log(1.5)
    lksyn <- log(0.01)
    lR0   <- log(12)

    # CORRECTION (b): initialised AT true omega^2, not below.
    eta_ka   ~ 0.3225
    eta_Vc   ~ 0.3225
    eta_Vp   ~ 0.0239
    eta_CL   ~ 0.0680
    eta_kint ~ 0.0062
    eta_Kss  ~ 0.2870
    eta_ksyn ~ 0.0492
    eta_R0   ~ 1.0195

    add.err  <- 0.72
    prop.err <- 0.07
  })

  model({
    ka    <- exp(lka + eta_ka)
    Vc    <- exp(lVc + eta_Vc)
    Vp    <- exp(lVp + eta_Vp)
    CL    <- exp(lCL + eta_CL)
    Q     <- exp(lQ)
    kint  <- exp(lkint + eta_kint)
    Kss   <- exp(lKss + eta_Kss)
    ksyn  <- exp(lksyn + eta_ksyn)
    R0    <- exp(lR0 + eta_R0)
    kdeg  <- ksyn / R0

    disc  <- (Ctot - Rtot - Kss)^2 + 4 * Kss * Ctot
    discP <- max(disc, 0)
    C     <- 0.5 * ((Ctot - Rtot - Kss) + sqrt(discP))
    Cfree <- max(C, 0)
    RC    <- Ctot - Cfree

    d/dt(depot) <- -ka * depot
    d/dt(Ctot)  <- (ka * depot) / Vc - (CL / Vc) * Cfree - (Q / Vc) * Cfree +
                   (Q / Vc) * Cp - kint * RC
    d/dt(Cp)    <- (Q / Vp) * Cfree - (Q / Vp) * Cp
    d/dt(Rtot)  <- ksyn - kdeg * (Rtot - RC) - kint * RC

    Rtot(0) <- R0

    Ctot ~ add(add.err) + prop(prop.err)
  })
}

cat("Starting SAEM fit (expect roughly 40-60 minutes)...\n\n")

fit2 <- nlmixr2(
  object  = qss_fit_model_v2,
  data    = pk_data,
  est     = "saem",
  control = saemControl(seed = 42, print = 50, nBurn = 300, nEm = 400)
)

saveRDS(fit2, "outputs/fit_v2.rds")
print(fit2)


# =============================================================================
# STAGE 6: COMPARE RECOVERED PARAMETERS AGAINST GROUND TRUTH
# =============================================================================
fe <- fixef(fit2)

comparison <- data.frame(
  parameter = c("ka","Vc","Vp","CL","kint","Kss","ksyn","R0"),
  true_value = c(true_params$ka, true_params$Vc, true_params$Vp,
                 true_params$CL, true_params$kint, true_params$Kss,
                 true_params$ksyn, true_params$R0),
  recovered = exp(c(fe[["lka"]], fe[["lVc"]], fe[["lVp"]], fe[["lCL"]],
                    fe[["lkint"]], fe[["lKss"]], fe[["lksyn"]], fe[["lR0"]]))
)
comparison$pct_error <- 100 * (comparison$recovered - comparison$true_value) /
                        comparison$true_value

cat("\n======================================================================\n")
cat("PARAMETER RECOVERY  (Q was fixed, so not shown)\n")
cat("======================================================================\n")
print(comparison, row.names = FALSE, digits = 4)
write.csv(comparison, "outputs/parameter_recovery_v2.csv", row.names = FALSE)

# --- the diagnostic that failed in v1 ---------------------------------------
ires <- fit2$DV - fit2$IPRED
cat("\nSYSTEMATIC BIAS CHECK (v1 was badly broken here):\n")
cat(sprintf("  mean(DV - IPRED) : %+8.4f   (v1: +5.014  | want near 0)\n",
            mean(ires, na.rm = TRUE)))
cat(sprintf("  %% positive residuals: %6.1f%%   (v1: >75%%    | want ~50%%)\n",
            100 * mean(ires > 0, na.rm = TRUE)))
cat(sprintf("  max IPRED        : %8.1f   (v1: 70.6    | data max %.1f)\n",
            max(fit2$IPRED, na.rm = TRUE),
            max(pk_data$DV[pk_data$EVID == 0], na.rm = TRUE)))
cat(sprintf("  add.err estimate : %8.4f   (v1: 2.683   | true 0.72)\n",
            fe[["add.err"]]))
cat(sprintf("  prop.err estimate: %8.4f   (v1: 0.589   | true 0.07)\n",
            fe[["prop.err"]]))

cat("\nSHRINKAGE (expected to stay HIGH for TMDD params -- this is the finding,\n")
cat("not a defect; Choi et al. report the same limitation):\n")
print(fit2$parFixedDf[, c("Back-transformed", "%RSE", "BSV(CV%)", "Shrink(SD)%")])


# =============================================================================
# STAGE 7: DIAGNOSTICS
# =============================================================================
xpdb <- xpose_data_nlmixr2(fit2)

p1 <- dv_vs_pred(xpdb)  + labs(title = "DV vs PRED (v2)")
p2 <- dv_vs_ipred(xpdb) + labs(title = "DV vs IPRED (v2)")
p3 <- res_vs_idv(xpdb, res = "CWRES") +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
      labs(title = "CWRES vs Time (v2)")
p4 <- res_vs_pred(xpdb, res = "CWRES") +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
      labs(title = "CWRES vs PRED (v2)")

print(p1); print(p2); print(p3); print(p4)

ggsave("outputs/v2_01_dv_vs_pred.png",  p1, width = 6, height = 5)
ggsave("outputs/v2_02_dv_vs_ipred.png", p2, width = 6, height = 5)
ggsave("outputs/v2_03_cwres_time.png",  p3, width = 6, height = 5)
ggsave("outputs/v2_04_cwres_pred.png",  p4, width = 6, height = 5)


# =============================================================================
# STAGE 8: VISUAL PREDICTIVE CHECK
# =============================================================================
vpc_v2 <- vpcPlot(fit2, n = 300, bins = 8,
                  pi = c(0.05, 0.95), ci = c(0.025, 0.975))
print(vpc_v2)
ggsave("outputs/v2_05_vpc.png", vpc_v2, width = 7, height = 5)

writeLines(capture.output(sessionInfo()), "outputs/session_info_v2.txt")

cat("\n=== v2 complete. Outputs in /outputs/ ===\n")
cat("NEXT: write the README. It must cover (1) the ground-truth source and\n")
cat("its CC-BY licence, (2) the simulation-based methodology stated plainly,\n")
cat("(3) the v1 -> v2 debugging trail, (4) parameter recovery, and\n")
cat("(5) the shrinkage finding versus Choi et al.'s reported limitation.\n")
