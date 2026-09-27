#!/usr/bin/env Rscript
# =============================================================================
# GROUND-TRUTH RECOVERY (Experiment C) -- ONE fit per invocation, Kss FIXED
# at its TRUE simulation value (1.56, NOT the profile grid centre 1.263), on
# the reproducibility INDEPENDENT 30-subject corrected dataset, with the SAME
# SAEM settings as the main profile analysis (nBurn=150, nEm=200).
#
# Usage:  Rscript ground_truth_recovery_repro.R <seed>
#
# Question this answers: conditional on Kss being known exactly, does the
# rest of the model (ka, Vc, Vp, CL, kint, ksyn, R0, add.err, prop.err)
# recover near its true simulation value? Multi-start (5 seeds, the same
# standard set used throughout: 42, 123, 999, 7, 2024) because a single
# recovery fit would itself be subject to the estimation fragility this
# paper documents -- report the lowest-OFV fit's recovery AND the spread
# across seeds, not a single clean number.
#
# True generating values, pulled directly from generate_corrected_pilot_data.R
# and generate_corrected_repro_data.R (both use the SAME generative model --
# only the simulated data realization differs):
#   ka=0.0078  Vc=1.58  Vp=6.06  CL=0.006  Q=0.20 (fixed everywhere, not
#   estimated -- recovery does not apply to it)  kint=0.022  Kss=1.56
#   ksyn=0.01  R0=15.23  add.err(true sim noise)=0.72  prop.err(true sim
#   noise)=0.07
# Kss=1.56 does not coincide with any existing profile grid point (grid
# points 4 and 5 are 1.263 and 2.082) -- this is a genuinely new fitted
# point, not a rerun of an already-reported grid fit.
#
# Does NOT touch gate_test_N30_Kss1263_repro.csv or any other existing
# results file. Writes ONLY to ground_truth_recovery_repro.csv (new file)
# and ground_truth_recovery_repro_logs/ (new directory).
#
# Model equations, init values, eta structure, and Q=0.20 are IDENTICAL to
# every other corrected-model script in this repo (verbatim copy of
# weighted_topup_corrected_repro.R's build_model_corrected()), changing
# only the fixed Kss value (1.56 instead of a profile grid point).
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript ground_truth_recovery_repro.R <seed>")
}
seed <- as.integer(args[1])
if (is.na(seed)) stop("seed must be an integer")

repo        <- path.expand("~/pharmacometrics/denosumab-tmdd-qss")
setwd(repo)
data_csv    <- file.path(repo, "simulated_denosumab_data_v2_corrected_repro.csv")
results_csv <- file.path(repo, "ground_truth_recovery_repro.csv")

# Kss fixed at the TRUE simulation value, not a profile grid point.
true_kss     <- 1.56
kss_log      <- log(true_kss)

# ---------------------------------------------------------------------------
# RESUME CHECK -- skip if this seed already has a row.
# ---------------------------------------------------------------------------
if (file.exists(results_csv)) {
  existing <- tryCatch(read.csv(results_csv, stringsAsFactors = FALSE),
                        error = function(e) NULL)
  if (!is.null(existing) && nrow(existing) > 0 &&
      "seed" %in% names(existing) && any(existing$seed == seed)) {
    cat(sprintf("SKIP seed=%d -- already recorded in %s\n",
                seed, basename(results_csv)))
    quit(save = "no", status = 0)
  }
}

if (!file.exists(data_csv)) {
  stop("Repro dataset not found: ", data_csv,
       " -- run generate_corrected_repro_data.R first.")
}

suppressPackageStartupMessages({
  library(nlmixr2)
  library(rxode2)
  library(dplyr)
})

pk_prof <- read.csv(data_csv)
if (length(unique(pk_prof$ID)) != 30) {
  stop("Expected 30 subjects in ", data_csv, ", found ", length(unique(pk_prof$ID)))
}

# ---------------------------------------------------------------------------
# MODEL BUILDER -- verbatim copy of weighted_topup_corrected_repro.R's
# build_model_corrected(), Kss fixed at the supplied value, CORRECTED
# transfer terms. Same init values, same eta structure, same Q=0.20 as every
# other corrected-model script in this repo.
# ---------------------------------------------------------------------------
build_model_corrected <- function(fixed_value) {

  init <- list(
    lka   = log(0.00876),
    lVc   = log(2.169),
    lVp   = log(6.190),
    lCL   = log(0.00706),
    lkint = log(0.01680),
    lKss  = fixed_value,
    lksyn = log(0.00765),
    lR0   = log(15.638)
  )

  ini_lines <- character(0)
  for (nm in names(init)) {
    if (nm == "lKss") {
      ini_lines <- c(ini_lines, sprintf("    %s <- fixed(%.10f)", nm, init[[nm]]))
    } else {
      ini_lines <- c(ini_lines, sprintf("    %s <- %.10f", nm, init[[nm]]))
    }
  }

  eta_map <- c(lka = "eta_ka", lVc = "eta_Vc", lVp = "eta_Vp", lCL = "eta_CL",
               lkint = "eta_kint", lKss = "eta_Kss", lksyn = "eta_ksyn",
               lR0 = "eta_R0")
  eta_init <- c(eta_ka = 0.3225, eta_Vc = 0.3225, eta_Vp = 0.0239,
                eta_CL = 0.0680, eta_kint = 0.0062, eta_Kss = 0.2870,
                eta_ksyn = 0.0492, eta_R0 = 1.0195)

  drop_eta <- "eta_Kss"
  for (e in names(eta_init)) {
    if (e != drop_eta) {
      ini_lines <- c(ini_lines, sprintf("    %s ~ %.6f", e, eta_init[[e]]))
    }
  }
  ini_lines <- c(ini_lines,
                 "    add.err  <- 7.88",
                 "    prop.err <- 0.298")

  par_expr <- function(lname) {
    eta <- eta_map[[lname]]
    base <- sub("^l", "", lname)
    if (lname == "lKss") {
      sprintf("    %s <- exp(%s)", base, lname)
    } else {
      sprintf("    %s <- exp(%s + %s)", base, lname, eta)
    }
  }

  model_lines <- c(
    par_expr("lka"), par_expr("lVc"), par_expr("lVp"), par_expr("lCL"),
    par_expr("lkint"), par_expr("lKss"), par_expr("lksyn"), par_expr("lR0"),
    "    Q    <- 0.20",
    "    kdeg <- ksyn / R0",
    "",
    "    disc  <- (Ctot - Rtot - Kss)^2 + 4 * Kss * Ctot",
    "    discP <- max(disc, 0)",
    "    C     <- 0.5 * ((Ctot - Rtot - Kss) + sqrt(discP))",
    "    Cfree <- max(C, 0)",
    "    RC    <- Ctot - Cfree",
    "",
    "    d/dt(depot) <- -ka * depot",
    "    d/dt(Ctot)  <- (ka * depot) / Vc - (CL / Vc) * Cfree - (Q / Vc) * Cfree +",
    "                   (Q / Vc) * Cp - kint * RC",
    "    d/dt(Cp)    <- (Q / Vp) * Cfree - (Q / Vp) * Cp",
    "    d/dt(Rtot)  <- ksyn - kdeg * (Rtot - RC) - kint * RC",
    "",
    "    Rtot(0) <- R0",
    "",
    "    Ctot ~ add(add.err) + prop(prop.err)"
  )

  code <- paste0(
    "function() {\n",
    "  ini({\n", paste(ini_lines, collapse = "\n"), "\n  })\n",
    "  model({\n", paste(model_lines, collapse = "\n"), "\n  })\n",
    "}"
  )

  eval(parse(text = code))
}

cat(sprintf("=== GROUND-TRUTH RECOVERY  Kss=%.5f (TRUE)  seed=%d  nBurn=150 nEm=200  (%d rows, %d subjects, CORRECTED model, repro dataset) ===\n",
            true_kss, seed, nrow(pk_prof), length(unique(pk_prof$ID))))

t0 <- Sys.time()

status   <- "OK"
ofv      <- NA_real_
add_err  <- NA_real_
prop_err <- NA_real_
lka_v <- lVc_v <- lVp_v <- lCL_v <- lkint_v <- lksyn_v <- lR0_v <- NA_real_

fit_result <- tryCatch({
  mod <- build_model_corrected(kss_log)
  fit <- nlmixr2(mod, pk_prof, est = "saem",
                 control = saemControl(seed = seed, print = 50,
                                       nBurn = 150, nEm = 200))

  pf <- fit$parFixed
  num <- function(x) as.numeric(sub("^\\s*([0-9.eE+-]+).*$", "\\1", as.character(x)))
  bt  <- function(nm) num(pf[nm, "Back-transformed(95%CI)"])

  list(
    ofv = fit$objf,
    add_err = bt("add.err"), prop_err = bt("prop.err"),
    ka = bt("lka"), Vc = bt("lVc"), Vp = bt("lVp"), CL = bt("lCL"),
    kint = bt("lkint"), ksyn = bt("lksyn"), R0 = bt("lR0")
  )
}, error = function(e) {
  cat("FAILED:", conditionMessage(e), "\n")
  NULL
})

t1 <- Sys.time()
elapsed_sec <- as.numeric(difftime(t1, t0, units = "secs"))

if (is.null(fit_result)) {
  status <- "FAILED"
} else {
  ofv      <- fit_result$ofv
  add_err  <- fit_result$add_err
  prop_err <- fit_result$prop_err
  lka_v    <- fit_result$ka
  lVc_v    <- fit_result$Vc
  lVp_v    <- fit_result$Vp
  lCL_v    <- fit_result$CL
  lkint_v  <- fit_result$kint
  lksyn_v  <- fit_result$ksyn
  lR0_v    <- fit_result$R0
}

row <- data.frame(
  seed        = seed,
  kss_fixed   = true_kss,
  nBurn       = 150,
  nEm         = 200,
  ofv         = ofv,
  add.err     = add_err,
  prop.err    = prop_err,
  ka          = lka_v,
  Vc          = lVc_v,
  Vp          = lVp_v,
  CL          = lCL_v,
  kint        = lkint_v,
  ksyn        = lksyn_v,
  R0          = lR0_v,
  status      = status,
  elapsed_sec = round(elapsed_sec, 1),
  timestamp   = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  check.names = FALSE
)

file_existed <- file.exists(results_csv)
write.table(row, file = results_csv, sep = ",", row.names = FALSE,
            col.names = !file_existed, append = file_existed)

cat(sprintf("DONE seed=%d status=%s ofv=%s add.err=%s prop.err=%s elapsed=%.1fs\n",
            seed, status,
            ifelse(is.na(ofv), "NA", sprintf("%.4f", ofv)),
            ifelse(is.na(add_err), "NA", sprintf("%.4f", add_err)),
            ifelse(is.na(prop_err), "NA", sprintf("%.4f", prop_err)),
            elapsed_sec))
