#!/usr/bin/env Rscript
# =============================================================================
# N=80 WARM-START DIAGNOSTIC -- ONE fit per invocation, at Kss FIXED = 1.263,
# on the N=80 CORRECTED dataset, warm-started from a KNOWN N=30 good-basin
# parameter vector (not default inits), with LONG chains (nBurn=500, nEm=1000).
#
# Usage:  Rscript warmstart_N80_Kss1263.R <source_seed>
#   <source_seed> selects WHICH good-basin repro vector to start from: one of
#   42, 999, 2024 -- these are read directly from gate_test_N30_Kss1263_repro.csv
#   (grid_point_index==4, seed==<source_seed>) at runtime, not retyped, so the
#   starting vector is guaranteed to match the actual recorded good-basin fit.
#
# Question this answers: does the ~2830-3000 good basin exist at N=80 (fit
# stays near it when started there), or does N=80 genuinely prefer the ~6840
# regime (fit drifts away even from a good-basin start)?
#
# Isolated test. Does NOT touch multistart_results_corrected_full.csv,
# multistart_results_corrected_pilot.csv, convergence_test_N80_Kss1263.csv,
# gate_test_N30_Kss1263_repro.csv (READ only, never written), or
# multistart_results_corrected_repro.csv. Writes ONLY to
# warmstart_N80_Kss1263.csv (new file) and warmstart_N80_Kss1263_logs/.
#
# Model equations embedded directly below, copied VERBATIM from the CORRECTED
# denosumab_QSS_TMDD_v2.R (commit 8fb40db) -- same source as every other
# corrected worker in this repo, not the stale stage9 copy.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript warmstart_N80_Kss1263.R <source_seed>  (one of 42, 999, 2024)")
}
source_seed <- as.integer(args[1])
if (is.na(source_seed) || !(source_seed %in% c(42, 999, 2024))) {
  stop("source_seed must be one of: 42, 999, 2024")
}

repo         <- path.expand("~/pharmacometrics/denosumab-tmdd-qss")
setwd(repo)
data_csv     <- file.path(repo, "simulated_denosumab_data_v2_corrected_full.csv")
source_csv   <- file.path(repo, "gate_test_N30_Kss1263_repro.csv")
results_csv  <- file.path(repo, "warmstart_N80_Kss1263.csv")

kss_log   <- log(1.263)
kss_value <- 1.263

# ---------------------------------------------------------------------------
# RESUME CHECK -- skip if this source_seed already has a row.
# ---------------------------------------------------------------------------
if (file.exists(results_csv)) {
  existing <- tryCatch(read.csv(results_csv, stringsAsFactors = FALSE),
                        error = function(e) NULL)
  if (!is.null(existing) && nrow(existing) > 0 && "source_seed" %in% names(existing) &&
      any(existing$source_seed == source_seed)) {
    cat(sprintf("SKIP source_seed=%d -- already recorded in %s\n", source_seed, basename(results_csv)))
    quit(save = "no", status = 0)
  }
}

if (!file.exists(data_csv)) {
  stop("N=80 corrected dataset not found: ", data_csv)
}
if (!file.exists(source_csv)) {
  stop("Repro results (source of warm-start vectors) not found: ", source_csv)
}

suppressPackageStartupMessages({
  library(nlmixr2)
  library(rxode2)
  library(dplyr)
})

# ---------------------------------------------------------------------------
# PULL THE WARM-START VECTOR -- read directly from the actual recorded repro
# results, grid_point_index==4 (Kss=1.263), seed==source_seed. NOT retyped.
# ---------------------------------------------------------------------------
src <- read.csv(source_csv)
src_row <- src[src$grid_point_index == 4 & src$seed == source_seed,]
if (nrow(src_row) != 1) {
  stop("Expected exactly 1 matching row in ", source_csv, " for grid_point_index=4, seed=",
       source_seed, ", found ", nrow(src_row))
}

cat(sprintf("=== WARM-START VECTOR (from repro seed=%d, grid_point=4, Kss=1.263, OFV=%.4f) ===\n",
            source_seed, src_row$ofv))
cat(sprintf("  ka=%.6f Vc=%.4f Vp=%.4f CL=%.6f kint=%.6f ksyn=%.6f R0=%.4f add.err=%.4f prop.err=%.4f\n",
            src_row$ka, src_row$Vc, src_row$Vp, src_row$CL, src_row$kint, src_row$ksyn, src_row$R0,
            src_row$add.err, src_row$prop.err))

pk_full <- read.csv(data_csv)
if (length(unique(pk_full$ID)) != 80) {
  stop("Expected 80 subjects in ", data_csv, ", found ", length(unique(pk_full$ID)))
}

# ---------------------------------------------------------------------------
# MODEL BUILDER -- ini() starting values set to the WARM-START VECTOR (not
# default inits), Kss fixed at 1.263, CORRECTED transfer terms.
# ---------------------------------------------------------------------------
build_model_warmstart <- function(fixed_kss_log, ws) {

  init <- list(
    lka   = log(ws$ka),
    lVc   = log(ws$Vc),
    lVp   = log(ws$Vp),
    lCL   = log(ws$CL),
    lkint = log(ws$kint),
    lKss  = fixed_kss_log,
    lksyn = log(ws$ksyn),
    lR0   = log(ws$R0)
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
  # Same standard eta-variance starting guesses used throughout this repo's
  # corrected-model workers -- the warm-start is on the FIXED EFFECTS and
  # residual error (what the repro CSV actually recorded), not on omega,
  # which was never saved anywhere.
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
                 sprintf("    add.err  <- %.10f", ws$add.err),
                 sprintf("    prop.err <- %.10f", ws$prop.err))

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

cat(sprintf("\n=== N=80 WARM-START FIT  Kss=%.5f (FIXED)  source_seed=%d  nBurn=500 nEm=1000  (%d rows, %d subjects, CORRECTED model) ===\n",
            kss_value, source_seed, nrow(pk_full), length(unique(pk_full$ID))))

t0 <- Sys.time()

status   <- "OK"
ofv      <- NA_real_
add_err  <- NA_real_
prop_err <- NA_real_
lka_v <- lVc_v <- lVp_v <- lCL_v <- lkint_v <- lksyn_v <- lR0_v <- NA_real_

fit_result <- tryCatch({
  mod <- build_model_warmstart(kss_log, src_row)
  fit <- nlmixr2(mod, pk_full, est = "saem",
                 control = saemControl(seed = source_seed, print = 50,
                                       nBurn = 500, nEm = 1000))

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
  source_seed     = source_seed,
  source_ofv      = src_row$ofv,
  kss_fixed       = kss_value,
  nBurn           = 500,
  nEm             = 1000,
  start_ka        = src_row$ka,
  start_Vc        = src_row$Vc,
  start_Vp        = src_row$Vp,
  start_CL        = src_row$CL,
  start_kint      = src_row$kint,
  start_ksyn      = src_row$ksyn,
  start_R0        = src_row$R0,
  start_add.err   = src_row$add.err,
  start_prop.err  = src_row$prop.err,
  ofv             = ofv,
  add.err         = add_err,
  prop.err        = prop_err,
  ka              = lka_v,
  Vc              = lVc_v,
  Vp              = lVp_v,
  CL              = lCL_v,
  kint            = lkint_v,
  ksyn            = lksyn_v,
  R0              = lR0_v,
  status          = status,
  elapsed_sec     = round(elapsed_sec, 1),
  timestamp       = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  check.names     = FALSE
)

file_existed <- file.exists(results_csv)
write.table(row, file = results_csv, sep = ",", row.names = FALSE,
            col.names = !file_existed, append = file_existed)

cat(sprintf("DONE source_seed=%d status=%s ofv=%s add.err=%s prop.err=%s elapsed=%.1fs\n",
            source_seed, status,
            ifelse(is.na(ofv), "NA", sprintf("%.4f", ofv)),
            ifelse(is.na(add_err), "NA", sprintf("%.4f", add_err)),
            ifelse(is.na(prop_err), "NA", sprintf("%.4f", prop_err)),
            elapsed_sec))
