#!/usr/bin/env Rscript
# =============================================================================
# REPRODUCIBILITY PROFILE -- ONE fit per invocation, at ANY of the 7 Kss grid points,
# on the fresh INDEPENDENT 30-subject corrected repro dataset (seed 2025), with
# the SAME SAEM settings as the pilot (nBurn=150, nEm=200).
#
# Usage:  Rscript gate_test_N30_repro.R <grid_point_index 1-7> <seed>
#
# Question this answers: does the pilot's N=30 two-sided finding (good basin
# reliably at ~2830) reproduce on an independent N=30 draw at the same Kss=1.263
# grid point, same seeds, same SAEM settings?
#
# Does NOT touch multistart_results_corrected_full.csv, multistart_results_corrected_pilot.csv, convergence_test_N80_Kss1263.csv, or multistart_results_corrected_repro.csv.
# Writes ONLY to gate_test_N30_Kss1263_repro.csv (the SAME file the 5 Kss=1.263
# gate fits already live in -- APPENDS the other 30 fits so the final file has
# all 35 rows, 7 grid points x 5 seeds) and gate_test_N30_repro_logs/.
#
# Model equations embedded directly below, copied VERBATIM from the CORRECTED
# denosumab_QSS_TMDD_v2.R (commit 8fb40db) -- same source as
# multistart_profile_Kss_corrected_full.R, not the stale stage9 copy.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript gate_test_N30_repro.R <grid_point_index 1-7> <seed>")
}
grid_point_index <- as.integer(args[1])
seed             <- as.integer(args[2])
if (is.na(seed)) stop("seed must be an integer")

repo        <- path.expand("~/pharmacometrics/denosumab-tmdd-qss")
setwd(repo)
data_csv    <- file.path(repo, "simulated_denosumab_data_v2_corrected_repro.csv")
results_csv <- file.path(repo, "gate_test_N30_Kss1263_repro.csv")

# Kss grid -- IDENTICAL construction to every other profile in this repo:
# centre 1.263, log_range 1.5, 7 points, evenly spaced in log-space.
centre_value <- 1.263
log_range    <- 1.5
n_points     <- 7
grid <- seq(log(centre_value) - log_range,
            log(centre_value) + log_range,
            length.out = n_points)

if (is.na(grid_point_index) || grid_point_index < 1 || grid_point_index > n_points) {
  stop(sprintf("grid_point_index must be an integer between 1 and %d", n_points))
}

kss_log   <- grid[grid_point_index]
kss_value <- exp(kss_log)

# ---------------------------------------------------------------------------
# RESUME CHECK -- skip if this EXACT (grid_point_index, seed) already has a row.
# ---------------------------------------------------------------------------
if (file.exists(results_csv)) {
  existing <- tryCatch(read.csv(results_csv, stringsAsFactors = FALSE),
                        error = function(e) NULL)
  if (!is.null(existing) && nrow(existing) > 0 &&
      "grid_point_index" %in% names(existing) && "seed" %in% names(existing) &&
      any(existing$grid_point_index == grid_point_index & existing$seed == seed)) {
    cat(sprintf("SKIP grid_point=%d (Kss=%.5f) seed=%d -- already recorded in %s\n",
                grid_point_index, kss_value, seed, basename(results_csv)))
    quit(save = "no", status = 0)
  }
}

if (!file.exists(data_csv)) {
  stop("Fresh repro dataset not found: ", data_csv,
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
# MODEL BUILDER -- same structure as multistart_profile_Kss_corrected_full.R's
# build_model_corrected(), Kss fixed at 1.263, CORRECTED transfer terms.
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

cat(sprintf("=== REPRO PROFILE  grid_point=%d/%d  Kss=%.5f  seed=%d  nBurn=150 nEm=200  (%d rows, %d subjects, CORRECTED model) ===\n",
            grid_point_index, n_points, kss_value, seed, nrow(pk_prof), length(unique(pk_prof$ID))))

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
  grid_point_index = grid_point_index,
  seed        = seed,
  kss_fixed   = kss_value,
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

cat(sprintf("DONE grid_point=%d seed=%d status=%s ofv=%s add.err=%s prop.err=%s elapsed=%.1fs\n",
            grid_point_index, seed, status,
            ifelse(is.na(ofv), "NA", sprintf("%.4f", ofv)),
            ifelse(is.na(add_err), "NA", sprintf("%.4f", add_err)),
            ifelse(is.na(prop_err), "NA", sprintf("%.4f", prop_err)),
            elapsed_sec))
