#!/usr/bin/env Rscript
# =============================================================================
# ONE (grid_point, seed) FIT -- FULL worker for the CORRECTED Kss profile (80 subj).
#
# Usage:  Rscript multistart_profile_Kss_corrected_full.R <grid_point_index 1-7> <seed>
#
# DIFFERS FROM multistart_profile_Kss.R:
#   - Model equations are embedded DIRECTLY below (verbatim from the CORRECTED
#     denosumab_QSS_TMDD_v2.R, commit 8fb40db), NOT extracted from
#     stage9_profile_likelihood.R's build_model() -- that file's embedded copy
#     still has the OLD BUGGY transfer terms and was never touched by the fix.
#   - Reads simulated_denosumab_data_v2_corrected_full.csv (fresh 80-subject
#     data from the corrected model), NOT simulated_denosumab_data_v2.csv (the
#     old superseded, buggy-model dataset).
#   - Writes to multistart_results_corrected_full.csv, a NEW file -- does not
#     touch or append to multistart_results.csv (the superseded results).
#
# Fits build_model_corrected("lKss", <grid value>) with SAEM (nBurn=150,
# nEm=200) on the fresh 80-subject full dataset, and appends ONE row to
# multistart_results_corrected_full.csv:
#   grid_point_index, kss_value, seed, ofv, add.err, prop.err, status,
#   elapsed_sec, timestamp
#
# CRASH-SAFETY / RESUMABLE: same pattern as multistart_profile_Kss.R -- checks
# the results CSV for an existing (grid_point_index, seed) row and skips if
# already recorded.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript multistart_profile_Kss_corrected_full.R <grid_point_index 1-7> <seed>")
}
grid_point_index <- as.integer(args[1])
seed             <- as.integer(args[2])

repo        <- path.expand("~/pharmacometrics/denosumab-tmdd-qss")
setwd(repo)
data_csv    <- file.path(repo, "simulated_denosumab_data_v2_corrected_full.csv")
results_csv <- file.path(repo, "multistart_results_corrected_full.csv")

# ---------------------------------------------------------------------------
# Kss grid -- IDENTICAL to the original: centre 1.263, log_range 1.5, 7 points,
# evenly spaced in log-space. (Unchanged -- the grid itself isn't part of the
# bug; only the ODEs were wrong.)
# ---------------------------------------------------------------------------
centre_value <- 1.263
log_range    <- 1.5
n_points     <- 7
grid <- seq(log(centre_value) - log_range,
            log(centre_value) + log_range,
            length.out = n_points)

if (is.na(grid_point_index) || grid_point_index < 1 || grid_point_index > n_points) {
  stop(sprintf("grid_point_index must be an integer between 1 and %d", n_points))
}
if (is.na(seed)) stop("seed must be an integer")

kss_log   <- grid[grid_point_index]
kss_value <- exp(kss_log)

# ---------------------------------------------------------------------------
# RESUME CHECK
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
  stop("Fresh full dataset not found: ", data_csv,
       " -- run generate_corrected_full_data.R first.")
}

suppressPackageStartupMessages({
  library(nlmixr2)
  library(rxode2)
  library(dplyr)
})

pk_prof <- read.csv(data_csv)
if (length(unique(pk_prof$ID)) != 80) {
  stop("Expected 80 subjects in ", data_csv, ", found ", length(unique(pk_prof$ID)))
}

# ---------------------------------------------------------------------------
# MODEL BUILDER -- CORRECTED transfer terms embedded directly.
# Structure (init values, eta map, Q fixed at 0.20, fixed() mechanism for the
# profiled parameter) copied from stage9_profile_likelihood.R's build_model();
# ONLY the d/dt(Ctot)/d/dt(Cp) lines are changed, to match denosumab_QSS_TMDD_v2.R
# commit 8fb40db exactly.
# ---------------------------------------------------------------------------
build_model_corrected <- function(which_par, fixed_value) {

  init <- list(
    lka   = log(0.00876),
    lVc   = log(2.169),
    lVp   = log(6.190),
    lCL   = log(0.00706),
    lkint = log(0.01680),
    lKss  = log(1.263),
    lksyn = log(0.00765),
    lR0   = log(15.638)
  )
  init[[which_par]] <- fixed_value

  ini_lines <- character(0)
  for (nm in names(init)) {
    if (nm == which_par) {
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

  drop_eta <- eta_map[[which_par]]
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
    if (lname == which_par) {
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

cat(sprintf("=== FULL grid_point=%d/%d  Kss=%.5f  seed=%d  (%d rows, %d subjects, CORRECTED model) ===\n",
            grid_point_index, n_points, kss_value, seed,
            nrow(pk_prof), length(unique(pk_prof$ID))))

t0 <- Sys.time()

status   <- "OK"
ofv      <- NA_real_
add_err  <- NA_real_
prop_err <- NA_real_

fit_result <- tryCatch({
  mod <- build_model_corrected("lKss", kss_log)
  fit <- nlmixr2(mod, pk_prof, est = "saem",
                 control = saemControl(seed = seed, print = 0,
                                       nBurn = 150, nEm = 200))

  ae <- fit$parFixed["add.err",  "Back-transformed(95%CI)"]
  pe <- fit$parFixed["prop.err", "Back-transformed(95%CI)"]
  num <- function(x) as.numeric(sub("^\\s*([0-9.eE+-]+).*$", "\\1", as.character(x)))

  list(ofv = fit$objf, add_err = num(ae), prop_err = num(pe))
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
}

row <- data.frame(
  grid_point_index = grid_point_index,
  kss_value        = kss_value,
  seed             = seed,
  ofv              = ofv,
  add.err          = add_err,
  prop.err         = prop_err,
  status           = status,
  elapsed_sec      = round(elapsed_sec, 1),
  timestamp        = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  check.names      = FALSE
)

file_existed <- file.exists(results_csv)
write.table(row, file = results_csv, sep = ",", row.names = FALSE,
            col.names = !file_existed, append = file_existed)

cat(sprintf("DONE grid_point=%d seed=%d status=%s ofv=%s elapsed=%.1fs\n",
            grid_point_index, seed, status,
            ifelse(is.na(ofv), "NA", sprintf("%.4f", ofv)), elapsed_sec))
