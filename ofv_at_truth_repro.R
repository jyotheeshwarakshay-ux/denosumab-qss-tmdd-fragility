#!/usr/bin/env Rscript
# =============================================================================
# REPRODUCIBILITY CHECK: OFV at the FULLY-TRUE parameter vector (all theta,
# omega, and residual-error terms fixed at the exact values used to generate
# simulated_denosumab_data_v2_corrected_repro.csv), evaluated via FOCEi's
# EBE-only mechanism (no population parameters free -> nlmixr2 auto-switches
# to EBE estimation; no SAEM burn-in/EM loop, no crash).
#
# Usage: Rscript ofv_at_truth_repro.R <rep_index>
#
# Question this answers: is the OFV-at-truth evaluation itself deterministic
# across independent, freshly-started R processes? Given this project's
# central finding (fixed-seed SAEM non-determinism), a number this
# consequential (true params fall in the high-OFV bucket, not the low-OFV
# one) must not be accepted from a single unverified evaluation. Each repeat
# is launched as a SEPARATE Rscript process (not repeated within one R
# session), same discipline as seed_determinism_forensic.R.
#
# Isolated -- writes ONLY to ofv_at_truth_repro.csv. Does not touch any
# existing results file.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript ofv_at_truth_repro.R <rep_index>")
}
rep_index <- as.integer(args[1])
if (is.na(rep_index)) stop("rep_index must be an integer")

repo        <- path.expand("~/pharmacometrics/denosumab-tmdd-qss")
setwd(repo)
data_csv    <- file.path(repo, "simulated_denosumab_data_v2_corrected_repro.csv")
results_csv <- file.path(repo, "ofv_at_truth_repro.csv")

suppressPackageStartupMessages({
  library(nlmixr2)
  library(rxode2)
  library(dplyr)
})

pk_prof <- read.csv(data_csv)
if (length(unique(pk_prof$ID)) != 30) {
  stop("Expected 30 subjects in ", data_csv, ", found ", length(unique(pk_prof$ID)))
}

mod <- function() {
  ini({
    lka   <- fixed(log(0.0078))
    lVc   <- fixed(log(1.58))
    lVp   <- fixed(log(6.06))
    lCL   <- fixed(log(0.006))
    lkint <- fixed(log(0.022))
    lKss  <- fixed(log(1.56))
    lksyn <- fixed(log(0.01))
    lR0   <- fixed(log(15.23))
    eta_ka   ~ fixed(0.2776)
    eta_Vc   ~ fixed(0.3225)
    eta_Vp   ~ fixed(0.0239)
    eta_CL   ~ fixed(0.0673)
    eta_kint ~ fixed(0.0062)
    eta_Kss  ~ fixed(0.2910)
    eta_ksyn ~ fixed(0.0492)
    eta_R0   ~ fixed(1.2544)
    add.err  <- fixed(0.72)
    prop.err <- fixed(0.07)
  })
  model({
    ka   <- exp(lka + eta_ka)
    Vc   <- exp(lVc + eta_Vc)
    Vp   <- exp(lVp + eta_Vp)
    CL   <- exp(lCL + eta_CL)
    kint <- exp(lkint + eta_kint)
    Kss  <- exp(lKss + eta_Kss)
    ksyn <- exp(lksyn + eta_ksyn)
    R0   <- exp(lR0 + eta_R0)
    Q    <- 0.20
    kdeg <- ksyn / R0
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

# Environment capture -- same fields as seed_determinism_forensic.R, so this
# check inherits the same environment-identity rigor.
env_info <- list(
  R_version       = R.version.string,
  nlmixr2_version = as.character(packageVersion("nlmixr2")),
  rxode2_version  = as.character(packageVersion("rxode2")),
  dplyr_version   = as.character(packageVersion("dplyr")),
  os              = Sys.info()[["sysname"]],
  machine         = Sys.info()[["machine"]],
  cores_detected  = parallel::detectCores(),
  rxCores         = tryCatch(rxode2::rxCores(), error = function(e) NA),
  rxThreads       = tryCatch(rxode2::getRxThreads(), error = function(e) NA),
  RNGkind         = paste(RNGkind(), collapse = "|"),
  OMP_NUM_THREADS = Sys.getenv("OMP_NUM_THREADS"),
  OMP_THREAD_LIMIT = Sys.getenv("OMP_THREAD_LIMIT"),
  MKL_NUM_THREADS = Sys.getenv("MKL_NUM_THREADS"),
  RXODE2_THREADS  = Sys.getenv("RXODE2_THREADS")
)

cat(sprintf("=== OFV-AT-TRUTH REPRODUCIBILITY  rep=%d ===\n", rep_index))

t0 <- Sys.time()
status <- "OK"
ofv <- NA_real_
fit_result <- tryCatch({
  fit <- nlmixr2(mod, pk_prof, est = "focei",
                 control = foceiControl(maxOuterIterations = 0, print = 0))
  fit$objf
}, error = function(e) {
  cat("FAILED:", conditionMessage(e), "\n")
  NULL
})
t1 <- Sys.time()
elapsed_sec <- as.numeric(difftime(t1, t0, units = "secs"))

if (is.null(fit_result)) {
  status <- "FAILED"
} else {
  ofv <- fit_result
}

row <- data.frame(
  rep_index        = rep_index,
  ofv              = ofv,
  status           = status,
  elapsed_sec      = round(elapsed_sec, 1),
  timestamp        = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  R_version        = env_info$R_version,
  nlmixr2_version  = env_info$nlmixr2_version,
  rxode2_version   = env_info$rxode2_version,
  dplyr_version    = env_info$dplyr_version,
  os               = env_info$os,
  machine          = env_info$machine,
  cores_detected   = env_info$cores_detected,
  rxCores          = env_info$rxCores,
  rxThreads        = env_info$rxThreads,
  RNGkind          = env_info$RNGkind,
  OMP_NUM_THREADS  = env_info$OMP_NUM_THREADS,
  OMP_THREAD_LIMIT = env_info$OMP_THREAD_LIMIT,
  MKL_NUM_THREADS  = env_info$MKL_NUM_THREADS,
  RXODE2_THREADS   = env_info$RXODE2_THREADS,
  check.names = FALSE
)

file_existed <- file.exists(results_csv)
write.table(row, file = results_csv, sep = ",", row.names = FALSE,
            col.names = !file_existed, append = file_existed)

cat(sprintf("DONE rep=%d status=%s ofv=%s elapsed=%.1fs\n",
            rep_index, status,
            ifelse(is.na(ofv), "NA", sprintf("%.6f", ofv)),
            elapsed_sec))
