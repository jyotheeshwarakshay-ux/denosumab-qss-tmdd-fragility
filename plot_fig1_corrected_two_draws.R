# =============================================================================
# FIGURE 1 REGENERATION (plot only, no caption, no paper_draft.md edit).
# Overlays two INDEPENDENT N=30 best-of-N Kss profiles, each ΔOFV relative to
# its OWN minimum (not pooled, not a shared minimum -- absolute OFVs are not
# comparable across independent random draws).
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
setwd("/Users/jyotheeshwarakshay/pharmacometrics/denosumab-tmdd-qss")

repro <- read.csv("gate_test_N30_Kss1263_repro.csv")
pilot <- read.csv("multistart_results_corrected_pilot.csv")

best_repro <- repro %>%
  group_by(grid_point_index, kss_fixed) %>%
  summarise(best_ofv = min(ofv), .groups = "drop") %>%
  arrange(kss_fixed) %>%
  mutate(delta_ofv = best_ofv - min(best_ofv), dataset = "reproducibility (N=30)")

best_pilot <- pilot %>%
  group_by(grid_point_index, kss_value) %>%
  summarise(best_ofv = min(ofv), .groups = "drop") %>%
  arrange(kss_value) %>%
  mutate(delta_ofv = best_ofv - min(best_ofv), dataset = "pilot (N=30)") %>%
  rename(kss_fixed = kss_value)

cat("=== REPRODUCIBILITY dataset (47 fits) -- best-of-N per grid point ===\n")
print(as.data.frame(best_repro[, c("kss_fixed", "best_ofv", "delta_ofv")]), row.names = FALSE)

cat("\n=== PILOT dataset (35 fits) -- best-of-N per grid point ===\n")
print(as.data.frame(best_pilot[, c("kss_fixed", "best_ofv", "delta_ofv")]), row.names = FALSE)

cat("\n=== Sanity checks ===\n")
cat("repro delta_ofv rounded to 2dp:", paste(sprintf("%.2f", best_repro$delta_ofv), collapse=", "), "\n")
cat("pilot max delta_ofv:", sprintf("%.2f", max(best_pilot$delta_ofv)), "\n")

combined <- bind_rows(
  best_repro[, c("kss_fixed", "delta_ofv", "dataset")],
  best_pilot[, c("kss_fixed", "delta_ofv", "dataset")]
)

p <- ggplot(combined, aes(x = kss_fixed, y = delta_ofv, colour = dataset, shape = dataset)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_x_log10() +
  scale_colour_manual(values = c("reproducibility (N=30)" = "#1F4E79", "pilot (N=30)" = "#B03A2E")) +
  labs(
    title = "Kss profile likelihood: two independent N=30 draws (corrected model)",
    x = "Kss (log scale)",
    y = expression(Delta*"OFV (relative to each dataset's own minimum)"),
    colour = "Dataset", shape = "Dataset"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("fig1_corrected_two_draws.png", p, width = 7, height = 5, dpi = 150)
cat("\nSaved: fig1_corrected_two_draws.png\n")
