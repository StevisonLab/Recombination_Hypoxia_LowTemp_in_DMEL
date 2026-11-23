##### Figure S1: triple-crossover (TCO) rate by observer
##
##  Redrawn for the revision (Reviewer 1 asked for the mean and SD per
##  observer). For each anonymized observer, the point is the mean TCO rate
##  across the vial-days they scored and the bar is +/- 1 SD; the open diamond
##  is the pooled rate (all TCO sons / all sons) that the exclusion rule uses.
##  Reference lines: expected rate from FlyBase map distances, and the
##  exclusion cutoff (mean + MAD of the observers' pooled rates), as in
##  scripts/01_Data_Clean-up.R. Excluded observers are shown in red.
##
##  Run from the repository root. Reads rawdata/merged.csv.
##  Writes images/FigureS1_observer_TCO.png / .pdf

library(dplyr)
library(ggplot2)

haplotypes <- c("wt.wt.wt.wt", "y.cv.v.f", "y.wt.wt.wt", "wt.cv.wt.wt", "wt.wt.v.wt", "wt.wt.wt.f",
                "y.cv.wt.wt", "y.cv.v.wt", "y.wt.v.wt", "y.wt.wt.f", "wt.wt.v.f", "wt.cv.v.f",
                "wt.cv.v.wt", "wt.cv.wt.f", "y.wt.v.f", "y.cv.wt.f")
merged <- read.csv("rawdata/merged.csv")
merged$observer <- as.numeric(as.factor(merged$observer))     # same anonymization as 01_Data_Clean-up.R
merged$sons <- rowSums(merged[, haplotypes])
merged$tco  <- merged$y.wt.v.wt + merged$wt.cv.wt.f

per_obs <- merged %>%
  filter(sons > 0) %>%
  group_by(observer) %>%
  summarise(mean_rate = mean(tco / sons), sd_rate = sd(tco / sons),
            pooled = sum(tco) / sum(sons), sons = sum(sons), vial_days = n(), .groups = "drop")

exp_tco <- 0.12 * 0.21 * 0.24                       # expected TCO rate (FlyBase map distances)
cutoff  <- mean(per_obs$pooled) + mad(per_obs$pooled)
per_obs$status <- ifelse(per_obs$pooled > cutoff, "Excluded (TCO rate)",
                         ifelse(per_obs$sons < 50, "Excluded (< 50 sons)", "Retained"))

figure_s1 <- ggplot(per_obs, aes(x = factor(observer), colour = status)) +
  geom_hline(yintercept = 100 * exp_tco, linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = 100 * cutoff, linetype = "dotted", colour = "grey40") +
  annotate("text", x = 5.5, y = 100 * exp_tco, label = "Expected", hjust = 0.5, vjust = -0.4, size = 3, colour = "grey30") +
  annotate("text", x = 5.5, y = 100 * cutoff, label = "Exclusion cutoff", hjust = 0.5, vjust = -0.4, size = 3, colour = "grey30") +
  geom_errorbar(aes(ymin = pmax(0, 100 * (mean_rate - sd_rate)), ymax = 100 * (mean_rate + sd_rate)), width = 0.15) +
  geom_point(aes(y = 100 * mean_rate), size = 2.5) +
  geom_point(aes(y = 100 * pooled), shape = 5, size = 2.5) +
  scale_x_discrete() +
  scale_colour_manual(values = c("Retained" = "grey20", "Excluded (TCO rate)" = "#c0392b",
                                 "Excluded (< 50 sons)" = "#e67e22"), name = NULL) +
  labs(x = "Anonymized observer", y = "Triple-crossover rate (%)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), legend.position = "bottom")

ggsave("images/FigureS1_observer_TCO.png", figure_s1, width = 7, height = 4.5, dpi = 600)
ggsave("images/FigureS1_observer_TCO.pdf", figure_s1, width = 7, height = 4.5)
