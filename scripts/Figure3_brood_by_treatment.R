##### Figure 3: total map length by brood, treatment and genetic background
##
##  Redrawn for the revision (Reviewer 1): brood (days post-mating) on the
##  x-axis, one panel per treatment, genetic backgrounds by color (lines connect
##  broods). Points are observed Kosambi map lengths (mean +/- SE across
##  backcross vials, calculated as for Table 2). Significance comes from the
##  total map length GLMM (Treatment x Stock x Brood, F1 cross as a random
##  effect): error bars are thick and dark where a background differs from
##  Control in that brood (Dunnett-adjusted).
##
##  Uses fit_full from scripts/recombination_updated.R if it exists; otherwise
##  refits the same model from rawdata/byvialday.csv. Run from the repository root.
##
##  Writes:
##    images/Figure3_brood_by_treatment.png / .pdf
##    output/figure3_emmeans.csv, output/figure3_contrasts.csv

library(dplyr)
library(ggplot2)
library(lme4)
library(emmeans)

## Brood labels, used in all figures (edit here if the day ranges change)
brood_labels <- c(A = "1–2", B = "3–4", C = "5–6", D = "7–8")
treat_levels <- c("Control", "Hypoxia", "ColdTemp", "Combined")
treat_labels <- c(Control = "Control", Hypoxia = "Hypoxia", ColdTemp = "Low temperature", Combined = "Combined")
treat_cols   <- c(Control = "#4d4d4d", Hypoxia = "#e08214", ColdTemp = "#4393c3", Combined = "#8073ac")
stock_levels <- c("42", "57", "391", "491", "508")

if (!exists("fit_full")) {
  bv <- read.csv("rawdata/byvialday.csv")
  sco <- c("y.wt.wt.wt", "wt.wt.wt.f", "y.cv.wt.wt", "y.cv.v.wt", "wt.wt.v.f", "wt.cv.v.f")
  dco <- c("wt.cv.wt.wt", "wt.wt.v.wt", "y.wt.wt.f", "wt.cv.v.wt", "y.wt.v.f", "y.cv.wt.f")
  tco <- c("y.wt.v.wt", "wt.cv.wt.f")
  bv$co_inds  <- rowSums(bv[, sco]) + 2 * rowSums(bv[, dco]) + 3 * rowSums(bv[, tco])
  bv$nco_inds <- bv$wt.wt.wt.wt + bv$y.cv.v.f
  bv$num_M    <- rowSums(bv[, c(sco, dco, tco, "wt.wt.wt.wt", "y.cv.v.f")])
  recomb_rate <- aggregate(cbind(co_inds, nco_inds, num_M) ~ vial_num + vial_letter + Stock + Treatment + F1_vial,
                           data = bv, FUN = sum)
  names(recomb_rate)[names(recomb_rate) %in% c("co_inds", "nco_inds", "num_M")] <- c("co_inds.sum", "nco_inds.sum", "num_M.sum")
  recomb_rate <- subset(recomb_rate, num_M.sum >= 5)
  recomb_rate$Stock <- factor(recomb_rate$Stock)
  fit_full <- glmer(cbind(co_inds.sum, nco_inds.sum) ~ (1 | F1_vial) + Treatment * Stock * vial_letter,
                    data = recomb_rate, family = binomial(link = "logit"),
                    control = glmerControl(optimizer = "optimx", optCtrl = list(method = "nlminb")))
}

## Estimated marginal means (response scale x 100 = crossovers per 100 sons)
emm_line <- emmeans(fit_full, ~ Treatment | Stock * vial_letter, type = "response")
emm_all  <- emmeans(fit_full, ~ Treatment | vial_letter, type = "response")   # averaged over backgrounds
fig_data <- bind_rows(
  as.data.frame(emm_line) %>% mutate(Stock = as.character(Stock)),
  as.data.frame(emm_all) %>% mutate(Stock = "All lines")) %>%
  mutate(across(any_of(c("prob", "asymp.LCL", "asymp.UCL")), ~ .x * 100))

## Treatment vs Control within each background x brood (Dunnett)
ctr_line <- as.data.frame(contrast(emmeans(fit_full, ~ Treatment | Stock * vial_letter),
                                   method = "trt.vs.ctrl", ref = "Control")) %>%
  mutate(Stock = as.character(Stock))
ctr_all  <- as.data.frame(contrast(emmeans(fit_full, ~ Treatment | vial_letter),
                                   method = "trt.vs.ctrl", ref = "Control")) %>%
  mutate(Stock = "All lines")
contrasts_fig3 <- bind_rows(ctr_line, ctr_all) %>%
  mutate(Treatment = sub(" - Control", "", contrast),
         sig = case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**", p.value < 0.05 ~ "*", TRUE ~ ""))

write.csv(fig_data, "output/figure3_emmeans.csv", row.names = FALSE)
write.csv(contrasts_fig3, "output/figure3_contrasts.csv", row.names = FALSE)

##### Observed map lengths, pooled as in Figure 4 (COI.R) -----------------
## Counts are summed across backcross vials (vial-broods with >= 5 sons) for
## each background x treatment x brood, and each interval is Kosambi-corrected
## and summed (y-cv + cv-v + v-f). Pooling weights each vial by its number of
## sons, as the GLMM does. Error bars: bootstrap SE (1,000 resamples of
## backcross vials within each group).
kosambi <- function(r) ifelse(r < 0.5, (100 / 4) * log((1 + 2 * r) / (1 - 2 * r)), NA)
bv <- read.csv("rawdata/byvialday.csv")
hap <- c("wt.wt.wt.wt", "y.cv.v.f", "y.wt.wt.wt", "wt.cv.wt.wt", "wt.wt.v.wt", "wt.wt.wt.f",
         "y.cv.wt.wt", "y.cv.v.wt", "y.wt.v.wt", "y.wt.wt.f", "wt.wt.v.f", "wt.cv.v.f",
         "wt.cv.v.wt", "wt.cv.wt.f", "y.wt.v.f", "y.cv.wt.f")
ycv <- c("y.wt.wt.wt", "wt.cv.wt.wt", "y.wt.v.wt", "y.wt.wt.f", "wt.cv.v.f", "wt.cv.v.wt", "wt.cv.wt.f", "y.wt.v.f")
cvv <- c("wt.cv.wt.wt", "wt.wt.v.wt", "y.cv.wt.wt", "y.wt.v.wt", "wt.wt.v.f", "wt.cv.wt.f", "y.wt.v.f", "y.cv.wt.f")
vf  <- c("wt.wt.v.wt", "wt.wt.wt.f", "y.cv.v.wt", "y.wt.v.wt", "y.wt.wt.f", "wt.cv.v.wt", "wt.cv.wt.f", "y.cv.wt.f")
bv$num_M <- rowSums(bv[, hap]); bv$ycv <- rowSums(bv[, ycv]); bv$cvv <- rowSums(bv[, cvv]); bv$vf <- rowSums(bv[, vf])
per_vial <- aggregate(cbind(num_M, ycv, cvv, vf) ~ vial_num + vial_letter + Stock + Treatment + F1_vial, data = bv, FUN = sum)
per_vial <- subset(per_vial, num_M >= 5)
per_vial$map_cM <- rowSums(cbind(kosambi(per_vial$ycv / per_vial$num_M), kosambi(per_vial$cvv / per_vial$num_M),
                                 kosambi(per_vial$vf / per_vial$num_M)), na.rm = TRUE)
pooled_map <- function(d) sum(kosambi(sum(d$ycv) / sum(d$num_M)), kosambi(sum(d$cvv) / sum(d$num_M)),
                              kosambi(sum(d$vf) / sum(d$num_M)), na.rm = TRUE)
set.seed(2026)
obs <- per_vial %>%
  group_by(Stock = as.character(Stock), Treatment, vial_letter) %>%
  group_modify(~ data.frame(mean_cM = pooled_map(.x),
                            se_cM = sd(replicate(1000, pooled_map(.x[sample(nrow(.x), replace = TRUE), ]))),
                            n_vials = nrow(.x), sons = sum(.x$num_M))) %>%
  ungroup() %>%
  left_join(contrasts_fig3[, c("Stock", "vial_letter", "Treatment", "p.value", "sig")],
            by = c("Stock", "vial_letter", "Treatment"))
write.csv(obs, "output/figure3_observed_map_length.csv", row.names = FALSE)

##### Figure -------------------------------------------------------------------
## Lines follow each DGRP line across broods. Error bars are +/- 1 SE across
## backcross vials; they are drawn in full color where that line differs from
## Control in that brood (GLMM contrast, Dunnett-adjusted p < 0.05) and faded
## otherwise.
stock_order <- c("508", "42", "57", "391", "491")
col_scheme2 <- c("#1b9e77", "#d95f02", "#7570b3", "#66a61e", "#e7298a")
obs <- obs %>%
  mutate(sig_pt = !is.na(sig) & sig != "",
         is_control = Treatment == "Control",
         Treatment = factor(Treatment, levels = treat_levels, labels = treat_labels[treat_levels]),
         Stock = factor(Stock, levels = stock_order))

pd <- position_dodge(width = 0.35)
figure3 <- ggplot(obs, aes(x = vial_letter, y = mean_cM, colour = Stock, group = Stock)) +
  geom_line(alpha = 0.75, linewidth = 1, position = pd) +
  geom_errorbar(aes(ymin = mean_cM - se_cM, ymax = mean_cM + se_cM, alpha = sig_pt),
                width = 0.2, linewidth = 0.7, position = pd) +
  geom_point(alpha = 0.75, size = 1.6, position = pd) +
  facet_wrap(~ Treatment, nrow = 1) +
  scale_x_discrete(labels = brood_labels) +
  scale_colour_manual(values = col_scheme2) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.3), guide = "none") +
  labs(x = "Days post-mating", y = "Genetic map length (cM)", colour = NULL) +
  theme_minimal(base_size = 10) +
  theme(strip.text = element_text(size = 10, face = "bold"),
        panel.grid.minor = element_blank(),
        legend.position = "top")

ggsave("images/Figure3_brood_by_treatment.png", figure3, width = 8, height = 4, dpi = 600)
ggsave("images/Figure3_brood_by_treatment.pdf", figure3, width = 8, height = 4)
