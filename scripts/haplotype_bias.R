##### Reciprocal class (haplotype) bias
##
##  Tests whether reciprocal phenotypic classes of backcross sons are recovered
##  in equal numbers. Each crossover class has a reciprocal class carrying the
##  opposite alleles at every marker (for example + + + + and y cv v f). Without
##  differential viability the two should be equally common. Adapted from
##  scripts/04.haplotype-bias.R in Novak et al. 2026
##  (StevisonLab/Diet-and-recombination-in-Dmel).
##
##  Run from the repository root. Reads rawdata/byvialday.csv (written by
##  scripts/01_Data_Clean-up.R).
##
##  Writes:
##    output/hap_bias_binomial_tests.csv   pooled binomial tests by treatment x background
##    output/hap_bias_glmm.csv             GLMM tests of treatment x background x brood
##    output/hap_bias_by_vialday.csv       per vial-day skew and binomial p-values
##    images/hap_bias_pooled.png           pooled skew, faceted by background
##    images/hap_bias_vialday.png          per vial-day skew (significant only), as in Novak et al. 2026

library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(car)

by_vialday <- read.csv("rawdata/byvialday.csv")
by_vialday$Stock <- factor(by_vialday$Stock, levels = c("42", "57", "391", "491", "508"))
by_vialday$Treatment <- factor(by_vialday$Treatment,
                               levels = c("Control", "ColdTemp", "Hypoxia", "Combined"))

## Reciprocal pairs. "left" is always the class with y+ (wild type at y) and
## "right" the class with y, matching the axis layout in Novak et al. 2026.
pairs <- data.frame(
  pair  = c("NCO", "SCO y-cv", "SCO cv-v", "SCO v-f",
            "DCO y-cv/cv-v", "DCO y-cv/v-f", "DCO cv-v/v-f", "TCO"),
  left  = c("wt.wt.wt.wt", "wt.cv.v.f", "wt.wt.v.f", "wt.wt.wt.f",
            "wt.cv.wt.wt", "wt.cv.v.wt", "wt.wt.v.wt", "wt.cv.wt.f"),
  right = c("y.cv.v.f", "y.wt.wt.wt", "y.cv.wt.wt", "y.cv.v.wt",
            "y.wt.v.f", "y.wt.wt.f", "y.cv.wt.f", "y.wt.v.wt"),
  stringsAsFactors = FALSE)

## Readable labels, e.g. "wt.cv.v.f" -> "+ cv v f"
hap_label <- function(x) gsub("\\.", " ", gsub("wt", "+", x))
pairs$left_label  <- hap_label(pairs$left)
pairs$right_label <- hap_label(pairs$right)

## Long format: one row per vial-day x pair
hap_long <- bind_rows(lapply(seq_len(nrow(pairs)), function(i) {
  data.frame(by_vialday[, c("vial_num", "vial_letter", "Stock", "Treatment", "F1_vial")],
             pair  = pairs$pair[i],
             left  = by_vialday[[pairs$left[i]]],
             right = by_vialday[[pairs$right[i]]])
}))
hap_long$pair <- factor(hap_long$pair, levels = pairs$pair)
hap_long$n    <- hap_long$left + hap_long$right
hap_long$skew <- ifelse(hap_long$n > 0, hap_long$right / hap_long$n, NA)   # 0.5 = no bias; > 0.5 = toward the y class
hap_long$pval <- mapply(function(l, r) if (l + r > 0) binom.test(c(r, l), p = 0.5)$p.value else NA,
                        hap_long$left, hap_long$right)

write.csv(hap_long, "output/hap_bias_by_vialday.csv", row.names = FALSE)


##### Pooled binomial tests by treatment x background -------------------------
hap_pooled <- hap_long %>%
  group_by(Treatment, Stock, pair) %>%
  summarise(left = sum(left), right = sum(right), .groups = "drop") %>%
  mutate(n    = left + right,
         skew = right / n,
         lcl  = mapply(function(r, n) if (n > 0) binom.test(r, n)$conf.int[1] else NA, right, n),
         ucl  = mapply(function(r, n) if (n > 0) binom.test(r, n)$conf.int[2] else NA, right, n),
         pval = mapply(function(l, r) if (l + r > 0) binom.test(c(r, l), p = 0.5)$p.value else NA, left, right))
hap_pooled$padj <- p.adjust(hap_pooled$pval, method = "BH")
hap_pooled <- left_join(hap_pooled, pairs[, c("pair", "left_label", "right_label")], by = "pair")

write.csv(hap_pooled, "output/hap_bias_binomial_tests.csv", row.names = FALSE)


##### Does the bias change with treatment? ------------------------------------
## Binomial GLMM per pair on vial-day counts, same fixed effects as the
## recombination models. DCO and TCO classes are too rare to model.
glmm_pairs <- c("NCO", "SCO y-cv", "SCO cv-v", "SCO v-f")
hap_glmm <- bind_rows(lapply(glmm_pairs, function(p) {
  dat <- subset(hap_long, pair == p & n > 0)
  fit <- glmer(cbind(right, left) ~ Treatment * Stock * vial_letter + (1 | F1_vial),
               data = dat, family = binomial,
               control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
  a <- as.data.frame(Anova(fit))
  data.frame(pair = p, term = rownames(a), a, row.names = NULL, check.names = FALSE)
}))

write.csv(hap_glmm, "output/hap_bias_glmm.csv", row.names = FALSE)


##### Figures -----------------------------------------------------------------
treat_cols <- c(Control = "#4d4d4d", ColdTemp = "#4393c3", Hypoxia = "#e08214", Combined = "#8073ac")

## Pooled skew with 95% binomial CI, one panel per background
pooled_plot_data <- subset(hap_pooled, n >= 10)
pooled_plot_data$y <- as.numeric(factor(pooled_plot_data$pair, levels = rev(pairs$pair))) +
  (as.numeric(pooled_plot_data$Treatment) - 2.5) * 0.18
## Significant deviations from 1:1 (BH-adjusted) are drawn solid; the rest are faded and open
pooled_plot_data$sig <- factor(ifelse(pooled_plot_data$padj < 0.05, "BH-adjusted p < 0.05", "Not significant"),
                               levels = c("BH-adjusted p < 0.05", "Not significant"))

hap_bias_pooled <- ggplot(pooled_plot_data, aes(x = skew, y = y, colour = Treatment, alpha = sig)) +
  geom_vline(xintercept = 0.5, colour = "grey40", linetype = "dashed") +
  geom_errorbar(aes(xmin = lcl, xmax = ucl), width = 0, orientation = "y") +
  geom_point(aes(fill = ifelse(sig == "Not significant", "white", as.character(Treatment))),
             shape = 21, size = 1.8, stroke = 0.6) +
  scale_fill_manual(values = c(treat_cols, white = "white"), guide = "none") +
  scale_alpha_manual(values = c("BH-adjusted p < 0.05" = 1, "Not significant" = 0.35), name = NULL) +
  scale_colour_manual(values = treat_cols,
                      labels = c("Control", "Low temperature", "Hypoxia", "Combined")) +
  guides(colour = guide_legend(override.aes = list(fill = unname(treat_cols), alpha = 1), order = 1),
         alpha  = guide_legend(override.aes = list(colour = "grey20", fill = c("grey20", "white")), order = 2)) +
  scale_y_continuous(breaks = seq_len(nrow(pairs)),
                     labels = rev(pairs$left_label),
                     sec.axis = dup_axis(labels = rev(pairs$right_label))) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  facet_wrap(~ Stock, nrow = 1, labeller = labeller(Stock = function(x) paste0("DGRP_", x))) +
  labs(x = "Proportion of pair in the y class", y = NULL, colour = NULL,
       caption = "Pairs with fewer than 10 flies in a treatment are not shown.") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

ggsave("images/hap_bias_pooled.png", hap_bias_pooled, width = 10, height = 3.6, dpi = 300)

## Per vial-day skew, significant vial-days only (layout of Novak et al. 2026, Figure S9)
vialday_sig <- subset(hap_long, !is.na(pval) & pval < 0.05)
vialday_sig$y <- as.numeric(factor(vialday_sig$pair, levels = rev(pairs$pair)))

hap_bias_vialday <- ggplot(vialday_sig, aes(x = skew, y = y)) +
  geom_vline(xintercept = 0.5, colour = "red") +
  geom_jitter(alpha = 0.3, width = 0, height = 0.15, size = 0.8) +
  scale_y_continuous(breaks = seq_len(nrow(pairs)), limits = c(0.5, nrow(pairs) + 0.5),
                     labels = rev(pairs$left_label),
                     sec.axis = dup_axis(labels = rev(pairs$right_label))) +
  scale_x_continuous(limits = c(-0.05, 1.05), breaks = c(0, 0.5, 1)) +
  facet_grid(Treatment ~ Stock,
             labeller = labeller(Stock = function(x) paste0("DGRP_", x),
                                 Treatment = c(Control = "Control", ColdTemp = "Low temp",
                                               Hypoxia = "Hypoxia", Combined = "Combined"))) +
  labs(x = "Proportion of pair in the y class", y = NULL) +
  theme_light(base_size = 8) +
  theme(panel.grid.major.y = element_blank())

ggsave("images/hap_bias_vialday.png", hap_bias_vialday, width = 10, height = 8, dpi = 300)
