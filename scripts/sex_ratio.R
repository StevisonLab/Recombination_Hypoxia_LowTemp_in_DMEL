##### Sex ratio and sample sizes
##
##  1. Backcross progeny (the generation scored for recombination). Sons carry
##     the maternal (marked or recombinant) X and daughters carry the paternal
##     DGRP X, so marker-associated mortality lowers the proportion of sons.
##  2. F1 progeny (from the companion study). These developed under treatment,
##     and all F1 sons carry the y cv v f X, so the F1 sex ratio shows whether
##     marker-associated mortality changes with treatment.
##
##  Also rebuilds the Table 1 sample sizes from the vial-day counts, with total
##  offspring = females + scored males.
##
##  Run from the repository root. Reads rawdata/byvialday.csv and the companion
##  raw files in rawdata/companion/ (copied from
##  StevisonLab/Combination-effects-of-hypoxia-and-low-temperature-in-Drosophila).
##
##  Writes:
##    output/table1_sample_sizes.csv
##    output/sex_ratio_backcross_glmm.csv
##    output/sex_ratio_F1_glmm.csv
##    output/sex_ratio_summary.csv
##    images/sex_ratio.png

library(dplyr)
library(ggplot2)
library(lme4)
library(car)
library(patchwork)

treat_levels <- c("Control", "ColdTemp", "Hypoxia", "Combined")
treat_labels <- c(Control = "Control", ColdTemp = "Low temp", Hypoxia = "Hypoxia", Combined = "Combined")
treat_cols   <- c(Control = "#4d4d4d", ColdTemp = "#4393c3", Hypoxia = "#e08214", Combined = "#8073ac")
stock_levels <- c("42", "57", "391", "491", "508")
dgrp_label   <- function(x) paste0("DGRP_", x)

haplotypes <- c("wt.wt.wt.wt", "y.cv.v.f", "y.wt.wt.wt", "wt.cv.wt.wt", "wt.wt.v.wt", "wt.wt.wt.f",
                "y.cv.wt.wt", "y.cv.v.wt", "y.wt.v.wt", "y.wt.wt.f", "wt.wt.v.f", "wt.cv.v.f",
                "wt.cv.v.wt", "wt.cv.wt.f", "y.wt.v.f", "y.cv.wt.f")


##### 1. Backcross progeny ----------------------------------------------------
bc <- read.csv("rawdata/byvialday.csv")
bc$males   <- rowSums(bc[, haplotypes])
bc$females <- bc$num_females
bc$Stock     <- factor(bc$Stock, levels = stock_levels)
bc$Treatment <- factor(bc$Treatment, levels = treat_levels)
bc <- subset(bc, males + females > 0)

## Table 1 sample sizes. Backcrosses with fewer than five scored sons are
## dropped from the recombination models, as in scripts/recombination_updated.R.
## num_offspring in byvialday.csv does not equal females + scored males for
## most vial-days, so total offspring is recalculated here.
males_per_bc <- aggregate(males ~ vial_num, data = bc, FUN = sum)
table1 <- bc %>%
  group_by(Stock, Treatment) %>%
  summarise(F1_crosses_with_scored_progeny = n_distinct(F1_vial),
            backcrosses = n_distinct(vial_num),
            backcrosses_5plus_males = sum(unique(vial_num) %in% males_per_bc$vial_num[males_per_bc$males >= 5]),
            females = sum(females),
            males = sum(males),
            total_offspring = females + males,
            num_offspring_column = sum(num_offspring),
            prop_male = round(males / total_offspring, 3),
            .groups = "drop")
write.csv(table1, "output/table1_sample_sizes.csv", row.names = FALSE)

## Binomial GLMM: proportion of sons by treatment x background x brood, with
## backcross vial nested in F1 cross and an observation-level random effect
## for overdispersion.
bc$obs <- factor(seq_len(nrow(bc)))
fit_bc <- glmer(cbind(males, females) ~ Treatment * Stock * vial_letter +
                  (1 | F1_vial / vial_num) + (1 | obs),
                data = bc, family = binomial,
                control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
anova_bc <- as.data.frame(Anova(fit_bc))
write.csv(anova_bc, "output/sex_ratio_backcross_glmm.csv")


##### 2. F1 progeny (companion study) -----------------------------------------
## Same steps as scripts/merge-data.R and scripts/Fertility.R in the companion
## repository: drop contaminated vials, keep counts from the first seven days.
f1_crosses <- read.csv("rawdata/companion/F1_Crosses.csv", check.names = FALSE)
names(f1_crosses) <- c("completion_time_f1cross", "crossed_by", "date_cross", "treatment", "stock",
                       "female_virgin_vial", "male_virgin_vial", "number_females", "number_males",
                       "F1_cross_vial_num", "comments_F1")
f1_counts <- read.csv("rawdata/companion/F1_Fecundity_Dup_rem.csv", check.names = FALSE)
names(f1_counts) <- c("completion_time_fec", "count_by", "sex", "number_flies", "vial_f1col_num",
                      "commentsFec", "F1_cross_vial_num")

f1 <- merge(f1_counts, f1_crosses, by = "F1_cross_vial_num")
f1 <- subset(f1, !(comments_F1 %in% "CONTAMINATED") & !(commentsFec %in% c("CONTAMINATED", "CONT")))
f1$date <- as.Date(f1$completion_time_fec, "%m/%d/%y")
f1 <- f1 %>% group_by(F1_cross_vial_num) %>% filter(date <= min(date) + 7) %>% ungroup()

f1$Treatment <- dplyr::recode(f1$treatment,
                       "Normoxia + Control Temp (21% O2 + 25 C)" = "Control",
                       "Normoxia + Low Temp (21% O2 + 18 C)"     = "ColdTemp",
                       "Hypoxia + Control Temp (8% O2 + 25 C)"   = "Hypoxia",
                       "Hypoxia + Low Temp (8% O2 + 18 C)"       = "Combined")
f1_vial <- f1 %>%
  group_by(F1_cross_vial_num, Treatment, stock) %>%
  summarise(females = sum(number_flies[sex == "Female"]),
            males   = sum(number_flies[sex == "Male"]), .groups = "drop") %>%
  filter(females + males > 0)
f1_vial$Stock     <- factor(f1_vial$stock, levels = stock_levels)
f1_vial$Treatment <- factor(f1_vial$Treatment, levels = treat_levels)

fit_f1 <- glmer(cbind(males, females) ~ Treatment * Stock + (1 | F1_cross_vial_num),
                data = f1_vial, family = binomial,
                control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
anova_f1 <- as.data.frame(Anova(fit_f1))
write.csv(anova_f1, "output/sex_ratio_F1_glmm.csv")


##### Summary table and figure ------------------------------------------------
prop_ci <- function(m, f) {
  ci <- binom.test(m, m + f)$conf.int
  c(prop = m / (m + f), lcl = ci[1], ucl = ci[2])
}
summarise_sr <- function(dat, generation, by) {
  dat %>%
    group_by(across(all_of(by))) %>%
    summarise(males = sum(males), females = sum(females), .groups = "drop") %>%
    rowwise() %>%
    mutate(prop = prop_ci(males, females)["prop"],
           lcl  = prop_ci(males, females)["lcl"],
           ucl  = prop_ci(males, females)["ucl"]) %>%
    ungroup() %>%
    mutate(generation = generation)
}
sr_bc_brood <- summarise_sr(bc, "Backcross", c("Stock", "Treatment", "vial_letter"))
sr_bc       <- summarise_sr(bc, "Backcross", c("Stock", "Treatment"))
sr_f1       <- summarise_sr(f1_vial, "F1", c("Stock", "Treatment"))
write.csv(bind_rows(sr_f1, sr_bc, sr_bc_brood), "output/sex_ratio_summary.csv", row.names = FALSE)

## A: F1 sons (developed under treatment); B: backcross sons by brood
sr_main <- bind_rows(sr_f1, sr_bc)
sr_main$generation <- factor(sr_main$generation, levels = c("F1", "Backcross"),
                             labels = c("F1 (reared under treatment)", "Backcross (scored for recombination)"))
pd <- position_dodge(width = 0.6)

panel_a <- ggplot(sr_main, aes(x = Stock, y = prop, colour = Treatment)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(ymin = lcl, ymax = ucl), position = pd, size = 0.25) +
  facet_wrap(~ generation) +
  scale_x_discrete(labels = dgrp_label) +
  scale_colour_manual(values = treat_cols, labels = treat_labels) +
  labs(x = NULL, y = "Proportion of sons", colour = NULL) +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

panel_b <- ggplot(sr_bc_brood, aes(x = vial_letter, y = prop, colour = Treatment, group = Treatment)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
  geom_line(position = pd) +
  geom_pointrange(aes(ymin = lcl, ymax = ucl), position = pd, size = 0.2) +
  facet_wrap(~ Stock, nrow = 1, labeller = labeller(Stock = dgrp_label)) +
  scale_x_discrete(labels = c(A = "2-3", B = "4-5", C = "6-7", D = "8-9")) +
  scale_colour_manual(values = treat_cols, labels = treat_labels) +
  labs(x = "Days post-mating (brood)", y = "Proportion of sons", colour = NULL) +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "none")

sex_ratio_fig <- panel_a / panel_b + plot_annotation(tag_levels = "A")
ggsave("images/sex_ratio.png", sex_ratio_fig, width = 8, height = 6, dpi = 300)
