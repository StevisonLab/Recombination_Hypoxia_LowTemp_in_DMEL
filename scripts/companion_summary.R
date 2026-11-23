##### Recombination alongside the companion study's reproductive and physiological traits
##
##  Builds one figure and one table that put the recombination results next to
##  the traits measured on the same crosses in the companion study
##  (Rivera-Rincón et al. 2026; StevisonLab/Combination-effects-of-hypoxia-and-low-temperature-in-Drosophila).
##  The figure follows Figure 4 of Novak et al. 2026: change of each treatment
##  relative to Control, for each trait and genetic background.
##
##  Trait values are rebuilt from the companion raw files in rawdata/companion/
##  with the same filters as the companion scripts, and each companion model is
##  refit here with the same formula. Because the models are full factorials of
##  oxygen x temperature (x the other factors), the four-level Treatment factor
##  gives the same fit. Each treatment is then compared with Control within each
##  background (emmeans, Dunnett adjustment), averaging over brood, stage or
##  timepoint where the model includes them. Significant differences are bold
##  in the figure.
##
##  Total map length uses the binomial GLMM of scripts/recombination_updated.R
##  with both parental classes (+ + + + and y cv v f) counted as non-crossovers.
##
##  Run from the repository root.
##
##  Writes:
##    output/companion_trait_means.csv       mean of each trait by background x treatment, with change vs Control
##    output/companion_contrasts.csv         each treatment vs Control within background (Dunnett-adjusted)
##    output/companion_model_summary.csv     p-values for oxygen, temperature, background and their interactions
##    images/companion_traits_heatmap.png

library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(glmmTMB)
library(emmeans)

treat_levels <- c("Control", "ColdTemp", "Hypoxia", "Combined")
stock_levels <- c("42", "57", "391", "491", "508")
treat_recode <- c("Normoxia + Control Temp (21% O2 + 25 C)" = "Control",
                  "Normoxia + Low Temp (21% O2 + 18 C)"     = "ColdTemp",
                  "Hypoxia + Control Temp (8% O2 + 25 C)"   = "Hypoxia",
                  "Hypoxia + Low Temp (8% O2 + 18 C)"       = "Combined")
read_companion <- function(file) read.csv(file.path("rawdata/companion", file),
                                          fileEncoding = "UTF-8-BOM", check.names = FALSE)
as_factors <- function(d) {
  d$Treatment <- factor(d$Treatment, levels = treat_levels)
  d$stock     <- factor(d$stock, levels = intersect(stock_levels, unique(as.character(d$stock))))
  d
}
bobyqa <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))


##### Key from F1 collection vial to treatment and background -----------------
## As in scripts/merge-data.R of the companion repository
f1_crosses <- read_companion("F1_Crosses.csv")
names(f1_crosses) <- c("completion_time_f1cross", "crossed_by", "date_cross", "treatment", "stock",
                       "female_virgin_vial", "male_virgin_vial", "number_females", "number_males",
                       "F1_cross_vial_num", "comments_F1")
f1_col <- read_companion("F1_Virgin_Collections.csv")
names(f1_col) <- c("completion_time_col", "collected_by", "sex", "number_flies_f1col",
                   "vial_f1col_num", "comments_f1col", "F1_cross_vial_num")
f1_key <- merge(f1_col, f1_crosses, by = "F1_cross_vial_num")
f1_key <- subset(f1_key, !(vial_f1col_num %in% 404:408))     # duplicated vial numbers
f1_key$Treatment <- unname(treat_recode[f1_key$treatment])
f1_key <- f1_key[, c("vial_f1col_num", "F1_cross_vial_num", "Treatment", "stock")]


##### Trait data (one data frame per model) -----------------------------------
## F1 fertility: offspring per F1 cross mother (Fertility.R)
f1_counts <- read_companion("F1_Fecundity_Dup_rem.csv")
names(f1_counts) <- c("completion_time_fec", "count_by", "sex", "number_flies", "vial_f1col_num",
                      "commentsFec", "F1_cross_vial_num")
f1_fert <- merge(f1_counts, f1_crosses, by = "F1_cross_vial_num")
f1_fert <- subset(f1_fert, !(comments_F1 %in% "CONTAMINATED") & !(commentsFec %in% c("CONTAMINATED", "CONT")))
f1_fert <- aggregate(number_flies ~ F1_cross_vial_num + treatment + stock + number_females, data = f1_fert, FUN = sum)
f1_fert$Treatment <- unname(treat_recode[f1_fert$treatment])
f1_fert$value <- as.integer(f1_fert$number_flies / f1_fert$number_females)
f1_fert <- as_factors(f1_fert)

## F2 fertility: offspring per backcross mother per brood, first 8 days of counts (Fertility.R)
merged <- read.csv("rawdata/merged.csv")
fer_bc <- merged[, c("vial", "day", "completion_time", "treatment", "stock", "num_offspring", "number_bc_females")]
fer_bc$completion_time <- as.Date(fer_bc$completion_time, "%m/%d/%y")
fer_bc <- aggregate(num_offspring ~ vial + day + treatment + stock + number_bc_females + completion_time,
                    data = fer_bc, FUN = sum)
fer_bc <- fer_bc %>% group_by(vial, day) %>% filter(completion_time <= min(completion_time) + 8) %>% ungroup()
f2_fert <- aggregate(num_offspring ~ vial + day + treatment + stock + number_bc_females, data = fer_bc, FUN = sum)
f2_fert$Treatment <- f2_fert$treatment
f2_fert$value <- as.integer(f2_fert$num_offspring / f2_fert$number_bc_females)
f2_fert <- as_factors(f2_fert)

## Body mass of F1 females, mg (Physiology.R)
mass <- subset(read_companion("Weight.csv")[, c("Vial_Number", "Weight", "Sex", "Treatment", "Stock")], Sex == "F" & !is.na(Weight))
mass$stock <- mass$Stock
mass$value <- mass$Weight * 1000
mass <- as_factors(mass)

## CTmax of F1 females: values <= 37 C ("looks dead") and DGRP_391 (no Control) excluded (Physiology.R)
ct <- merge(read_companion("CTmax.csv")[, c("Vial_Number", "Sex", "Ctmax_Temp")], f1_key, by.x = "Vial_Number", by.y = "vial_f1col_num")
ct <- subset(ct, Sex == "F" & Ctmax_Temp > 37 & stock != "391")
ct$value <- ct$Ctmax_Temp
ct <- as_factors(ct)

## Respiratory quotient of F1 females: RQ >= 1.1 and DGRP_391 excluded (Physiology.R).
## The volume and time corrections cancel, so RQ = VCO2 / VO2.
rq <- merge(read_companion("RQ.csv")[, c("Vial_Number", "Sex", "VCO2", "VO2")], f1_key,
            by.x = "Vial_Number", by.y = "vial_f1col_num")
rq$value <- rq$VCO2 / rq$VO2
rq <- subset(rq, Sex == "F" & value < 1.1 & stock != "391")
rq <- as_factors(rq)

## Ovaries (Ovaries.R): add missing stages as zero counts, then drop DGRP_57 and
## timepoint 2, as in the companion models
ova <- read_companion("Ovarioles.csv")
ova$Ovary_cor <- paste(ova$Ovary, ova$Vial_Number, ova$Timepoint, sep = "_")
stages <- c("S1-7", "S8-10", "S11", "S12-14", "Germarium")
ovas <- ova %>%
  distinct(Ovary_cor, .keep_all = TRUE) %>%
  dplyr::select(Ovary_cor, Ovary, Ovariole, Vial_Number, Treatment, Stock, Timepoint) %>%
  tidyr::crossing(Stage = stages) %>%
  left_join(ova[, c("Ovary_cor", "Stage", "Number_Oocyte", "TUNEL_Cell")], by = c("Ovary_cor", "Stage")) %>%
  mutate(Number_Oocyte = replace_na(Number_Oocyte, 0), TUNEL_Cell = replace_na(TUNEL_Cell, 0)) %>%
  filter(Stock != "57", Timepoint != "2")
ovas$Treatment <- unname(treat_recode[ovas$Treatment])
ovas$stock <- ovas$Stock
ovas$Timepoint <- factor(ovas$Timepoint)
ovas$Stage <- factor(ovas$Stage)
ovas$Stage_grouped <- factor(ifelse(ovas$Stage %in% c("Germarium", "S1-7"), "Early", "Late"))
ovas$TUNEL_negative <- ovas$Number_Oocyte - ovas$TUNEL_Cell
ovas <- as_factors(ovas)

## Per-ovary values for the figure: oocytes per ovariole and % TUNEL-positive oocytes
ova_by_ovary <- ovas %>%
  group_by(Ovary_cor, Treatment, stock) %>%
  summarise(oocytes = sum(Number_Oocyte), tunel = sum(TUNEL_Cell), ovarioles = first(Ovariole),
            .groups = "drop") %>%
  filter(ovarioles > 0, oocytes > 0)

## Total X map length: per vial-day counts with >= 5 sons (recombination_updated.R)
bv <- read.csv("rawdata/byvialday.csv")
sco <- c("y.wt.wt.wt", "wt.wt.wt.f", "y.cv.wt.wt", "y.cv.v.wt", "wt.wt.v.f", "wt.cv.v.f")
dco <- c("wt.cv.wt.wt", "wt.wt.v.wt", "y.wt.wt.f", "wt.cv.v.wt", "y.wt.v.f", "y.cv.wt.f")
tco <- c("y.wt.v.wt", "wt.cv.wt.f")
bv$co  <- rowSums(bv[, sco]) + 2 * rowSums(bv[, dco]) + 3 * rowSums(bv[, tco])
bv$nco <- bv$wt.wt.wt.wt + bv$y.cv.v.f
bv$num_M <- rowSums(bv[, c(sco, dco, tco, "wt.wt.wt.wt", "y.cv.v.f")])
recomb <- aggregate(cbind(co, nco, num_M) ~ vial_num + vial_letter + Stock + Treatment + F1_vial, data = bv, FUN = sum)
recomb <- subset(recomb, num_M >= 5)
recomb$stock <- recomb$Stock
recomb <- as_factors(recomb)
map_means <- read.csv("output/recombination_summary_Kosambi.csv")   # Kosambi map lengths as in Table 2


##### Models and treatment vs Control contrasts -------------------------------
fits <- list(
  "Map length (y-f)" = glmer(cbind(co, nco) ~ Treatment * stock * vial_letter + (1 | F1_vial),
                             data = recomb, family = binomial,
                             control = glmerControl(optimizer = "optimx", optCtrl = list(method = "nlminb"))),
  "F1 fertility" = glmer(value ~ Treatment * stock + (1 | F1_cross_vial_num),
                         data = f1_fert, family = poisson, control = bobyqa),
  "F2 fertility" = glmer(value ~ Treatment * stock + (1 | day) + (1 | vial),
                         data = f2_fert, family = poisson, control = bobyqa),
  "Oocytes per ovariole" = glmmTMB(Number_Oocyte ~ Stage + Treatment * stock * Timepoint + (1 | Ovary_cor),
                                   ziformula = ~ Stage, family = nbinom2, data = ovas,
                                   control = glmmTMBControl(optimizer = nlminb)),
  "TUNEL-positive oocytes" = glmer(cbind(TUNEL_Cell, TUNEL_negative) ~ Timepoint + Treatment * stock * Stage_grouped +
                                     (1 | Ovary_cor), data = ovas, family = binomial, control = bobyqa),
  "Body mass" = lm(value ~ Treatment * stock, data = mass),
  "RQ"        = lmer(value ~ Treatment * stock + (1 | F1_cross_vial_num), data = rq),
  "CTmax"     = lmer(value ~ Treatment * stock + (1 | F1_cross_vial_num), data = ct))

contrasts_vs_control <- bind_rows(lapply(names(fits), function(t) {
  emm <- emmeans(fits[[t]], ~ Treatment | stock)
  ctr <- as.data.frame(contrast(emm, method = "trt.vs.ctrl", ref = 1))   # Dunnett within background
  data.frame(trait = t, stock = as.character(ctr$stock),
             Treatment = sub(" - Control", "", ctr$contrast),
             estimate = ctr$estimate, SE = ctr$SE, p.value = ctr$p.value)
}))
write.csv(contrasts_vs_control, "output/companion_contrasts.csv", row.names = FALSE)


##### Means and change vs Control ---------------------------------------------
mean_table <- function(d, trait, value = "value") {
  d %>% group_by(stock, Treatment) %>%
    summarise(mean = mean(.data[[value]], na.rm = TRUE), sd = sd(.data[[value]], na.rm = TRUE),
              n = sum(!is.na(.data[[value]])), .groups = "drop") %>%
    mutate(trait = trait)
}
trait_means <- bind_rows(
  data.frame(stock = as.character(map_means$Stock), Treatment = map_means$Treatment,
             mean = map_means$recomb_rate, trait = "Map length (y-f)"),
  mean_table(f1_fert, "F1 fertility"),
  mean_table(f2_fert, "F2 fertility"),
  mean_table(mutate(ova_by_ovary, value = oocytes / ovarioles), "Oocytes per ovariole"),
  mean_table(mutate(ova_by_ovary, value = 100 * tunel / oocytes), "TUNEL-positive oocytes"),
  mean_table(mass, "Body mass"),
  mean_table(rq, "RQ"),
  mean_table(ct, "CTmax")) %>%
  mutate(stock = as.character(stock), Treatment = as.character(Treatment)) %>%
  group_by(trait, stock) %>%
  mutate(control_mean = mean[Treatment == "Control"][1],
         log2FC = log2(mean / control_mean),
         pct_change = 100 * (mean / control_mean - 1)) %>%
  ungroup() %>%
  left_join(contrasts_vs_control[, c("trait", "stock", "Treatment", "p.value")],
            by = c("trait", "stock", "Treatment"))
write.csv(trait_means, "output/companion_trait_means.csv", row.names = FALSE)


##### Figure ------------------------------------------------------------------
trait_order <- c("Map length (y-f)",
                 "F1 fertility", "F2 fertility", "Oocytes per ovariole", "TUNEL-positive oocytes",
                 "Body mass", "RQ", "CTmax")
heat <- subset(trait_means, Treatment != "Control")
heat <- complete(heat, trait = trait_order, stock = stock_levels,
                 Treatment = c("ColdTemp", "Hypoxia", "Combined"))
heat$trait     <- factor(heat$trait, levels = rev(trait_order))
heat$stock     <- factor(heat$stock, levels = stock_levels)
heat$Treatment <- factor(heat$Treatment, levels = c("Hypoxia", "ColdTemp", "Combined"),
                         labels = c("Hypoxia", "Low temperature", "Combined"))
heat$label <- ifelse(is.na(heat$pct_change), "",
                     paste0(ifelse(heat$pct_change > 0, "+", ""), round(heat$pct_change), "%"))
heat$sig <- !is.na(heat$p.value) & heat$p.value < 0.05

## Blue = increase, orange = decrease (as Reviewer 2 suggested);
## boxed and bold = differs from Control
companion_heatmap <- ggplot(heat, aes(x = stock, y = trait, fill = log2FC)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_tile(data = subset(heat, sig), fill = NA, colour = "black", linewidth = 0.5,
            width = 0.88, height = 0.88) +
  geom_text(aes(label = label, fontface = ifelse(sig, "bold", "plain")), size = 2.3) +
  geom_hline(yintercept = length(trait_order) - 0.5, linewidth = 0.4) +
  scale_fill_gradient2(low = "#b35806", mid = "#f7f7f7", high = "#2166ac", midpoint = 0,
                       limits = c(-1, 1), oob = scales::squish, na.value = "grey85",
                       name = expression(log[2]*"FC vs Control")) +
  facet_wrap(~ Treatment, nrow = 1) +
  labs(x = "DGRP line", y = NULL,
       caption = paste("Boxed: differs from Control within that line (Dunnett-adjusted p < 0.05).",
                       "Grey: not analysed in the companion study.")) +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(),
        legend.position = "bottom", plot.caption = element_text(hjust = 0))

ggsave("images/companion_traits_heatmap.png", companion_heatmap, width = 8, height = 3.8, dpi = 300)


##### Table: model results ----------------------------------------------------
## Keep the oxygen, temperature and background terms. Oocyte and TUNEL models
## also include stage and timepoint terms, which are left out here.
model_files <- c("F1 fertility" = "F1_fertility_anova.csv", "F2 fertility" = "F2_fertility_anova.csv",
                 "Body mass" = "Body_Mass_Anova.csv", "CTmax" = "CTmax_Anova.csv", "RQ" = "RQ_Anova.csv",
                 "Oocytes per ovariole" = "OocyteCount_anova.csv",
                 "TUNEL-positive oocytes" = "TUNEL_anova.csv")
read_model <- function(trait, path) {
  m <- read.csv(path, check.names = FALSE)
  names(m)[1] <- "term"
  data.frame(trait = trait, term = m$term, p = m[[grep("^Pr", names(m))]])
}
models <- bind_rows(
  read_model("Map length (y-f)", "output/recomb-stats_anova-full_splitTreatment.csv"),
  bind_rows(lapply(names(model_files), function(t)
    read_model(t, file.path("rawdata/companion/model_tables", model_files[[t]])))))

## Standardize term names: lower case, fixed order within interactions
term_order <- c("oxygen", "temperature", "stock")
models$term <- sapply(strsplit(tolower(models$term), ":"), function(x) {
  if (!all(x %in% term_order)) return(NA)
  paste(term_order[term_order %in% x], collapse = " x ")
})
model_summary <- models %>%
  filter(!is.na(term)) %>%
  mutate(p = signif(p, 3)) %>%
  pivot_wider(names_from = term, values_from = p)
model_summary <- model_summary[match(trait_order, model_summary$trait), ]
write.csv(model_summary, "output/companion_model_summary.csv", row.names = FALSE)
