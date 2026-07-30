# Panel A: Chromosome map with markers (from Novak et al 2026)
chromosome_map <- ggplot() +
  # Chromosome line
  geom_segment(aes(x = 0, xend = 70, y = 0.5, yend = 0.5), 
               size = 3, color = "gray30") +
  # Centromere indicator
  geom_point(aes(x = 68.5, y = 0.5), size = 5, shape = 16, color = "red") +
  # Marker positions
  geom_point(aes(x = c(0, 12, 33, 57), y = 0.5), 
             size = 4, color = "black") +
  geom_text(aes(x = c(0, 12, 33, 57), 
                y = 0.7,
                label = c("y (0)", "cv (12)", "v (33)", "f (57)")),
            size = 3.5) +
  # Interval brackets
  annotate("segment", x = 0, xend = 12, y = 0.3, yend = 0.3, 
           arrow = arrow(ends = "both", length = unit(0.1, "inches"))) +
  annotate("segment", x = 12, xend = 33, y = 0.25, yend = 0.25,
           arrow = arrow(ends = "both", length = unit(0.1, "inches"))) +
  annotate("segment", x = 33, xend = 57, y = 0.2, yend = 0.2,
           arrow = arrow(ends = "both", length = unit(0.1, "inches"))) +
  ylim(0, 1) +
  xlim(-5, 75) +
  labs(x = "Map Position (cM)", y = "") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank())
chromosome_map

########################################################
# Panel B: Interval-level map length by Treatment,
#          faceted by Stock (background) x Timepoint
#
# Design goals:
#  - Treatment is the "star" comparison -> x axis, so treatment
#    effects are directly comparable within every panel
#  - Stock (rows) x Timepoint (columns) faceting keeps each
#    individual panel dead simple (just 4 stacked bars),
#    so the complexity lives in the grid layout, not in any
#    single plot element
#  - Interval order (y-cv, cv-v, v-f) mirrors the chromosome
#    order from Panel A, and uses the SAME color scale, so
#    readers can connect the two panels visually
#  - Raw cM stacked bars preserve both total map length
#    (bar height) AND composition (segment sizes) at once
########################################################

## ---- 1. Reshape recom_merge to long format ----------------------------
# Assumes recom_merge has columns: Treatment, vial_letter, Stock,
# ycv_rr, cvv_rr, vf_rr (per-interval Kosambi-corrected cM)

recom_long <- recom_merge %>%
  dplyr::select(Treatment, vial_letter, Stock, ycv_rr, cvv_rr, vf_rr) %>%
  pivot_longer(cols = c(ycv_rr, cvv_rr, vf_rr),
               names_to = "interval",
               values_to = "cM") %>%
  mutate(
    interval = case_match(interval,
                      "ycv_rr" ~ "y-cv",
                      "cvv_rr" ~ "cv-v",
                      "vf_rr"  ~ "v-f"),
    interval = factor(interval, levels = c("y-cv", "cv-v", "v-f")),
    Treatment = factor(Treatment,
                       levels = c("Control", "Hypoxia", "ColdTemp", "Combined")),
    treat_order = as.numeric(Treatment)
  )

## ---- 2. Summarise: mean + SE per Treatment x Timepoint x Stock x Interval ----
# If recom_merge is already one row per individual/replicate, this collapses
# to means for plotting. If it's already summarised, this is a no-op.

recom_summary <- recom_long %>%
  group_by(Treatment, treat_order, vial_letter, Stock, interval) %>%
  summarise(
    mean_cM = mean(cM, na.rm = TRUE),
    se_cM   = sd(cM, na.rm = TRUE) / sqrt(sum(!is.na(cM))),
    .groups = "drop"
  ) %>%
  # cumulative height needed to place asterisks/outlines at the right
  # vertical position within each stacked bar
  group_by(Treatment, treat_order, vial_letter, Stock) %>%
  arrange(interval, .by_group = TRUE) %>%
  mutate(
    ymax = cumsum(mean_cM),
    ymin = ymax - mean_cM,
    y_mid = ymin + mean_cM / 2
  ) %>%
  ungroup()

## ---- 3. Extract significance from fit_ycv / fit_cvv / fit_vf ----------
# Each fit_* is a glmer binomial model with a THREE-WAY interaction:
#   cbind(count, nco_inds) ~ (1|F1_vial) + Treatment * Stock * vial_letter
#
# Because Treatment's effect can differ by Stock and by vial_letter (that's
# the whole point of fitting the interaction), a single "is Treatment
# significant" p-value doesn't exist -- what we actually want is:
# "was Treatment X different from Control, specifically within THIS
# Stock x vial_letter panel?" That's a simple-effects question, answered
# with emmeans: get estimated marginal means of Treatment within each
# Stock x vial_letter cell, then test each Treatment level against a
# reference (Control) within that cell.
#
# This also has the nice side effect of returning Stock and vial_letter
# as their own columns already, so no regex parsing of interaction term
# names is needed.

extract_sig <- function(model, interval_name, ref_level = "Control") {
  emm <- tryCatch(
    emmeans::emmeans(model, ~ Treatment | Stock + vial_letter),
    error = function(e) NULL
  )
  if (is.null(emm)) return(tibble())
  
  contr <- emmeans::contrast(emm, method = "trt.vs.ctrl", ref = ref_level)
  
  broom::tidy(contr) %>%
    mutate(
      interval = interval_name,
      # contrast labels look like "Hypoxia - Control" -- pull out the
      # non-reference treatment name
      Treatment = trimws(sub(paste0("- ", ref_level), "", contrast)),
      sig = case_when(
        adj.p.value < 0.001 ~ "***",
        adj.p.value < 0.01  ~ "**",
        adj.p.value < 0.05  ~ "*",
        TRUE ~ ""
      )
    ) %>%
    filter(sig != "") %>%
    dplyr::select(interval, Treatment, Stock, vial_letter, p.value = adj.p.value, sig)
}

extract_sig2 <- function(model, interval_name, ref_level = "Control") {
  emm <- tryCatch(
    emmeans::emmeans(model, ~ Treatment | Stock ),
    error = function(e) NULL
  )
  if (is.null(emm)) return(tibble())
  
  contr <- emmeans::contrast(emm, method = "trt.vs.ctrl", ref = ref_level)
  
  broom::tidy(contr) %>%
    mutate(
      interval = interval_name,
      # contrast labels look like "Hypoxia - Control" -- pull out the
      # non-reference treatment name
      Treatment = trimws(sub(paste0("- ", ref_level), "", contrast)),
      sig = case_when(
        adj.p.value < 0.001 ~ "***",
        adj.p.value < 0.01  ~ "**",
        adj.p.value < 0.05  ~ "*",
        TRUE ~ ""
      )
    ) %>%
    filter(sig != "") %>%
    dplyr::select(interval, Treatment, Stock, p.value = adj.p.value, sig)
}

sig_df <- bind_rows(
  extract_sig(fit_ycv, "y-cv"),
  extract_sig(fit_cvv, "cv-v"),
  extract_sig(fit_vf,  "v-f")
) %>%
  mutate(interval = factor(interval, levels = c("y-cv", "cv-v", "v-f")))

# Join significance flags onto the summary data -- now matching on
# Stock and vial_letter too, since significance is cell-specific here,
# not just a blanket "Treatment X was significant" across the board.
recom_summary <- recom_summary %>%
  left_join(sig_df %>% dplyr::select(interval, Treatment, Stock, vial_letter, sig),
            by = c("interval", "Treatment", "Stock", "vial_letter"))

## ---- 4. Shared color scale (match Panel A intervals) -------------------
interval_colors <- c("y-cv" = "#4477AA", "cv-v" = "#CCBB44", "v-f" = "#AA3377")
recom_summary$int_order=ifelse(recom_summary$interval=="y-cv",1,ifelse(recom_summary$interval=="cv-v",2,3))

treat_labels <- c("Control" = "25°C+21% Oxygen", "Hypoxia" = "8% Oxygen",
                  "ColdTemp" = "18°C", "Combined" = "18°C+8% Oxygen")

## ---- 5A. Option 1: Asterisks above significant segments ---------------
interval_map_plot_asterisk <- ggplot(recom_summary,
                                     aes(x = treat_order, y = mean_cM, fill = interval)) +
  geom_col(position = "stack", width = 0.7, color = "white", linewidth = 0.2) +
  geom_text(
    data = filter(recom_summary, !is.na(sig)),
    aes(x = treat_order, y = ymax + 0.5, label = sig),
    inherit.aes = FALSE, size = 3, vjust = 0
  ) +
   facet_grid(Stock ~ vial_letter,
              labeller = labeller(vial_letter = c(
                "A" = "1-2 DPM", "B" = "3-4 DPM",
                "C" = "5-6 DPM", "D" = "7-8 DPM"))) +
  scale_x_continuous(breaks = 1:4, labels = treat_labels) +
  scale_fill_manual(values = interval_colors, name = "Interval") +
  labs(x = "Treatment", y = "Map distance (cM)") +
    theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(size = 8),
    legend.position = "bottom"
  )
#interval_map_plot_asterisk

#Get color scales for the five stocks to plot for the three intervals at different hues!
n_treat <- 3
shade_seq <- seq(0.6, -0.3, length.out = n_treat)

col_pal_list <- lapply(col_scheme2, function(base_col) {
  lighten(base_col, amount = shade_seq)
})
names(col_pal_list) <- stock_order

col_pal_full <- unlist(col_pal_list)

# match interaction(Stock, treat_order)'s "Stock.treatnum" format exactly
names(col_pal_full) <- paste(rep(stock_order, each = n_treat), 
                             rep(1:n_treat, times = 5), sep = ".")

col_pal_full
scales::show_col(col_pal_full, ncol = n_treat)

#Wrap by stock instead of treatment:
interval_map_plot_gb <- ggplot(recom_summary,
                                     aes(x = factor(Stock,levels = stock_order), y = mean_cM, fill = interaction(Stock,int_order))) +
  geom_col(position = "stack", width = 0.7, color = "white", linewidth = 0.2) +
  geom_text(
    data = filter(recom_summary, !is.na(sig)),
    aes(x = factor(Stock,levels = stock_order), y = mean_cM+ymin/2+7, label = sig),
    inherit.aes = FALSE, size = 5, vjust = 0,color="white") +
  facet_grid(treat_order ~ vial_letter,
             labeller = labeller(vial_letter = c(
               "A" = "1-2 DPM", "B" = "3-4 DPM",
               "C" = "5-6 DPM", "D" = "7-8 DPM"),
             treat_order = c(
               "1" = "Control", "2" = "Hypoxia",
               "3" = "ColdTemp", "4" = "Combined"))) +
  #scale_x_continuous(breaks = 1:5) +
  scale_fill_manual(values = col_pal_full, breaks = paste(stock_order, 2, sep = "."), name = "Genotype",labels=stock_order) +
  labs(x = "Genotype", y = "Genetic map length (cM)") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(size = 8),
    legend.position = "bottom"
  )
interval_map_plot_gb
ggsave("images/Interval_map_byGB.png", interval_map_plot_gb,
       width = 10, height = 11, units = "in")


#Wrap by stock instead of treatment and remove day:
recom_summary2 <- recom_long %>%
  group_by(Treatment, treat_order, Stock, interval) %>%
  summarise(
    mean_cM = mean(cM, na.rm = TRUE),
    se_cM   = sd(cM, na.rm = TRUE) / sqrt(sum(!is.na(cM))),
    .groups = "drop"
  ) %>%
  # cumulative height needed to place asterisks/outlines at the right
  # vertical position within each stacked bar
  group_by(Treatment, treat_order, Stock) %>%
  arrange(interval, .by_group = TRUE) %>%
  mutate(
    ymax = cumsum(mean_cM),
    ymin = ymax - mean_cM,
    y_mid = ymin + mean_cM / 2
  ) %>%
  ungroup()

sig_df2 <- bind_rows(
  extract_sig2(fit_ycv, "y-cv"),
  extract_sig2(fit_cvv, "cv-v"),
  extract_sig2(fit_vf,  "v-f")
) %>%
  mutate(interval = factor(interval, levels = c("y-cv", "cv-v", "v-f")))

recom_summary2 <- recom_summary2 %>%
  left_join(sig_df2 %>% dplyr::select(interval, Treatment, Stock, sig),
            by = c("interval", "Treatment", "Stock"))

#Now plot new summary without Day
recom_summary2$int_order=ifelse(recom_summary2$interval=="y-cv",1,ifelse(recom_summary2$interval=="cv-v",2,3))
interval_map_gb_simple <- ggplot(recom_summary2,
                               aes(x = factor(Stock,levels = stock_order), y = mean_cM, fill = interaction(Stock,int_order))) +
  geom_col(position = "stack", width = 0.7, color = "white", linewidth = 0.2) +
  geom_text(
    data = filter(recom_summary2, !is.na(sig)),
    aes(x = factor(Stock,levels = stock_order), y = mean_cM+ymin/2+7, label = sig),
    inherit.aes = FALSE, size = 5, vjust = 0,color="white") +
  facet_wrap(~treat_order,
             labeller = labeller(treat_order = c(
                 "1" = "Control", "2" = "Hypoxia",
                 "3" = "ColdTemp", "4" = "Combined"))) +
  #scale_x_continuous(breaks = 1:5) +
  scale_fill_manual(values = col_pal_full, breaks = paste(stock_order, 2, sep = "."), name = "Genotype",labels=stock_order) +
  labs(x = "Genotype", y = "Genetic map length (cM)") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(size = 8),
    legend.position = "bottom"
  )
interval_map_gb_simple
ggsave("images/Interval_map_simple.png", interval_map_gb_simple,
       width = 10, height = 11, units = "in")

## ---- 6. View / save -----------------------------------------------------
#interval_map_plot_asterisk

#ggsave("images/Interval_map_asterisk.png", interval_map_plot_asterisk,
#       width = 10, height = 11, units = "in")


## ---- 7. Optional: combine with Panel A using patchwork ------------------
 library(patchwork)
 chromosome_map / interval_map_plot_asterisk +
   plot_layout(heights = c(1, 4)) +
   plot_annotation(tag_levels = "A")
