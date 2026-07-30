##### Load required libraries
library(lmerTest)
library(dplyr)
library(reshape2)
library(ggplot2)
library(emmeans)
library(lme4)
library(car)
library(ggthemes)
library(doBy)
library(multcomp)
library(optimx)
library(readr)
library(doBy)
library(patchwork)
library(colorspace)
library(gghighlight)
library(ggpubr)
library(tidyr)
library(broom)         
library(broom.mixed)

##### RUN ANALYSIS
##
### Load custom functions
## Assign stars based on level of significance
eval.significance <- function(input) {
  ifelse(input$p.value < 0.001, "***",
         ifelse(input$p.value < 0.01, "**",
                ifelse(input$p.value < 0.05, "*", "")))
}

#Run Data clean-up:
source("scripts/01_Data_Clean-up.R")

# Examine TCO Plot:
tco_check

#Run recombination analysis:
source("scripts/recombination_updated.R")

#Numbers for main text:
sum(recomb$num_M)
sum(recomb$num_offspring)
length(unique(sort(recomb$F1_vial)))
