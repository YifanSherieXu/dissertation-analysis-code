## DVs:
##  - Binary: comforts, physical, verbal, object (glmer, binomial)
##  - Score:  comforting (lmer, Gaussian)
##
## Structure per DV:
##  full (interaction) -> reduced (main effects) -> null (control only)
##  + model comparisons + drop1 + summary + emmeans + rough plot
############################################################

rm(list = ls())

## 1) Working directory & data
library(bruceR)

setwd("~/Documents/comforting")
xdata <- read.csv("comforting_final.csv", header = TRUE, stringsAsFactors = TRUE)

## 2) Packages
install.packages("lmerTest")
library(lmerTest)
library(lme4)
library(emmeans)
library(ggeffects)

library(car)

## 3) Controls for optimisers (improves convergence)
ctrl_glmer <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
ctrl_lmer  <- lmerControl(optimizer = "bobyqa",  optCtrl = list(maxfun = 2e5))

## 4) Quick checks
str(xdata)
table(xdata$site, useNA = "ifany")
table(xdata$timepoint, useNA = "ifany")
table(xdata$sex, useNA = "ifany")

## 5) Ensure binary outcomes are coded consistently
xdata$comforts <- factor(xdata$comforts, levels = c("N","Y"))
xdata$physical <- factor(xdata$physical, levels = c("N","Y"))
xdata$verbal   <- factor(xdata$verbal,   levels = c("N","Y"))
xdata$object   <- factor(xdata$object,   levels = c("N","Y"))



xdata$comforts_num <- ifelse(xdata$comforts == "N", 0, 1)
xdata$physical_num <- ifelse(xdata$physical == 'N',0,1)
xdata$verbal_num <- ifelse(xdata$verbal == 'N',0,1)
xdata$object_num <- ifelse(xdata$object =='N',0,1)

table(xdata$comforts, useNA = "ifany")
table(xdata$physical, useNA = "ifany")
table(xdata$verbal,   useNA = "ifany")
table(xdata$object,   useNA = "ifany")

## 6) Precentage of different comforting behaviour
xdata %>% group_by(site, timepoint) %>%
  summarise(
    total_N = n(),
    valid_N = sum(!is.na(comforts)),
    empathy_count = sum (comforts == 'Y', na.rm=T),
    empathy_percentage = mean(comforts == 'Y', na.rm=T)*100,
    phy_empathy = sum(physical == 'Y', na.rm=T),
    phy_percentage = mean(physical == 'Y', na.rm=T)*100,
    verbal_empathy = sum(verbal == 'Y',na.rm = T),
    verbal_percentage = mean(verbal == 'Y',na.rm = T)*100,
    object_empathy = sum(object == 'Y',na.rm=T),
    object_percentage = mean(object == 'Y',na.rm = T)*100,
    conforting_mean = mean(comforting,na.rm=T),
    conforting_sd = sd(comforting, na.rm=T)
  )

############################################################
## Model 1: Comforting occurrence (binary Y/N)  [DV = comforts]
############################################################
cat("\n==================== Model 1: Comforting occurrence (comforts) ====================\n")

# Full model (with interaction)
occurrence.full <- glmer(
  comforts ~ site * timepoint + sex + (1 | participant),
  data = xdata, family = binomial, control = ctrl_glmer
)


# Reduced model (main effects only)
occurrence.reduced <- glmer(
  comforts ~ site + timepoint + sex + (1 | participant),
  data = xdata, family = binomial, control = ctrl_glmer
)


# Null model (control only)
occurrence.null <- glmer(
  comforts ~ sex + (1 | participant),
  data = xdata, family = binomial, control = ctrl_glmer
)

## Model comparisons
anova(occurrence.null, occurrence.full, test = "chisq")      # overall predictors
anova(occurrence.reduced, occurrence.full, test = "chisq")   # interaction test

## Main effects (if interaction not supported)
drop1(occurrence.reduced, test = "Chisq")
summary(occurrence.reduced)

## Post-hoc contrasts (main effects model)
emmeans(occurrence.reduced, pairwise ~ site)
emmeans(occurrence.reduced, pairwise ~ timepoint)

## Rough visualisation
plot(allEffects(occurrence.reduced), ylab = "outcome")

## Formal visualisation
library(cowplot)
library(lemon)
library(ggsci)
plot_data_site <- as.data.frame(emmeans(occurrence.reduced, ~ site, type = "response"))
plot_data_time <- as.data.frame(emmeans(occurrence.reduced, ~ timepoint, type = "response"))


p1 <- ggplot(plot_data_site, aes(x = site, y = prob,color= site)) +
  geom_jitter(data = xdata,aes(x=site, y=comforts_num, color=site),alpha=0.3)+
  scale_color_npg()+
  geom_point(size = 3) +  
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) + 
  theme_classic() +      
  labs(x = "", y = "Estimated Probability of Occurrence") +
  theme(text = element_text(size = 12))+
  scale_y_continuous(breaks = c( 0,0.25,0.5, 0.75,1))+
  coord_capped_cart(bottom='none', left='both', ylim = c(0,1)) +theme(legend.position="None")
p1

p2 <- ggplot(plot_data_time, aes(x = timepoint, y = prob, color= timepoint)) +
  geom_jitter(data = xdata,aes(x=timepoint, y=comforts_num), alpha=0.5)+
  geom_point(size = 3,color = "#2A5783") +
  geom_line(aes(group = 1), color = "#2A5783", size = 1) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, color = "#2A5783") +
  theme_classic() +
  ylim(0, 1) +       
  labs(x = "", y = "") + 
  theme(text = element_text(size = 12))+
  scale_y_continuous(breaks = c( 0,0.25,0.5, 0.75,1))+
  coord_capped_cart(bottom='none', left='both', ylim = c(0,1))+
  scale_color_manual(values = c('#ACD4EC', '#99C5E3','#7CADD2'))  +theme(legend.position="None")
p2

plot_grid(p1, p2, labels = c("A", "B"), rel_widths = c(1, 1))




############################################################
## Model 2: Physical comforting (binary Y/N)  [DV = physical]
############################################################
cat("\n==================== Model 2: Physical comforting (physical) ====================\n")

## Subset to non-missing rows for this DV
x_phys <- droplevels(subset(xdata, !is.na(physical)))
x_phys <- xdata
physical.full <- glmer(
  physical ~ site * timepoint + sex + (1 | participant),
  data = x_phys, family = binomial, control = ctrl_glmer
)

physical.reduced <- glmer(
  physical ~ site + timepoint + sex + (1 | participant),
  data = x_phys, family = binomial, control = ctrl_glmer
)

physical.null <- glmer(
  physical ~ sex + (1 | participant),
  data = x_phys, family = binomial, control = ctrl_glmer
)

anova(physical.null, physical.full, test = "chisq")
anova(physical.reduced, physical.full, test = "chisq")

drop1(physical.full, test = "Chisq")
summary(physical.full)

emmeans(physical.full, pairwise ~ site | timepoint)
emmeans(physical.full, pairwise ~ timepoint | site)

emmip(physical.full, site ~ timepoint, type = "response", CIs = TRUE) +
  theme_classic() +
  labs(
    title = "Developmental changes in physical comforting across sites",
    x = "Timepoint",
    y = "Predicted probability of physical comforting"
  )



# formal plots
emm_inter <- emmeans(physical.full, ~ site * timepoint, type = "response")

plot_data_inter <- as.data.frame(emm_inter)

p_inter <- ggplot(plot_data_inter, aes(x = timepoint, y = prob, color = site, group = site)) +
  geom_jitter(data = xdata, 
              aes(x = timepoint, y = physical_num, color = site), 
              alpha = 0.15,      
              width = 0.2,       
              height = 0.03,     
              size = 1) +        
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), 
                width = 0.15, 
                size = 1, 
                position = position_dodge(0.15)) + 
  geom_line(size = 1.2, position = position_dodge(0.15)) +
  geom_point(size = 3.5, position = position_dodge(0.15)) +
  scale_color_npg() + 
  theme_classic() +  
  coord_cartesian(ylim = c(-0.05, 1.05)) + 
  scale_y_continuous(labels = scales::percent_format(), breaks = seq(0, 1, 0.25)) +
  
  labs(x = "", 
       y = "Probability of Physical Comforting", 
       color = "Site") +
  theme(text = element_text(size = 12),
        legend.position = "bottom") +
  scale_y_continuous(breaks = c( 0,0.25,0.5, 0.75,1))+
  coord_capped_cart(bottom='none', left='both', ylim = c(0,1))

p_inter

############################################################
## Model 3: Verbal comforting (binary Y/N)  [DV = verbal]
############################################################
cat("\n==================== Model 3: Verbal comforting (verbal) ====================\n")

x_ver <- droplevels(subset(xdata, !is.na(verbal)))

verbal.full <- glmer(
  verbal ~ site * timepoint + sex + (1 | participant),
  data = x_ver, family = binomial, control = ctrl_glmer
)

verbal.reduced <- glmer(
  verbal ~ site + timepoint + sex + (1 | participant),
  data = x_ver, family = binomial, control = ctrl_glmer
)

verbal.null <- glmer(
  verbal ~ sex + (1 | participant),
  data = x_ver, family = binomial, control = ctrl_glmer
)

anova(verbal.null, verbal.full, test = "chisq")
anova(verbal.reduced, verbal.full, test = "chisq")

drop1(verbal.reduced, test = "Chisq")
summary(verbal.reduced)

emmeans(verbal.reduced, pairwise ~ site)
emmeans(verbal.reduced, pairwise ~ timepoint)

plot(allEffects(verbal.reduced), ylab = "outcome")

# formal plots
plot_data_time <- as.data.frame(emmeans(verbal.reduced, ~ timepoint, type = "response"))


ggplot(plot_data_time, aes(x = timepoint, y = prob, color= timepoint)) +
  geom_jitter(data = xdata,aes(x=timepoint, y=verbal_num), alpha=0.5)+
  geom_point(size = 3,color = "#2A5783") +
  geom_line(aes(group = 1), color = "#2A5783", size = 1) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, color = "#2A5783") +
  theme_classic() +
  ylim(0, 1) +       
  labs(x = "", y = "Probability of Verbal Comforting") + 
  theme(text = element_text(size = 12))+
  scale_y_continuous(breaks = c( 0,0.1,0.2, 0.3, 0.4))+
  coord_capped_cart(bottom='none', left='both', ylim = c(0,0.4))+
  scale_color_manual(values = c('#ACD4EC', '#99C5E3','#7CADD2'))  +theme(legend.position="None")


############################################################
## Model 4: Object-mediated comforting (binary Y/N)  [DV = object]
############################################################
cat("\n==================== Model 4: Object-mediated comforting (object) ====================\n")

x_obj <- droplevels(subset(xdata, !is.na(object)))

object.full <- glmer(
  object ~ site * timepoint + sex + (1 | participant),
  data = x_obj, family = binomial, control = ctrl_glmer
)

object.reduced <- glmer(
  object ~ site + timepoint + sex + (1 | participant),
  data = x_obj, family = binomial, control = ctrl_glmer
)

object.null <- glmer(
  object ~ sex + (1 | participant),
  data = x_obj, family = binomial, control = ctrl_glmer
)

anova(object.null, object.full, test = "chisq")
anova(object.reduced, object.full, test = "chisq")

drop1(object.reduced, test = "Chisq")
summary(object.reduced)

emmeans(object.reduced, pairwise ~ site)
emmeans(object.reduced, pairwise ~ timepoint)

plot(allEffects(object.reduced), ylab = "outcome")



############################################################
## Model 5: Comforting score (0–3)  [DV = comforting; Gaussian]
############################################################
cat("\n==================== Model 5: Comforting score (comforting) ====================\n")

## Fit models (ML for comparisons)
score.full <- lmer(
  comforting ~ site * timepoint + sex + (1 | participant),
  data = xdata, REML = FALSE, control = ctrl_lmer
)

score.reduced <- lmer(
  comforting ~ site + timepoint + sex + (1 | participant),
  data = xdata, REML = FALSE, control = ctrl_lmer
)

score.null <- lmer(
  comforting ~ sex + (1 | participant),
  data = xdata, REML = FALSE, control = ctrl_lmer
)

## LRT model comparisons (p values in Pr(>Chisq))
anova(score.null, score.full)
anova(score.reduced, score.full)

## Term tests for FULL model (LRT p values in Pr(Chi))
drop1(score.full, test = "Chisq")

## FULL model estimates
summary(score.full)

## p values for fixed effects (Type III; p in Pr(>F))
library(lmerTest)
anova(score.full, type = 3)

## Post-hoc from FULL model (simple effects because interaction is supported)
emmeans(score.full, pairwise ~ site | timepoint)
emmeans(score.full, pairwise ~ timepoint | site)

## Plot full model effects
plot(allEffects(score.full), ylab = "outcome")

### formal plots
emm_inter_score <- as.data.frame(emmeans(score.full, ~ site * timepoint))


p_diversity <- ggplot(emm_inter_score, aes(x = timepoint, y = emmean, color = site, group = site)) +
  geom_jitter(data = xdata, 
              aes(x = timepoint, y = comforting, color = site), 
              alpha = 0.2,       
              width = 0.2,       
              height = 0.05,     
              size = 1.5) +      

  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                width = 0.15, 
                size = 1, 
                position = position_dodge(0.2)) + 
  
  geom_line(size = 1.2, position = position_dodge(0.2)) +
  geom_point(size = 3.5, position = position_dodge(0.2)) +

  scale_color_npg() + 
  theme_classic() +  
  labs(x = "", 
       y = "Comforting Behaviour Diversity Score", 
       color = "Site") +
  theme(text = element_text(size = 13),
        legend.position = "bottom")+
  scale_y_continuous(breaks = c( 0,1,2,3))+
  coord_capped_cart(bottom='none', left='both', ylim = c(0,3))


print(p_diversity)

############################################################
## Model 6: Longitudinal models.
##
## Strategy:
## 1. Convert long-format data to wide format
## 2. Run separate GLMs for each timepoint pair:
##    - 16 -> 22
##    - 22 -> 36
##    - 16 -> 36
##
## DVs:
## - Occurrence (comforts): binomial glm
## - Score (comforting): Gaussian glm
############################################################

cat("\n==================== Model 6: Longitudinal models ====================\n")

## Packages needed for reshaping
library(dplyr)
library(tidyr)

############################################################
## 6.1 Convert long data to wide data
############################################################
## One row per participant, with separate columns for each timepoint

wide <- xdata %>%
  select(participant, site, sex, timepoint, comforts, comforting) %>%
  pivot_wider(
    id_cols = c(participant, site, sex),
    names_from = timepoint,
    values_from = c(comforts, comforting)
  )

## Check automatically created names
names(wide)

## Rename columns to simpler labels
names(wide) <- c(
  "participant",
  "site",
  "sex",
  "comforts16",
  "comforts22",
  "comforts36",
  "comforting16",
  "comforting22",
  "comforting36"
)

## Quick check
str(wide)

############################################################
## 6.2 Recode wide-format binary variables
############################################################

wide$comforts16 <- factor(wide$comforts16, levels = c("N", "Y"))
wide$comforts22 <- factor(wide$comforts22, levels = c("N", "Y"))
wide$comforts36 <- factor(wide$comforts36, levels = c("N", "Y"))

## Optional checks
table(wide$comforts16, useNA = "ifany")
table(wide$comforts22, useNA = "ifany")
table(wide$comforts36, useNA = "ifany")

summary(wide$comforting16)
summary(wide$comforting22)
summary(wide$comforting36)

############################################################
## 6.3 Sample size checks for each longitudinal pair
############################################################

cat("\n==================== Longitudinal sample sizes ====================\n")

cat("Occurrence 16 -> 22 N = ",
    sum(complete.cases(wide[, c("comforts16", "comforts22", "site", "sex")])),
    "\n")

cat("Occurrence 22 -> 36 N = ",
    sum(complete.cases(wide[, c("comforts22", "comforts36", "site", "sex")])),
    "\n")

cat("Occurrence 16 -> 36 N = ",
    sum(complete.cases(wide[, c("comforts16", "comforts36", "site", "sex")])),
    "\n")

cat("Score 16 -> 22 N = ",
    sum(complete.cases(wide[, c("comforting16", "comforting22", "site", "sex")])),
    "\n")

cat("Score 22 -> 36 N = ",
    sum(complete.cases(wide[, c("comforting22", "comforting36", "site", "sex")])),
    "\n")

cat("Score 16 -> 36 N = ",
    sum(complete.cases(wide[, c("comforting16", "comforting36", "site", "sex")])),
    "\n")

############################################################
## 6.4 Longitudinal occurrence models (binary)
############################################################
## Outcome = later comforting occurrence
## Predictor = earlier comforting occurrence
## Interaction = whether this predictive relation varies by site

cat("\n==================== Model 6A: Longitudinal occurrence ====================\n")
## 16 -> 22
dat_16_22_occ <- wide[complete.cases(wide[, c("comforts16", "comforts22", "site", "sex")]), ]

occ_16_22_null <- glm(
  comforts22 ~ sex,
  data = dat_16_22_occ,
  family = binomial
)

occ_16_22_reduced <- glm(
  comforts22 ~ comforts16 + site + sex,
  data = dat_16_22_occ,
  family = binomial
)

occ_16_22_full <- glm(
  comforts22 ~ comforts16 * site + sex,
  data = dat_16_22_occ,
  family = binomial
)

anova(occ_16_22_null, occ_16_22_full, test = "Chisq")
anova(occ_16_22_reduced, occ_16_22_full, test = "Chisq")

drop1(occ_16_22_reduced, test = "Chisq")
emm_site_occ <- emmeans(occ_16_22_reduced, ~ site)
pairs(emm_site_occ, adjust = "tukey")

## 22 -> 36
dat_22_36_occ <- wide[complete.cases(wide[, c("comforts22", "comforts36", "site", "sex")]), ]

occ_22_36_null <- glm(
  comforts36 ~ sex,
  data = dat_22_36_occ,
  family = binomial
)

occ_22_36_reduced <- glm(
  comforts36 ~ comforts22 + site + sex,
  data = dat_22_36_occ,
  family = binomial
)

occ_22_36_full <- glm(
  comforts36 ~ comforts22 * site + sex,
  data = dat_22_36_occ,
  family = binomial
)
anova(occ_22_36_null, occ_22_36_full, test = "Chisq")
anova(occ_22_36_reduced, occ_22_36_full, test = "Chisq")
drop1(occ_22_36_reduced, test = "Chisq")
emm_site_occ <- emmeans(occ_22_36_reduced, ~ site)
pairs(emm_site_occ, adjust = "tukey")
## 16 -> 36
dat_16_36_occ <- wide[complete.cases(wide[, c("comforts16", "comforts36", "site", "sex")]), ]
occ_16_36_null <- glm(
  comforts36 ~ sex,
  data = dat_16_36_occ,
  family = binomial
)
occ_16_36_reduced <- glm(
  comforts36 ~ comforts16 + site + sex,
  data = dat_16_36_occ,
  family = binomial
)
occ_16_36_full <- glm(
  comforts36 ~ comforts16 * site + sex,
  data = dat_16_36_occ,
  family = binomial
)

anova(occ_16_36_null, occ_16_36_full, test = "Chisq")
anova(occ_16_36_reduced, occ_16_36_full, test = "Chisq")

drop1(occ_16_36_reduced, test = "Chisq")
emm_site_occ <- emmeans(occ_16_36_reduced, ~ site)
pairs(emm_site_occ, adjust = "tukey")

##plot
library(emmeans)
library(ggplot2)

emm_occ_16_22 <- as.data.frame(
  emmeans(occ_16_22_reduced, ~ site, type = "response")
)
emm_occ_16_22$interval <- "16 to 22 months"

emm_occ_22_36 <- as.data.frame(
  emmeans(occ_22_36_reduced, ~ site, type = "response")
)
emm_occ_22_36$interval <- "22 to 36 months"

emm_occ_16_36 <- as.data.frame(
  emmeans(occ_16_36_reduced, ~ site, type = "response")
)
emm_occ_16_36$interval <- "16 to 36 months"

plotdat_occ <- rbind(
  emm_occ_16_22,
  emm_occ_22_36,
  emm_occ_16_36
)
plotdat_occ

library(dplyr)

raw_16_22 <- dat_16_22_occ %>%
  group_by(site) %>%
  summarise(prob = mean(comforts22)) %>%
  mutate(interval = "16 to 22 months")

raw_22_36 <- dat_22_36_occ %>%
  group_by(site) %>%
  summarise(prob = mean(comforts36)) %>%
  mutate(interval = "22 to 36 months")

raw_16_36 <- dat_16_36_occ %>%
  group_by(site) %>%
  summarise(prob = mean(comforts36)) %>%
  mutate(interval = "16 to 36 months")

raw_occ <- rbind(raw_16_22, raw_22_36, raw_16_36)

ggplot() +
  
  # raw proportions
  geom_point(
    data = raw_occ,
    aes(x = site, y = prob, colour = interval),
    size = 2,
    alpha = 0.4,
    position = position_dodge(width = 0.3)
  ) +
  
  # model predictions
  geom_point(
    data = plotdat_occ,
    aes(x = site, y = prob, colour = interval, shape = interval),
    size = 3,
    position = position_dodge(width = 0.3)
  ) +
  
  geom_line(
    data = plotdat_occ,
    aes(x = site, y = prob, colour = interval, group = interval),
    position = position_dodge(width = 0.3)
  ) +
  
  geom_errorbar(
    data = plotdat_occ,
    aes(x = site, ymin = asymp.LCL, ymax = asymp.UCL, colour = interval),
    width = 0.15,
    position = position_dodge(width = 0.3)
  ) +
  
  scale_colour_manual(values = c(
    "16 to 22 months" = "#1b9e77",
    "22 to 36 months" = "#d95f02",
    "16 to 36 months" = "#7570b3"
  )) +
  
  ylab("Predicted probability of comforting") +
  xlab("Site") +
  ylim(0,1) +
  theme_classic()
############################################################
## 6.5 Longitudinal score models (Gaussian)
############################################################
## Outcome = later comforting score
## Predictor = earlier comforting score
## Interaction = whether this predictive relation varies by site

cat("\n==================== Model 6B: Longitudinal score ====================\n")

## 16 -> 22
dat_16_22 <- wide[complete.cases(wide[, c("comforting16", "comforting22", "site", "sex")]), ]

score_16_22_null <- glm(
  comforting22 ~ sex,
  data = dat_16_22,
  family = gaussian
)

score_16_22_reduced <- glm(
  comforting22 ~ comforting16 + site + sex,
  data = dat_16_22,
  family = gaussian
)

score_16_22_full <- glm(
  comforting22 ~ comforting16 * site + sex,
  data = dat_16_22,
  family = gaussian
)

anova(score_16_22_null, score_16_22_full, test = "F")
anova(score_16_22_reduced, score_16_22_full, test = "F")

drop1(score_16_22_reduced, test="F")

emmeans(score_16_22_reduced, ~ site)
pairs(emmeans(score_16_22_reduced, ~ site), adjust = "tukey")

## 22 -> 36
dat_22_36 <- wide[complete.cases(wide[, c("comforting22", "comforting36", "site", "sex")]), ]

score_22_36_null <- glm(
  comforting36 ~ sex,
  data = dat_22_36,
  family = gaussian
)

score_22_36_reduced <- glm(
  comforting36 ~ comforting22 + site + sex,
  data = dat_22_36,
  family = gaussian
)

score_22_36_full <- glm(
  comforting36 ~ comforting22 * site + sex,
  data = dat_22_36,
  family = gaussian
)

anova(score_22_36_null, score_22_36_full, test = "F")
anova(score_22_36_reduced, score_22_36_full, test = "F")

drop1(score_22_36_reduced, test = "F")
pairs(emmeans(score_22_36_reduced, ~ site), adjust = "tukey")

## 16 -> 36
dat_16_36 <- wide[complete.cases(wide[, c("comforting16", "comforting36", "site", "sex")]), ]

score_16_36_null <- glm(
  comforting36 ~ sex,
  data = dat_16_36,
  family = gaussian
)

score_16_36_reduced <- glm(
  comforting36 ~ comforting16 + site + sex,
  data = dat_16_36,
  family = gaussian
)

score_16_36_full <- glm(
  comforting36 ~ comforting16 * site + sex,
  data = dat_16_36,
  family = gaussian
)

anova(score_16_36_null, score_16_36_full, test = "F")
anova(score_16_36_reduced, score_16_36_full, test = "F")

drop1(score_16_36_reduced, test = "F")
pairs(emmeans(score_16_36_reduced, ~ site), adjust = "tukey")

##plot 
library(dplyr)
library(ggplot2)

plot_16_22 <- dat_16_22 %>%
  transmute(
    site = site,
    earlier = comforting16,
    later = comforting22,
    interval = "16 to 22 months"
  )

plot_22_36 <- dat_22_36 %>%
  transmute(
    site = site,
    earlier = comforting22,
    later = comforting36,
    interval = "22 to 36 months"
  )

plot_16_36 <- dat_16_36 %>%
  transmute(
    site = site,
    earlier = comforting16,
    later = comforting36,
    interval = "16 to 36 months"
  )

plotdat_score <- bind_rows(plot_16_22, plot_22_36, plot_16_36)

ggplot(plotdat_score, aes(x = earlier, y = later, colour = site)) +
  geom_jitter(width = 0.08, height = 0.08, alpha = 0.35, size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~ interval) +
  labs(
    title = "Longitudinal prediction of comforting scores across sites",
    x = "Earlier comforting score",
    y = "Later comforting score",
    colour = "Site"
  ) +
  theme_classic()

##table
xdata %>% group_by(site, timepoint) %>%
  summarise(
    total_N = n(),
    valid_N = sum(!is.na(comforts)),
    empathy_count = sum (comforts == 'Y', na.rm=T),
    empathy_percentage = mean(comforts == 'Y', na.rm=T)*100,
    phy_empathy = sum(physical == 'Y', na.rm=T),
    phy_percentage = mean(physical == 'Y', na.rm=T)*100,
    verbal_empathy = sum(verbal == 'Y',na.rm = T),
    verbal_percentage = mean(verbal == 'Y',na.rm = T)*100,
    object_empathy = sum(object == 'Y',na.rm=T),
    object_percentage = mean(object == 'Y',na.rm = T)*100,
    conforting_mean = mean(comforting,na.rm=T),
    conforting_sd = sd(comforting, na.rm=T)
  )
Table1 <- xdata %>% group_by(site, timepoint) %>%
  summarise(
    total_N = n(),
    valid_N = sum(!is.na(comforts)),
    empathy_count = sum (comforts == 'Y', na.rm=T),
    empathy_percentage = mean(comforts == 'Y', na.rm=T)*100,
    phy_empathy = sum(physical == 'Y', na.rm=T),
    phy_percentage = mean(physical == 'Y', na.rm=T)*100,
    verbal_empathy = sum(verbal == 'Y',na.rm = T),
    verbal_percentage = mean(verbal == 'Y',na.rm = T)*100,
    object_empathy = sum(object == 'Y',na.rm=T),
    object_percentage = mean(object == 'Y',na.rm = T)*100,
    conforting_mean = mean(comforting,na.rm=T),
    conforting_sd = sd(comforting, na.rm=T)
  )
View(Table1)

############################################################
## Participant table
############################################################

library(dplyr)

dob_var  <- "DOB"
test_var <- "DOT"

participant_table <- xdata %>%
  mutate(
    dob = as.Date(.data[[dob_var]]),
    test_date = as.Date(.data[[test_var]]),
    
dob = as.Date(.data[[dob_var]], format = "%d/%m/%Y"),
test_date = as.Date(.data[[test_var]], format = "%d/%m/%Y"),
    
    age_months = as.numeric(test_date - dob) / 30.4375,
    
    timepoint = factor(
      timepoint,
      levels = c("16 months", "22 months", "36 months")
    )
  ) %>%
  group_by(site, timepoint) %>%
  summarise(
    N = n_distinct(participant),
    
    female_n = sum(sex %in% c("F", "Female"), na.rm = TRUE),
    male_n   = sum(sex %in% c("M", "Male"), na.rm = TRUE),
    
    female_pct = female_n / N * 100,
    
    avg_age = mean(age_months, na.rm = TRUE),
    sd_age = sd(age_months, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  mutate(
    site = ifelse(site == "UK", "Durham (UK)", as.character(site)),
    
    `Female (n, %)` = sprintf("%d (%.1f%%)", female_n, female_pct),
    `Male (n)` = male_n,
    
    `Mean age (months, SD)` = sprintf("%.2f (%.2f)", avg_age, sd_age)
  ) %>%
  select(
    Site = site,
    Timepoint = timepoint,
    N,
    `Female (n, %)`,
    `Male (n)`,
    `Mean age (months, SD)`
  ) %>%
  arrange(Site, Timepoint)

participant_table
View(participant_table)

write.csv(participant_table, "participant_table.csv", row.names = FALSE)
