# BC : Association of Untargeted Metabolomics and Breast Cancer Risk

library(broom)
library(survival)
library(haven)
library(tidyverse)
library(ggrepel)
library(data.table)
library(readr)

library(kableExtra)
library(qwraps2)
options(qwraps2_markup = "markdown")
library(dplyr)
library(MASS)

#___SABC data ___##
SABC_data <- read_sas("sabc_match_confl.sas7bdat")

# scale to unit variance
# dataset contains metabolomics data (positive mode), which were previously scaled and log2-transformed.
scale_meta_pos <- comb_all_ %>% select(1:941)

# Merge metabolomics data with breast cancer risk factor covariates
data_pos_all <- scale_meta_pos %>% left_join(SABC_data, by = "PERSONID")

# Convert metabolomics data columns to numeric format
scale_meta_pos_ <- as.data.frame(sapply(data_pos_all[, 2:941], as.numeric))

#Extract covariates and lifestyle factors
covariate_all_ <- data_pos_all %>% select(1, 942:1004)


#Selection of metabolites names
library(reticulate)
library(janitor)

meta_names <- scale_meta_pos_[1,]
meta_names_ <- t(meta_names)
meta_names_0 <- as.data.frame(meta_names_)
meta_names_1 <- setDT(meta_names_0, keep.rownames=TRUE)[]
meta_names_1 <- meta_names_1 %>% rename (Compound=rn)
meta_names_1 <- meta_names_1 %>% select(Compound)
meta_names_1 <- meta_names_1 %>% mutate(display_name=Compound, description=Compound)

################################INTERACTION TERM#########################

subtype_ER <- covariate_all_ %>% mutate(BC_subtype=as.factor(case_when (
  Status==1 & ER_statusIndex == 1 ~ "ER_Positive",
  Status==1 & ER_statusIndex == 0 ~ "ER_Negative",
  Status==0 ~"Control"
)))

subtype_PR <- covariate_all_ %>% mutate(BC_subtype=as.factor(case_when (
  Status==1 & PR_statusIndex == 1 ~ "PR_Positive",
  Status==1 & PR_statusIndex == 0 ~ "PR_Negative",
  Status==0 ~"Control"
)))

subtype_HER2 <- covariate_all_ %>% mutate(BC_subtype=as.factor(case_when (
  Status==1 & HER2_statusIndex == 1 ~ "HER2_Positive",
  Status==1 & HER2_statusIndex == 0 ~ "HER2_Negative",
  Status==1 & HER2_statusIndex == 888 ~"Unknown",
  Status==0 ~"Control"
)))

ER_pos_tab_0 <- ER_pos_tab_ %>% mutate(BC_subtype=1)
ER_neg_tab_0 <- ER_neg_tab_ %>% mutate(BC_subtype=2)

PR_pos_tab_0 <- PR_pos_tab_ %>% mutate(BC_subtype=1)
PR_neg_tab_0 <- PR_neg_tab_ %>% mutate(BC_subtype=2)

HER2_pos_tab_0 <- HER2_pos_tab_ %>% mutate(BC_subtype=1)
HER2_neg_tab_0 <- HER2_neg_tab_ %>% mutate(BC_subtype=2)

ER_tab <- rbind(ER_pos_tab_0, ER_neg_tab_0)
PR_tab <- rbind(PR_pos_tab_0, PR_neg_tab_0)
HER2_tab <- rbind(HER2_pos_tab_0, HER2_neg_tab_0)

##Definition of interaction function
#ER
interaction_pos_ER <- function(x) {
  mod1 <- clogit(as.numeric(Status) ~ x + ETH + WC + V443_Menopausal + Bmi + OCEver_1 + breastfed + 
                   ever_alcohol + parous_1 + VIH + strata(MATCH), data=ER_tab)
  mod2 <- clogit(as.numeric(Status) ~ x + ETH + WC + V443_Menopausal + Bmi + OCEver_1 + breastfed + 
                   ever_alcohol + parous_1 + VIH + strata(MATCH) + BC_subtype + x*BC_subtype, data=ER_tab)
  interaction_test <- tidy(anova(mod2, mod1)) %>% select(p.value) }

p_interaction_pos_ER <- apply(scale_meta_pos_, 2, interaction_pos_ER) %>% map_df(tidy, exponentiate = T) %>% 
  bind_cols(meta_names_1) %>% select(Compound, metabolite_name, mean) %>% mutate(p.value = mean) %>% mutate(p.adj = p.adjust(p.value, method = "fdr"))

#PR
interaction_pos_PR <- function(x) {
  mod1 <- clogit(as.numeric(Status) ~ x + ETH + WC + V443_Menopausal + Bmi + OCEver_1 + breastfed + 
                   ever_alcohol + parous_1 + VIH + strata(MATCH), data=PR_tab)
  mod2 <- clogit(as.numeric(Status) ~ x + ETH + WC + V443_Menopausal + Bmi + OCEver_1 + breastfed + 
                   ever_alcohol + parous_1 + VIH + strata(MATCH) + BC_subtype + x*BC_subtype, data=PR_tab)
  interaction_test <- tidy(anova(mod2, mod1)) %>% select(p.value) }

p_interaction_pos_PR <- apply(scale_meta_pos_, 2, interaction_pos_PR) %>% map_df(tidy, exponentiate = T) %>% 
  bind_cols(meta_names_1) %>% select(Compound, metabolite_name, mean) %>% mutate(p.value = mean) %>% mutate(p.adj = p.adjust(p.value, method = "fdr"))


#HER2
HER2_all <- rbind(HER2_pos_tab_, HER2_neg_tab_)
dim(HER2_all)

scale_meta_pos_HER2_0  <- HER2_all %>% select(65:1004)

interaction_pos_HER2 <- function(x) {
  mod1 <- clogit(as.numeric(Status) ~ x + ETH + WC + V443_Menopausal + Bmi + OCEver_1 + breastfed + 
                   ever_alcohol + parous_1 + VIH + strata(MATCH), data=HER2_tab)
  mod2 <- clogit(as.numeric(Status) ~ x + ETH + WC + V443_Menopausal + Bmi + OCEver_1 + breastfed + 
                   ever_alcohol + parous_1 + VIH + strata(MATCH) + BC_subtype + x*BC_subtype, data=HER2_tab)
  interaction_test <- tidy(anova(mod2, mod1)) %>% select(p.value) }

p_interaction_pos_HER2 <- apply(scale_meta_pos_HER2_0, 2, interaction_pos_HER2) %>% map_df(tidy, exponentiate = T) %>% 
  bind_cols(meta_names_1) %>% select(Compound, metabolite_name, mean) %>% mutate(p.value = mean) %>% mutate(p.adj = p.adjust(p.value, method = "fdr"))

write_csv(p_interaction_pos_ER, path = "XXXX/Tables - revised/P_interaction_ER.csv")
write_csv(p_interaction_pos_PR, path = "XXXX/Tables - revised/P_interaction_PR.csv")
write_csv(p_interaction_pos_HER2, path = "XXXX/Tables - revised/P_interaction_HER2.csv")

# Subroup analyses by VIH positivity
table(covariate_all_$VIH)

HIV_yes <- covariate_all_$VIH == 1
HIV_no <- covariate_all_$VIH == 0

HIV_yes_tab <- covariate_all_[HIV_yes, ]
HIV_no_tab <- covariate_all_[HIV_no, ]

# VIH YES
fits_vih_y <- apply(scale_meta_pos_[HIV_yes, ], 2, function(x) clogit(as.numeric(Status) ~ x + ETH + WC + V443_Menopausal + Bmi + 
                                                                        ever_alcohol + parous_1 + strata(MATCH), data=covariate_all_[HIV_yes, ]))
# VIH NOS
fits_vih_n <- apply(scale_meta_pos_[HIV_no, ], 2, function(x) clogit(as.numeric(Status) ~ x + ETH + WC + V443_Menopausal + Bmi + breastfed + 
                                                                       ever_alcohol + parous_1 + strata(MATCH), data=covariate_all_[HIV_no, ]))

all_vih_y <- tidy.output(fits_vih_y)
all_vih_n <- tidy.output(fits_vih_n)

# WITH HIV Retain only metabolite groups with at least P.value < 0.05
tab_vih_yes <- bind_rows("HIV YES" = all_vih_y, .id = "Group") %>%
  arrange(Group, -OR) %>% 
  group_by(Compound) %>%
  select(Compound, description, everything()) 
tab_vih_yes

# Retain only metabolite groups with at least one FDR < 0.05
tab_vih_yes_ <- bind_rows("HIV YES" = all_vih_y, .id = "Group") %>%
  arrange(Group, -OR) %>% 
  group_by(Compound) %>%
  filter(min(FDR) < 0.05) %>%
  select(Compound, description, everything()) 
tab_vih_yes_

 
# WITHOUT HIV Retain only metabolite groups with at least P.value < 0.05
tab_vih_no <- bind_rows("HIV NO" = all_vih_n, .id = "Group") %>%
  arrange(Group, -OR) %>% 
  group_by(Compound) %>%
  select(Compound, description, everything()) 
tab_vih_no

# Retain only metabolite groups with at least one FDR < 0.05
tab_vih_no_ <- bind_rows("HIV NO" = all_vih_n, .id = "Group") %>%
  arrange(Group, -OR) %>% 
  group_by(Compound) %>%
  filter(min(FDR) < 0.05) %>%
  select(Compound, description, everything()) 
tab_vih_no_

write_csv(tab_vih_yes, path = "XXXX/HIV_YES.csv")
write_csv(tab_vih_no, path = "XXXX/HIV_NO.csv")

# DIABETES
table(covariate_all_$DMT)

DMT_yes <- covariate_all_$DMT == 1
DMT_no <- covariate_all_$DMT == 0

DMT_yes_tab <- covariate_all_[DMT_yes, ]
DMT_no_tab <- covariate_all_[DMT_no, ]


# DIABETES YES
fits_DMT_y <- apply(scale_meta_pos_[DMT_yes, ], 2, function(x) clogit(as.numeric(Status) ~ x + ETH + WC + V443_Menopausal + Bmi + OCEver_1 + breastfed + 
                                                                        ever_alcohol + parous_1 + VIH + strata(MATCH), data=covariate_all_[DMT_yes, ]))
# DIABETES NO
fits_DMT_n <- apply(scale_meta_pos_[DMT_no, ], 2, function(x) clogit(as.numeric(Status) ~ x + ETH + WC + V443_Menopausal + Bmi + OCEver_1 + breastfed + 
                                                                        ever_alcohol + parous_1 + VIH + strata(MATCH), data=covariate_all_[DMT_no, ]))

all_DMT_y <- tidy.output(fits_DMT_y)
all_DMT_n <- tidy.output(fits_DMT_n)

# WITH DIABETES Retain only metabolite groups with at least P.value < 0.05
tab_DMT_yes <- bind_rows("DMT YES" = all_DMT_y, .id = "Group") %>%
  arrange(Group, -OR) %>% 
  group_by(Compound) %>%
  select(Compound, description, everything()) 
tab_DMT_yes

# Retain only metabolite groups with at least one FDR < 0.05
tab_DMT_yes_ <- bind_rows("DMT YES" = all_DMT_y, .id = "Group") %>%
  arrange(Group, -OR) %>% 
  group_by(Compound) %>%
  filter(min(FDR) < 0.05) %>%
  select(Compound, description, everything()) 
tab_DMT_yes_

# WITHOUT DIABETES Retain only metabolite groups with at least P.value < 0.05
tab_DMT_no <- bind_rows("DMT NO" = all_DMT_n, .id = "Group") %>%
  arrange(Group, -OR) %>% 
  group_by(Compound) %>%
  select(Compound, description, everything()) 
tab_DMT_no

# Retain only metabolite groups with at least one FDR < 0.05
tab_DMT_no_ <- bind_rows("DMT NO" = all_DMT_n, .id = "Group") %>%
  arrange(Group, -OR) %>% 
  group_by(Compound) %>%
  filter(min(FDR) < 0.05) %>%
  select(Compound, description, everything()) 
tab_DMT_no_


write_csv(tab_DMT_yes, path = "XXXX/DMT_YES.csv")
write_csv(tab_DMT_no, path = "XXXX/DMT_NO.csv")


# ORAL CONTRACEPTIVE
table(covariate_all_$OCEver_1)

OC_ne <- covariate_all_$OCEver_1 == 1
OC_ev <- covariate_all_$OCEver_1 == 2

table(OC_ne)
table(OC_ev)

# ORAL CONTRACEPTIVE NEVER
fits_OC_n <- apply(scale_meta_pos_[OC_ne, ], 2, function(x) clogit(as.numeric(Status) ~ x + ETH + WC + V443_Menopausal + Bmi + breastfed + 
                                                                     ever_alcohol + parous_1 + VIH + strata(MATCH), data=covariate_all_[OC_ne, ]))
# ORAL CONTRACEPTIVE EVER
fits_OC_ev <- apply(scale_meta_pos_[OC_ev, ], 2, function(x) clogit(as.numeric(Status) ~ x + ETH + WC + V443_Menopausal + Bmi + breastfed + 
                                                                      ever_alcohol + parous_1 + VIH + strata(MATCH), data=covariate_all_[OC_ev, ]))
all_OC_ne <- tidy.output(fits_OC_n)
all_OC_ev <- tidy.output(fits_OC_ev)

# NEVER ORAL CONTRACEPTIVE USE Retain only metabolite groups with at least P.value < 0.05
tab_OC_ne <- bind_rows("OC NEVER" = all_OC_ne, .id = "Group") %>%
  arrange(Group, -OR) %>% 
  group_by(Compound) %>%
  select(Compound, description, everything()) 
tab_OC_ne

# Retain only metabolite groups with at least one FDR < 0.05
tab_OC_ne_ <- bind_rows("OC NEVER" = all_OC_ne, .id = "Group") %>%
  arrange(Group, -OR) %>% 
  group_by(Compound) %>%
  filter(min(FDR) < 0.05) %>%
  select(Compound, description, everything()) 
tab_OC_ne_


#EVER ORAL CONTRACEPTIVE USE Retain only metabolite groups with at least P.value < 0.05
tab_OC_ev <- bind_rows("OC EVER" = all_OC_ev, .id = "Group") %>%
  arrange(Group, -OR) %>% 
  group_by(Compound) %>%
  select(Compound, description, everything()) 
tab_OC_ev

# Retain only metabolite groups with at least one FDR < 0.05
tab_OC_ev_ <- bind_rows("OC EVER" = all_OC_ev, .id = "Group") %>%
  arrange(Group, -OR) %>% 
  group_by(Compound) %>%
  filter(min(FDR) < 0.05) %>%
  select(Compound, description, everything()) 
tab_OC_ev_


write_csv(tab_OC_ne, path = "XXXX/OC_NEVER.csv")
write_csv(tab_OC_ev, path = "XXXX/OC_EVER.csv")

#Menopausal status
interaction_pos_meno <- function(x) {
  mod1 <- clogit(as.numeric(Status) ~ x + ETH + WC + Bmi + OCEver_1 + breastfed + 
                   ever_alcohol + parous_1 + VIH + strata(MATCH), data=covariate_all_)
  mod2 <- clogit(as.numeric(Status) ~ x + ETH + WC + Bmi + OCEver_1 + breastfed + 
                   ever_alcohol + parous_1 + VIH + strata(MATCH) + V443_Menopausal + x*V443_Menopausal, data=covariate_all_)
  interaction_test <- tidy(anova(mod2, mod1)) %>% select(p.value) }

p_interaction_pos_meno <- apply(scale_meta_pos_, 2, interaction_pos_meno) %>% map_df(tidy, exponentiate = T) %>% 
  bind_cols(meta_names_1) %>% select(Compound, metabolite_name, mean) %>% mutate(p.value = mean) %>% mutate(p.adj = p.adjust(p.value, method = "fdr"))

#HIV
interaction_pos_hiv <- function(x) {
  mod1 <- clogit(as.numeric(Status) ~ x + ETH + WC + V443_Menopausal + Bmi + OCEver_1 + breastfed + 
                   ever_alcohol + parous_1 + strata(MATCH), data=covariate_all_)
  mod2 <- clogit(as.numeric(Status) ~ x + ETH + WC + V443_Menopausal + Bmi + OCEver_1 + breastfed + 
                   ever_alcohol + parous_1 + strata(MATCH) + VIH + x*VIH, data=covariate_all_)
  interaction_test <- tidy(anova(mod2, mod1)) %>% select(p.value) }

p_interaction_pos_hiv <- apply(scale_meta_pos_, 2, interaction_pos_hiv) %>% map_df(tidy, exponentiate = T) %>% 
  bind_cols(meta_names_1) %>% select(Compound, metabolite_name, mean) %>% mutate(p.value = mean) %>% mutate(p.adj = p.adjust(p.value, method = "fdr"))

#DMT
interaction_pos_DMT <- function(x) {
  mod1 <- clogit(as.numeric(Status) ~ x + ETH + WC + V443_Menopausal + Bmi + OCEver_1 + breastfed + 
                   ever_alcohol + parous_1 + VIH + strata(MATCH), data=covariate_all_)
  mod2 <- clogit(as.numeric(Status) ~ x + ETH + WC + V443_Menopausal + Bmi + OCEver_1 + breastfed + 
                   ever_alcohol + parous_1 + VIH + strata(MATCH) + DMT + x*DMT, data=covariate_all_)
  interaction_test <- tidy(anova(mod2, mod1)) %>% select(p.value) }

p_interaction_pos_DMT <- apply(scale_meta_pos_, 2, interaction_pos_DMT) %>% map_df(tidy, exponentiate = T) %>% 
  bind_cols(meta_names_1) %>% select(Compound, metabolite_name, mean) %>% mutate(p.value = mean) %>% mutate(p.adj = p.adjust(p.value, method = "fdr"))

#OC
interaction_pos_OC <- function(x) {
  mod1 <- clogit(as.numeric(Status) ~ x + ETH + WC + V443_Menopausal + Bmi + breastfed + 
                   ever_alcohol + parous_1 + VIH + strata(MATCH), data=covariate_all_)
  mod2 <- clogit(as.numeric(Status) ~ x + ETH + WC + V443_Menopausal + Bmi + breastfed + 
                   ever_alcohol + parous_1 + VIH + strata(MATCH) + OCEver_1 + x*OCEver_1, data=covariate_all_)
  interaction_test <- tidy(anova(mod2, mod1)) %>% select(p.value) }

p_interaction_pos_OC <- apply(scale_meta_pos_, 2, interaction_pos_OC) %>% map_df(tidy, exponentiate = T) %>% 
  bind_cols(meta_names_1) %>% select(Compound, metabolite_name, mean) %>% mutate(p.value = mean) %>% mutate(p.adj = p.adjust(p.value, method = "fdr"))


write_csv(p_interaction_pos_meno, path = "XXXX/P_interaction_meno.csv")
write_csv(p_interaction_pos_hiv, path = "XXXX/P_interaction_HIV.csv")
write_csv(p_interaction_pos_DMT, path = "XXXX/P_interaction_DMT.csv")
write_csv(p_interaction_pos_OC, path = "XXXX/P_interaction_OC.csv")

#WAIST CIRCUNFERENCE
mean(covariate_all_$WC, na.rm = TRUE)
mean(covariate_all_$WC, na.rm = TRUE)

#WC
covariate_all_$WC_cat[covariate_all_$WC<=80]=1
covariate_all_$WC_cat[covariate_all_$WC>80]=2

table(covariate_all_$WC_cat)


WC_80_l <- covariate_all_$WC_cat == 1
WC_80_m <- covariate_all_$WC_cat == 2

table(WC_80_l)
table(WC_80_m)

# WC < 80 cm
fits_WC_l <- apply(scale_meta_pos_[WC_80_l, ], 2, function(x) clogit(as.numeric(Status) ~ x + ETH + V443_Menopausal + Bmi + OCEver_1 + breastfed + 
                                                                      ever_alcohol + parous_1 + strata(MATCH), data=covariate_all_[WC_80_l, ]))
# WC > 80 cm
fits_WC_m <- apply(scale_meta_pos_[WC_80_m, ], 2, function(x) clogit(as.numeric(Status) ~ x + ETH + V443_Menopausal + Bmi + OCEver_1 + breastfed + 
                                                                       ever_alcohol + parous_1 + VIH + strata(MATCH), data=covariate_all_[WC_80_m, ]))

all_WC_l <- tidy.output(fits_WC_l)
all_WC_m <- tidy.output(fits_WC_m)


# WC < 80 cm retain only metabolite groups with at least P.value < 0.05
tab_WC_l <- bind_rows("WC <= 80 cm" = all_WC_l, .id = "Group") %>%
  arrange(Group, -OR) %>% 
  group_by(Compound) %>%
  select(Compound, description, everything()) 
tab_WC_l

# Retain only metabolite groups with at least one FDR < 0.05
tab_WC_l_ <- bind_rows("WC <= 80 cm" = all_WC_l, .id = "Group") %>%
  arrange(Group, -OR) %>% 
  group_by(Compound) %>%
  filter(min(FDR) < 0.05) %>%
  select(Compound, description, everything()) 
tab_WC_l_


# WC > 80 cm retain only metabolite groups with at least P.value < 0.05
tab_WC_m <- bind_rows("WC > 80 cm" = all_WC_m, .id = "Group") %>%
  arrange(Group, -OR) %>% 
  group_by(Compound) %>%
  select(Compound, description, everything()) 
tab_WC_m

# Retain only metabolite groups with at least one FDR < 0.05
tab_WC_m_ <- bind_rows("WC > 80 cm" = all_WC_m, .id = "Group") %>%
  arrange(Group, -OR) %>% 
  group_by(Compound) %>%
  filter(min(FDR) < 0.05) %>%
  select(Compound, description, everything()) 
tab_WC_m_


write_csv(tab_WC_l, path = "XXXX/WC_80_less.csv")
write_csv(tab_WC_m, path = "XXXX/WC_80_more.csv")

#BMI
mean(covariate_all_$Bmi, na.rm = TRUE)
mean(covariate_all_$Bmi, na.rm = TRUE)

#BMI
covariate_all_$Bmi_c[covariate_all_$Bmi<=30]=1
covariate_all_$Bmi_c[covariate_all_$Bmi>30]=2

table(covariate_all_$Bmi_c)

BMI_25_l <- covariate_all_$Bmi_c == 1
BMI_25_m <- covariate_all_$Bmi_c == 2

table(BMI_25_l)
table(BMI_25_m)

# WC < 80 cm
fits_bmi_l <- apply(scale_meta_pos_[BMI_25_l, ], 2, function(x) clogit(as.numeric(Status) ~ x + ETH + WC + V443_Menopausal + OCEver_1 + breastfed + 
                                                                         ever_alcohol + parous_1 + VIH + strata(MATCH), data=covariate_all_[BMI_25_l, ]))
# WC > 80 cm
fits_bmi_m <- apply(scale_meta_pos_[BMI_25_m, ], 2, function(x) clogit(as.numeric(Status) ~ x + ETH + WC + V443_Menopausal + OCEver_1 + breastfed + 
                                                                         ever_alcohol + parous_1 + VIH + strata(MATCH), data=covariate_all_[BMI_25_m, ]))

all_bmi_l <- tidy.output(fits_bmi_l)
all_bmi_m <- tidy.output(fits_bmi_m)


#  retain only metabolite groups with at least P.value < 0.05
tab_bmi_l <- bind_rows("BMI <= 30 " = all_bmi_l, .id = "Group") %>%
  arrange(Group, -OR) %>% 
  group_by(Compound) %>%
  select(Compound, description, everything()) 
tab_bmi_l

# Retain only metabolite groups with at least one FDR < 0.05
tab_bmi_l_ <- bind_rows("BMI <= 30 " = all_bmi_l, .id = "Group") %>%
  arrange(Group, -OR) %>% 
  group_by(Compound) %>%
  filter(min(FDR) < 0.05) %>%
  select(Compound, description, everything()) 
tab_bmi_l_


#  retain only metabolite groups with at least P.value < 0.05
tab_bmi_m <- bind_rows("BMI > 30 " = all_bmi_m, .id = "Group") %>%
  arrange(Group, -OR) %>% 
  group_by(Compound) %>%
  select(Compound, description, everything()) 
tab_bmi_m

# Retain only metabolite groups with at least one FDR < 0.05
tab_bmi_m_ <- bind_rows("BMI > 30 " = all_bmi_m, .id = "Group") %>%
  arrange(Group, -OR) %>% 
  group_by(Compound) %>%
  filter(min(FDR) < 0.05) %>%
  select(Compound, description, everything()) 
tab_bmi_m_

write_csv(tab_bmi_l, path = "XXXX/tab_bmi_l.csv")
write_csv(tab_bmi_m, path = "XXXX/tab_bmi_m.csv")


#WC
interaction_pos_WC <- function(x) {
  mod1 <- clogit(as.numeric(Status) ~ x + ETH + V443_Menopausal + Bmi + OCEver_1 + breastfed + 
                   ever_alcohol + parous_1 + VIH + strata(MATCH), data=covariate_all_)
  mod2 <- clogit(as.numeric(Status) ~ x + ETH + V443_Menopausal + Bmi + OCEver_1 + breastfed + 
                   ever_alcohol + parous_1 + VIH + strata(MATCH) + WC_cat + x*WC_cat, data=covariate_all_)
  interaction_test <- tidy(anova(mod2, mod1)) %>% select(p.value) }

p_interaction_pos_WC <- apply(scale_meta_pos_, 2, interaction_pos_WC) %>% map_df(tidy, exponentiate = T) %>% 
  bind_cols(meta_names_1) %>% select(Compound, metabolite_name, mean) %>% mutate(p.value = mean) %>% mutate(p.adj = p.adjust(p.value, method = "fdr"))

#BMI
interaction_pos_bmi <- function(x) {
  mod1 <- clogit(as.numeric(Status) ~ x + ETH + WC + V443_Menopausal + OCEver_1 + breastfed + 
                   ever_alcohol + parous_1 + VIH + strata(MATCH), data=covariate_all_)
  mod2 <- clogit(as.numeric(Status) ~ x + ETH + WC + V443_Menopausal + OCEver_1 + breastfed + 
                   ever_alcohol + parous_1 + VIH + strata(MATCH) + Bmi_c + x*Bmi_c, data=covariate_all_)
  interaction_test <- tidy(anova(mod2, mod1)) %>% select(p.value) }

p_interaction_pos_bmi <- apply(scale_meta_pos_, 2, interaction_pos_bmi) %>% map_df(tidy, exponentiate = T) %>% 
  bind_cols(meta_names_1) %>% select(Compound, metabolite_name, mean) %>% mutate(p.value = mean) %>% mutate(p.adj = p.adjust(p.value, method = "fdr"))


write_csv(p_interaction_pos_WC, path = "XXXX/P_interaction_WC.csv")
write_csv(p_interaction_pos_bmi, path = "XXXX/P_interaction_bmi.csv")



#Definition of other subtypes of breast cancer

table(covariate_all_$ER_statusIndex)
table(covariate_all_$PR_statusIndex)
table(covariate_all_$HER2_statusIndex)

baseline_all_ <- covariate_all_ %>%filter(Status==1) %>% group_by(ER_statusIndex, PR_statusIndex, HER2_statusIndex) %>% count()

baseline_all_1 <- covariate_all_ %>% mutate(BC_subtype=as.factor(case_when(
  Status==1 & ER_statusIndex==0 & PR_statusIndex==0 & HER2_statusIndex==0 ~ "ER-PR-HER2-",
  Status==1 & ER_statusIndex==1 & PR_statusIndex==1 & HER2_statusIndex==0 ~ "HER2-",
  Status==1 & ER_statusIndex==0 & PR_statusIndex==1 & HER2_statusIndex==0 ~ "HER2-",
  Status==1 & ER_statusIndex==1 & PR_statusIndex==0 & HER2_statusIndex==0 ~ "HER2-",
  Status==1 & ER_statusIndex==1 & PR_statusIndex==1 & HER2_statusIndex==1 ~ "HER2+",
  Status==1 & ER_statusIndex==0 & PR_statusIndex==1 & HER2_statusIndex==1 ~ "HER2+",
  Status==1 & ER_statusIndex==1 & PR_statusIndex==0 & HER2_statusIndex==1 ~ "HER2+",
  Status==1 & ER_statusIndex==0 & PR_statusIndex==0 & HER2_statusIndex==1 ~ "HER2+",
  Status==1 ~ "Other",
  Status==0 ~ "Control"
)))

table(baseline_all_1$BC_subtype)
table(covariate_all_$ER_statusIndex, covariate_all_$PR_statusIndex, baseline_all_1$HER2_statusIndex)

CasesetA <- baseline_all_1 %>% filter(BC_subtype=="HER2-") %>% dplyr::select(MATCH) %>% pull()
CasesetB <- baseline_all_1 %>% filter(BC_subtype=="ER-PR-HER2-") %>% dplyr::select(MATCH) %>% pull()
CasesetC <- baseline_all_1 %>% filter(BC_subtype=="HER2+") %>% dplyr::select(MATCH) %>% pull()

#ER+PR+HER2-
#ER+PR-HER2-
#ER-PR+HER2-
#HER2+
#TRIPLE NEGATIVE

#Create the appropriate variable for grouping cases with their controls, and exclude those with unknown combination
D01_BC_subtype <- baseline_all_1 %>% mutate(BC_subtype=case_when(BC_subtype=="Control" & MATCH %in% CasesetA ~ "HER2-",
                                                                 BC_subtype=="Control" & MATCH %in% CasesetB ~ "ER-PR-HER2-",
                                                                 BC_subtype=="Control" & MATCH %in% CasesetC ~ "HER2+",
                                                                 TRUE ~ as.character(BC_subtype)))%>% filter(!BC_subtype=="Other")
table(D01_BC_subtype$BC_subtype)

# Check 
D01_ <- D01_BC_subtype %>% count(BC_subtype)

baseline_all <- D01_BC_subtype

save(baseline_all, file = "XXXX/baseline_all.RData")

# Oestrogen amd Progesterone receptor positive

baseline_all$subtype[baseline_all$BC_subtype=="HER2-"]=1
baseline_all$subtype[baseline_all$BC_subtype=="ER-PR-HER2-"]=2
baseline_all$subtype[baseline_all$BC_subtype=="HER2+"]=3

baseline_all <- baseline_all %>% drop_na(subtype)

table(baseline_all$subtype)
table(baseline_all$Status, baseline_all$subtype)


ERPR_pos <- (baseline_all$ER_statusIndex == 1 & baseline_all$PR_statusIndex == 1) | baseline_all$Status == 0
ERPR_neg <- (baseline_all$ER_statusIndex == 0 & baseline_all$PR_statusIndex == 0) | baseline_all$Status == 0

ERPR_pos_tab <- baseline_all[ERPR_pos, ]
ERPR_neg_tab <- baseline_all[ERPR_neg, ]

ERPR_pos_tab <- ERPR_pos_tab[duplicated(ERPR_pos_tab$MATCH) | duplicated(ERPR_pos_tab$MATCH, fromLast = T), ]
ERPR_neg_tab <- ERPR_neg_tab[duplicated(ERPR_neg_tab$MATCH) | duplicated(ERPR_neg_tab$MATCH, fromLast = T), ]

table(ERPR_pos_tab$Status)
table(ERPR_neg_tab$Status)


# Tester l'interaction 
library(car)
library(tidyverse)
library(broom)
library(emmeans)
library(skimr)
library(Hmisc)

ERPR_ <- (baseline_all$ER_statusIndex == 1 & baseline_all$PR_statusIndex == 1) | (baseline_all$ER_statusIndex == 0 & baseline_all$PR_statusIndex == 0) | baseline_all$Status == 0
ERPR_tab_ <- baseline_all[ERPR_, ]

ERPR_tab_ <- ERPR_tab_[duplicated(ERPR_tab_$MATCH) | duplicated(ERPR_tab_$MATCH, fromLast = T), ]

table(ERPR_tab_$Status)
table(ERPR_tab_$Status)

ERPR_tab_all <- ERPR_tab_ %>% left_join(scale_meta_pos, by = "PERSONID")
dim(ERPR_tab_all)
scale_meta_pos_s_ <- as.data.frame(sapply(ERPR_tab_all[, 67:1006], as.numeric))
dim(scale_meta_pos_s_)
class(scale_meta_pos_s_)


baseline_all_1 <- baseline_all %>% left_join(scale_meta_pos, by = "PERSONID")
dim(baseline_all_1)
scale_meta_pos_s1 <- as.data.frame(sapply(baseline_all_1[, 67:1006], as.numeric))
dim(scale_meta_pos_s1)
class(scale_meta_pos_s1)

baseline_all_s <- baseline_all_1 %>% select(1:66)


interaction_pos_subtype <- function(x) {
  
  mod1 <- clogit(as.numeric(Status) ~ x + ETH + WC + V443_Menopausal + Bmi + OCEver_1 + breastfed + 
                   ever_alcohol + parous_1 + VIH + strata(MATCH), data=baseline_all_1)
  
  mod2 <- clogit(as.numeric(Status) ~ x + ETH + WC + V443_Menopausal + Bmi + OCEver_1 + breastfed + 
                   ever_alcohol + parous_1 + VIH + strata(MATCH) + subtype + x*subtype, data=baseline_all_1)
  interaction_test <- tidy(anova(mod2, mod1)) %>% select(p.value)}

p_interaction_pos_subtype <- apply(scale_meta_pos_s1, 2, interaction_pos_subtype) %>% map_df(tidy, exponentiate = T) %>% 
  bind_cols(meta_names_1) %>% select(Compound, metabolite_name, mean) %>% mutate(p.value = mean) %>% mutate(p.adj = p.adjust(p.value, method = "fdr"))


dim(scale_meta_pos_s1)
dim(baseline_all_1)

table(baseline_all_1$Status, baseline_all_1$subtype)

write_csv(p_interaction_pos_subtype, path = "XXXX/pintdat_pvalue_pos.csv")

HER2_n  <- baseline_all$BC_subtype=="HER2-"
TNBC  <- baseline_all$BC_subtype=="ER-PR-HER2-"
HER2_p  <- baseline_all$BC_subtype=="HER2+"


#ER POSITIVE AND NEGATIVE #

fits_HER2_n <- apply(scale_meta_pos_s1[HER2_n, ], 2, function(x) clogit(as.numeric(Status) ~ x + ETH + WC + V443_Menopausal + Bmi + OCEver_1 + breastfed + 
                                                                               ever_alcohol + parous_1 + VIH + strata(MATCH), data=baseline_all_1[HER2_n, ]))


fits_TNBC <- apply(scale_meta_pos_s1[TNBC, ], 2, function(x) clogit(as.numeric(Status) ~ x + ETH + WC + V443_Menopausal + Bmi + OCEver_1 + breastfed + 
                                                                                ever_alcohol + parous_1 + VIH + strata(MATCH), data=baseline_all_1[TNBC, ]))

fits_HER2_p <- apply(scale_meta_pos_s1[HER2_p, ], 2, function(x) clogit(as.numeric(Status) ~ x + ETH + WC + V443_Menopausal + Bmi + OCEver_1 + breastfed + 
                                                                      ever_alcohol + parous_1 + VIH + strata(MATCH), data=baseline_all_1[HER2_p, ]))

all_TNBC  <- tidy.output(fits_TNBC)
all_HER2_p  <- tidy.output(fits_HER2_p)
all_HER2_n  <- tidy.output(fits_HER2_n)

# ADJUSTED FOR ETHNICITY, WC, AND MENOPAUSAL STATUS
# Retain only metabolite groups with at least P.value < 0.05
#HER2-
table_HER2_n <- bind_rows("HER2-" = all_HER2_n, .id = "Group") %>%
  arrange(Group, -OR) %>% 
  group_by(Compound) %>%
  select(Compound, description, everything()) 
table_HER2_n


#ER-PR-HER2-
table_TNBC <- bind_rows("ER-PR-HER2-" = all_TNBC, .id = "Group") %>%
  arrange(Group, -OR) %>% 
  group_by(Compound) %>%
  select(Compound, description, everything()) 
table_TNBC

#HER2+
table_HER2_p <- bind_rows("HER2+" = all_HER2_p, .id = "Group") %>%
  arrange(Group, -OR) %>% 
  group_by(Compound) %>%
  select(Compound, description, everything()) 
table_HER2_p

write_csv(table_TNBC, path = "XXXX/table_TNBC.csv")
write_csv(table_HER2_n, path = "XXXX/table_HER2_n.csv")
write_csv(table_HER2_p, path = "XXXX/table_HER2_p.csv")
