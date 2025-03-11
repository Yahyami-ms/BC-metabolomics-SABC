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

#For subsetting
pre_ <- covariate_all_$V443_Menopausal == 2
post_ <- covariate_all_$V443_Menopausal == 1

# Oestrogen receptor positive
ER_pos <- covariate_all_$ER_statusIndex == 1 | covariate_all_$Status == 0
ER_neg <- covariate_all_$ER_statusIndex == 0 | covariate_all_$Status == 0

ER_pos_tab <- covariate_all_[ER_pos, ]
ER_neg_tab <- covariate_all_[ER_neg, ]

ER_pos_tab <- ER_pos_tab[duplicated(ER_pos_tab$MATCH) | duplicated(ER_pos_tab$MATCH, fromLast = T), ]
ER_neg_tab <- ER_neg_tab[duplicated(ER_neg_tab$MATCH) | duplicated(ER_neg_tab$MATCH, fromLast = T), ]
ER_pos_tab_ <- ER_pos_tab %>% left_join(scale_meta_pos, by = "PERSONID")

scale_meta_pos_ER <- as.data.frame(sapply(ER_pos_tab_[, 64:1004], as.numeric))
ER_neg_tab_ <- ER_neg_tab %>% left_join(scale_meta_pos, by = "PERSONID")
scale_meta_pos_ER_ <- as.data.frame(sapply(ER_neg_tab_[, 64:1004], as.numeric))

# PROSTEGRONE receptor positive
PR_pos <- covariate_all_$PR_statusIndex == 1 | covariate_all_$Status == 0
PR_neg <- covariate_all_$PR_statusIndex == 0 | covariate_all_$Status == 0

PR_pos_tab <- covariate_all_[PR_pos, ]
PR_neg_tab <- covariate_all_[PR_neg, ]

PR_pos_tab <- PR_pos_tab[duplicated(PR_pos_tab$MATCH) | duplicated(PR_pos_tab$MATCH, fromLast = T), ]
PR_neg_tab <- PR_neg_tab[duplicated(PR_neg_tab$MATCH) | duplicated(PR_neg_tab$MATCH, fromLast = T), ]

PR_pos_tab_ <- PR_pos_tab %>% left_join(scale_meta_pos, by = "PERSONID")
scale_meta_pos_PR <- as.data.frame(sapply(PR_pos_tab_[, 64:1004], as.numeric))

PR_neg_tab_ <- PR_neg_tab %>% left_join(scale_meta_pos, by = "PERSONID")
scale_meta_pos_PR_ <- as.data.frame(sapply(PR_neg_tab_[, 64:1004], as.numeric))

# HER2 receptor positive
HER2_pos <- covariate_all_$HER2_statusIndex == 1 | covariate_all_$Status == 0
HER2_neg <- covariate_all_$HER2_statusIndex == 0 | covariate_all_$Status == 0

HER2_pos_tab <- covariate_all_[HER2_pos, ]
HER2_neg_tab <- covariate_all_[HER2_neg, ]

HER2_pos_tab <- HER2_pos_tab[duplicated(HER2_pos_tab$MATCH) | duplicated(HER2_pos_tab$MATCH, fromLast = T), ]
HER2_neg_tab <- HER2_neg_tab[duplicated(HER2_neg_tab$MATCH) | duplicated(HER2_neg_tab$MATCH, fromLast = T), ]

HER2_pos_tab_ <- HER2_pos_tab %>% left_join(scale_meta_pos, by = "PERSONID")
scale_meta_pos_HER2 <- as.data.frame(sapply(HER2_pos_tab_[, 64:1004], as.numeric))

HER2_neg_tab_ <- HER2_neg_tab %>% left_join(scale_meta_pos, by = "PERSONID")
scale_meta_pos_HER2_ <- as.data.frame(sapply(HER2_neg_tab_[, 64:1004], as.numeric))

# BC risk models for untargated metabolites 
# CLR models to get odds ratios for metabolites: all, pre-menopausal only and post-menopausal only, and

mean_sd(scale_meta_pos_$"165.0791@2.0293136")
hist(scale_meta_pos_$"165.0791@2.0293136")

table(covariate_all_$Status, covariate_all_$V443_Menopausal)

# All subjects
fits1 <- apply(scale_meta_pos_, 2, function(x) clogit(as.numeric(Status) ~ x + strata(MATCH), data=covariate_all_))
fits2 <- apply(scale_meta_pos_, 2, function(x) clogit(as.numeric(Status) ~ x + ETH + WC + V443_Menopausal + Bmi + OCEver_1 + breastfed + 
                                                      ever_alcohol + parous_1 + VIH + strata(MATCH), data=covariate_all_))
fits3 <- apply(scale_meta_pos_, 2, function(x) clogit(as.numeric(Status) ~ x + ETH + WC + V443_Menopausal + Bmi + OCEver_1 + breastfed + 
                                                        ever_alcohol + parous_1 + eduCat_1 + smokingEver_1 + DMT + VIH + strata(MATCH), data=covariate_all_))

#Pre-menopausal
fits_prem <- apply(scale_meta_pos_[pre_, ], 2, function(x) clogit(as.numeric(Status) ~ x + ETH + WC + Bmi + OCEver_1 + breastfed + 
                                                                    ever_alcohol + parous_1 + VIH + strata(MATCH), data=covariate_all_[pre_, ]))
#Post-menopausal
fits_post <- apply(scale_meta_pos_[post_, ], 2, function(x) clogit(as.numeric(Status) ~ x + ETH + WC + Bmi + OCEver_1 + breastfed + 
                                                                     ever_alcohol + parous_1 + VIH + strata(MATCH), data=covariate_all_[post_, ]))
#ER Postive and Negative #
fits_0_ER_p <- apply(scale_meta_pos_ER, 2, function(x) clogit(as.numeric(Status) ~ x + ETH + WC + V443_Menopausal + Bmi + OCEver_1 + breastfed + 
                                                                ever_alcohol + parous_1 + VIH + strata(MATCH), data=ER_pos_tab_))
fits_0_ER_n <- apply(scale_meta_pos_ER_, 2, function(x) clogit(as.numeric(Status) ~ x + ETH + WC + V443_Menopausal + Bmi + OCEver_1 + breastfed + 
                                                                 ever_alcohol + parous_1 + VIH + strata(MATCH), data=ER_neg_tab_))

#PR Postive and Negative #
fits_0_PR_p <- apply(scale_meta_pos_PR, 2, function(x) clogit(as.numeric(Status) ~ x + ETH + WC + V443_Menopausal + Bmi + OCEver_1 + breastfed + 
                                                                ever_alcohol + parous_1 + VIH + strata(MATCH), data=PR_pos_tab_))
fits_0_PR_n <- apply(scale_meta_pos_PR_, 2, function(x) clogit(as.numeric(Status) ~ x + ETH + WC + V443_Menopausal + Bmi + OCEver_1 + breastfed + 
                                                                 ever_alcohol + parous_1 + VIH + strata(MATCH), data=PR_neg_tab_))

fits_0_HER2_p <- apply(scale_meta_pos_HER2, 2, function(x) clogit(as.numeric(Status) ~ x + ETH + WC + V443_Menopausal + Bmi + OCEver_1 + breastfed + 
                                                                    ever_alcohol + parous_1 + VIH + strata(MATCH), data=HER2_pos_tab_))
fits_0_HER2_n <- apply(scale_meta_pos_HER2_, 2, function(x) clogit(as.numeric(Status) ~ x + ETH + WC + V443_Menopausal + Bmi + OCEver_1 + breastfed + 
                                                                     ever_alcohol + parous_1 + VIH + strata(MATCH), data=HER2_neg_tab_))

# Tables for manuscript
# Generate tidy output table from models
tidy.output <- function(mod) {
  library(broom) 
  # Function to exponentiate and round
  df <- map_df(mod, tidy) %>% filter(term == "x") %>% cbind(Compound = names(mod)) %>%
    left_join(meta_names_1, by  = "Compound") %>%
    mutate(P.value = round(p.value, 5), FDR = round(p.adjust(p.value, method = "fdr"), 5))
  # Calculate FDR adjusted confidence intervals  
  fdr.hits <- sum(p.adjust(df$p.value, method = "fdr") < 0.05)
  alpha.raw <- 0.05
  alpha.fdr <- (fdr.hits * 0.05)/nrow(df)
  
  df$OR <- exp(df$estimate)
  df$ci.low <- exp( df$estimate - (df$std.error * qnorm(1 - (alpha.raw/2))) )
  df$ci.high <- exp( df$estimate + (df$std.error * qnorm(1 - (alpha.raw/2))) )
  
  
  df$ci.low.fdr <- exp( df$estimate - (df$std.error * qnorm(1 - (alpha.fdr/2))) )
  df$ci.high.fdr <- exp( df$estimate + (df$std.error * qnorm(1 - (alpha.fdr/2))) )
  
  if(fdr.hits > 0) {
    df$ci.low.fdr <- exp( df$estimate - (df$std.error * qnorm(1 - (alpha.fdr/2))) )
    df$ci.high.fdr <- exp( df$estimate + (df$std.error * qnorm(1 - (alpha.fdr/2))) )
  } else {
    df$ci.low.fdr <- df$ci.low
    df$ci.high.fdr <- df$ci.high
  }
  
  df <- df %>% mutate_at(c("OR", "ci.low", "ci.high", "ci.low.fdr", "ci.high.fdr"), ~round(., 2))
  
  # Paste CIs together
  df$ci95  <- paste("(", df$ci.low, "-", df$ci.high, ")", sep = "")
  df$ciFDR <- paste("(", df$ci.low.fdr, "-", df$ci.high.fdr, ")", sep = "")
  
  # Select columns and order
  output <- df %>% select(Compound = "display_name", "description", "OR",
                          "ci95", "ciFDR", "P.value", "FDR", "metabolite_name") %>% arrange(description)
  
}


all_1 <- tidy.output(fits1) # model 1
all_2 <- tidy.output(fits2) # model 2
all_3 <- tidy.output(fits3) # model 3

prem <- tidy.output(fits_prem) # premenopausal women only
postm <- tidy.output(fits_post)# postmenopausal women only

all_ER_p <- tidy.output(fits_0_ER_p) # ER positive
all_ER_n <- tidy.output(fits_0_ER_n)# ER negative
all_PR_p <- tidy.output(fits_0_PR_p) # PR positive
all_PR_n <- tidy.output(fits_0_PR_n) # PR negative
all_HER2_p <- tidy.output(fits_0_HER2_p) # HER2 positive
all_HER2_n <- tidy.output(fits_0_HER2_n) # HER2 negative


# Retain only metabolite groups with at least P.value < 0.05

table_final_1 <- bind_rows("All" = all_1, .id = "Group") %>%
  arrange(Group, -OR) %>% 
  group_by(Compound) %>% 
  filter(P.value < 0.05) %>%
  select(Compound, description, everything())
table_final_1

table_final_2 <- bind_rows("All" = all_2, .id = "Group") %>%
  arrange(Group, -OR) %>% 
  group_by(Compound) %>% 
  filter(P.value < 0.05) %>% 
  select(Compound, description, everything()) 
table_final_2

table_final_3 <- bind_rows("All" = all_3, .id = "Group") %>%
  arrange(Group, -OR) %>% 
  group_by(Compound) %>%
  filter(P.value < 0.05) %>% 
  select(Compound, description, everything()) 
table_final_3

# Retain only metabolite groups with at least one FDR < 0.05

tab_fdr_1 <- bind_rows("All" = all_1, .id = "Group") %>%
  arrange(Group, -OR) %>% 
  group_by(Compound) %>%
  filter(min(FDR) < 0.05) %>%
  select(Compound, description, everything()) 
tab_fdr_1

tab_fdr_2 <- bind_rows("All" = all_2, .id = "Group") %>%
  arrange(Group, -OR) %>% 
  group_by(Compound) %>%
  filter(min(FDR) < 0.05) %>%
  select(Compound, description, everything()) 
tab_fdr_2

tab_fdr_3 <- bind_rows("All" = all_3, .id = "Group") %>%
  arrange(Group, -OR) %>% 
  group_by(Compound) %>%
  filter(min(FDR) < 0.05) %>%
  select(Compound, description, everything()) 
tab_fdr_3

# Pre-menopausal
tab_prem_all <- bind_rows("PREMO" = prem, .id = "Group") %>%
  arrange(Group, -OR) %>% 
  group_by(Compound) %>%
  select(Compound, description, everything()) 
tab_prem_all

tab_fdr_prem <- bind_rows("PREMO" = prem, .id = "Group") %>%
  arrange(Group, -OR) %>% 
  group_by(Compound) %>%
  filter(min(FDR) < 0.05) %>%
  select(Compound, description, everything()) 
tab_fdr_prem

# Post-menopausal
tab_post_all <- bind_rows("POST" = postm, .id = "Group") %>%
  arrange(Group, -OR) %>% 
  group_by(Compound) %>%
  select(Compound, description, everything()) 
tab_post_all

tab_fdr_post <- bind_rows("POST" = postm, .id = "Group") %>%
  arrange(Group, -OR) %>% 
  group_by(Compound) %>%
  filter(min(FDR) < 0.05) %>%
  select(Compound, description, everything()) 
tab_fdr_post


# save only metabolite groups with at least one p-value < 0.05

library(readr)

write_csv(table_final_1, path = "/Tables - revised/TABLE_POS_M1.csv")
write_csv(table_final_2, path = "/Tables - revised/TABLE_POS_M2.csv")
write_csv(tab_fdr_2_, path = "/Tables - revised/TABLE_POS_M2_P.value.csv")
write_csv(table_final_3, path = "/Tables - revised/TABLE_POS_M3.csv")

write_csv(tab_post_all, path = "/Tables - revised/TABLE_POS_M2_POST.csv")
write_csv(tab_prem_all, path = "/Tables - revised/TABLE_POS_M2_PREM.csv")

# ER POSITIVE Retain only metabolite groups with at least P.value < 0.05
tab_ER_p <- bind_rows("ER Positive_all" = all_ER_p, .id = "Group") %>%
  arrange(Group, -OR) %>% 
  group_by(Compound) %>%
  select(Compound, description, everything()) 
tab_ER_p

# ER NEGATIVE Retain only metabolite groups with at least P.value < 0.05
tab_ER_n <- bind_rows("ER Negative_all" = all_ER_n, .id = "Group") %>%
  arrange(Group, -OR) %>% 
  group_by(Compound) %>%
  select(Compound, description, everything()) 
tab_ER_n

# PR POSITIVE Retain only metabolite groups with at least P.value < 0.05
tab_PR_p <- bind_rows("PR Positive_all" = all_PR_p, .id = "Group") %>%
  arrange(Group, -OR) %>% 
  group_by(Compound) %>%
  select(Compound, description, everything()) 
tab_PR_p

# PR NEGATIVE Retain only metabolite groups with at least P.value < 0.05
tab_PR_n <- bind_rows("PR Negative_all" = all_PR_n, .id = "Group") %>%
  arrange(Group, -OR) %>% 
  group_by(Compound) %>%
  select(Compound, description, everything()) 
tab_PR_n

# HER2 POSITIVE Retain only metabolite groups with at least P.value < 0.05
tab_HER2_p <- bind_rows("HER2 Positive_all" = all_HER2_p, .id = "Group") %>%
  arrange(Group, -OR) %>% 
  group_by(Compound) %>%
  select(Compound, description, everything()) 
tab_HER2_p

# HER2 NEGATIVE Retain only metabolite groups with at least P.value < 0.05
tab_HER2_n <- bind_rows("HER2 Negative_all" = all_HER2_n, .id = "Group") %>%
  arrange(Group, -OR) %>% 
  group_by(Compound) %>%
  select(Compound, description, everything()) 
tab_HER2_n

write_csv(tab_ER_p, path = "/Tables - revised/TABLE_POS_M2_ER_p.csv")
write_csv(tab_ER_n, path = "/Tables - revised/TABLE_POS_M2_ER_n.csv")

write_csv(tab_PR_p, path = "/Tables - revised/TABLE_POS_M2_PR_p.csv")
write_csv(tab_PR_n, path = "/Tables - revised/TABLE_POS_M2_PR_n.csv")

write_csv(tab_HER2_p, path = "/Tables - revised/TABLE_POS_M2_HER2_p.csv")
write_csv(tab_HER2_n, path = "/Tables - revised/TABLE_POS_M2_HER2_n.csv")

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

write_csv(p_interaction_pos_ER, path = "/Tables - revised/P_interaction_ER.csv")
write_csv(p_interaction_pos_PR, path = "/Tables - revised/P_interaction_PR.csv")
write_csv(p_interaction_pos_HER2, path = "/Tables - revised/P_interaction_HER2.csv")

