## THE SOUTH AFRICA BREAST CANCER (SABC) STUDY ##
## UNTARGETED METABOLOMICS AND BREAST CANCER RISK IN BLACK SOUTH AFRICAN WOMEN ##

## ________________LIST OF LIBRARY________________________##

## ---- eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------
library(broom)
library(survival)
library(haven)
library(tidyverse)
library(ggrepel)
library(data.table)

library(kableExtra)
library(qwraps2)
options(qwraps2_markup = "markdown")

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
covariate_all <- data_pos_all %>% select(1, 942:1004)


#CHECK AND RECODE COVARIABLES FOR LOGISTIC REGRESSION
#BODY MASS INDEX (BMI)

mean(covariate_all$Bmi, na.rm = TRUE)
str(covariate_all$Bmi)
which(is.na(covariate_all$Bmi)) #0

covariate_all$bmi_cat[covariate_all$Bmi<18.50]=1
covariate_all$bmi_cat[covariate_all$Bmi>=18.50 & covariate_all$Bmi<25.00]=2
covariate_all$bmi_cat[covariate_all$Bmi>=25.00 & covariate_all$Bmi<=29.99]=3
covariate_all$bmi_cat[covariate_all$Bmi>=30.0]=4

table(covariate_all$bmi_cat, covariate_all$Status)

covariate_all$bmi_cat_1[covariate_all$Bmi<18.50]=1
covariate_all$bmi_cat_1[covariate_all$Bmi>=18.50 & covariate_all$Bmi<25.00]=2
covariate_all$bmi_cat_1[covariate_all$Bmi>=25.00 & covariate_all$Bmi<=29.99]=3
covariate_all$bmi_cat_1[covariate_all$Bmi>=30.00 & covariate_all$Bmi<35.00]=4
covariate_all$bmi_cat_1[covariate_all$Bmi>=35.00]=5

table(covariate_all$bmi_cat_1, covariate_all$Status)

# PARITY
table(covariate_all$Parity)
covariate_all$Parity_1[covariate_all$Parity==0]=0
covariate_all$Parity_1[covariate_all$Parity==1 | covariate_all$Parity==2]=1
covariate_all$Parity_1[covariate_all$Parity==3 | covariate_all$Parity==4 | covariate_all$Parity==5]=2
covariate_all$Parity_1[covariate_all$Parity==6 | covariate_all$Parity==7 | covariate_all$Parity==8]=3
covariate_all$Parity_1[covariate_all$Parity==9 | covariate_all$Parity==10]=3

table(covariate_all$Parity_1, covariate_all$Status)

covariate_allp = covariate_all %>% filter(Parity==!0)
median(covariate_allp$Parity)


#HEIGHT
mean(covariate_all$Height, na.rm = TRUE)
which(is.na(covariate_all$Height)) # 0

mean(covariate_all$Height, na.rm = TRUE)
which(is.na(covariate_all$Height)) # 0 

#SMOKING
table(covariate_all$smokingEver, covariate_all$Status)
table(covariate_all$smokingEver, useNA = "ifany")
covariate_all$smokingEver_1[covariate_all$smokingEver==0]=1
covariate_all$smokingEver_1[covariate_all$smokingEver==1]=2
covariate_all$smokingEver_1[covariate_all$smokingEver==2]=2

table(covariate_all$smokingEver_1, covariate_all$Status)

#EDUCATION LEVEL
table(covariate_all$eduCat, useNA = "ifany")
table(covariate_all$eduCat, covariate_all$Status)

covariate_all$eduCat_1[covariate_all$eduCat<=1]=1
covariate_all$eduCat_1[covariate_all$eduCat==2]=2
covariate_all$eduCat_1[covariate_all$eduCat>2]=3
table(covariate_all$eduCat_1, covariate_all$Status)

#EVER BREAST FED IN PAROUS WOMEN 
table(covariate_all$breastfed, useNA = "ifany")
table(covariate_all$breastfed, covariate_all$Status)

#ORAL CONTRACEPTIVE 
table(covariate_all$OCEver, useNA = "ifany")
table(covariate_all$OCEver, covariate_all$Status)

covariate_all$OCEver_1[covariate_all$OCEver==0]=1
covariate_all$OCEver_1[covariate_all$OCEver==1]=2

table(covariate_all$OCEver_1, covariate_all$Status)

#HRT use (NOT USE IN THE MODEL BECAUSE LOW NUMBER)
table(covariate_all$HRTEver, useNA = "ifany")
table(covariate_all$HRTEver, covariate_all$Status)

#PAROUS
table(covariate_all$parous, useNA = "ifany")
table(covariate_all$parous, covariate_all$Status)
covariate_all$parous_1[covariate_all$parous==0]=1
covariate_all$parous_1[covariate_all$parous==1]=2
table(covariate_all$parous_1, covariate_all$Status)

#ETHNICITY
table(covariate_all$ETHNICITY, useNA = "ifany")
table(covariate_all$ETHNICITY, covariate_all$Status)
str(covariate_all$ETHNICITY)

covariate_all$ETH[covariate_all$ETHNICITY==9]=4
covariate_all$ETH[covariate_all$ETHNICITY==6]=3
covariate_all$ETH[covariate_all$ETHNICITY==8]=3
covariate_all$ETH[covariate_all$ETHNICITY==3]=2
covariate_all$ETH[covariate_all$ETHNICITY==1]=1
covariate_all$ETH[covariate_all$ETHNICITY==2]=1
covariate_all$ETH[covariate_all$ETHNICITY==5]=1
covariate_all$ETH[covariate_all$ETHNICITY==7]=1
covariate_all$ETH[covariate_all$ETHNICITY==88]=1


#COMORBIDITIES
#HYPERTENSION
table(covariate_all$HP_diag, useNA = "ifany")
table(covariate_all$HP_diag, covariate_all$Status)

#DIABETES
table(covariate_all$DMT, useNA = "ifany")
table(covariate_all$DMT, covariate_all$Status)

#HIV STATUS
table(covariate_all$V188_Hiv_Positive, useNA = "ifany")
table(covariate_all$V188_Hiv_Positive, covariate_all$Status)

covariate_all$VIH[covariate_all$V188_Hiv_Positive==0]=0
covariate_all$VIH[covariate_all$V188_Hiv_Positive==1]=1
covariate_all$VIH[is.na(covariate_all$V188_Hiv_Positive)]=0
table(covariate_all$VIH, covariate_all$Status)

#ALCOHOL CONSUMPTION
table(covariate_all$ever_alcohol, useNA = "ifany")
table(covariate_all$ever_alcohol, covariate_all$Status)

#BREAST CANCER CARACTERISTIC
#Family history
table(covariate_all$FamHist, useNA = "ifany")
table(covariate_all$FamHist, covariate_all$Status)

covariate_all$mensAgeLast
median(covariate_all$mensAgeLast)
which(is.na(covariate_all$mensAgeLast)) # 

#HORMONE RECEPTOR STATUS
table(covariate_all$ER_statusIndex, useNA = "ifany")
table(covariate_all$PR_statusIndex, useNA = "ifany")
table(covariate_all$HER2_statusIndex, useNA = "ifany")
table(covariate_all$HER2_statusIndex, covariate_all$ER_statusIndex, covariate_all$PR_statusIndex)
table(covariate_all$ER_statusIndex, covariate_all$PR_statusIndex, useNA = "ifany")

# AGE AT MENARCHE
mean_sd(covariate_all$AgeMenarche)
median(covariate_all$AgeMenarche)

tapply(covariate_all$AgeMenarche, covariate_all$Status, mean)
quantile(covariate_all$AgeMenarche[covariate_all$Status==1], probs = c(.25, .5, .75))
quantile(covariate_all$AgeMenarche[covariate_all$Status==0], probs = c(.25, .5, .75))

# SUBSET FOR BREASTFED AND POSTMENOPAUSAL WOMEN
covariate_allbreastfed <- covariate_all %>% filter(breastfed == 1)
covariate_allpostm <- covariate_all %>% filter(MenoStat==2)

tapply(covariate_allpostm$mensAgeLast, covariate_allpostm$Status, median)
tapply(covariate_allbreastfed$breastMos, covariate_allbreastfed$Status, median)

quantile(covariate_allbreastfed$breastMos[covariate_allbreastfed$Status==0], probs = c(.25, .5, .75))
quantile(covariate_allbreastfed$breastMos[covariate_allbreastfed$Status==1], probs = c(.25, .5, .75))
quantile(covariate_all$mensAgeLast[covariate_all$Status==0], probs = c(.25, .5, .75))
quantile(covariate_all$mensAgeLast[covariate_all$Status==1], probs = c(.25, .5, .75))


sum <-
  list(#"Total subjects"                  = list("N"   =   ~ n()),
    
    "Age at blood collection (years)" = list("Mean"     =  ~ mean_sd(AgeInt, digits = 1)),
    "Menopausal status at blood collection" =
      
      list("Pre-menopausal"   =  ~ n_perc0(MenoStat == 1, na_rm = T, digits = 1),
           "Post-menopausal"  =  ~ n_perc0(MenoStat == 2, na_rm = T, digits = 1),
           "Unknown"  =  ~ n_perc0(MenoStat == 888, na_rm = T, digits = 1)),
    "Age at menarche" =
      
      list("Mean (SD)" = ~ mean_sd(AgeMenarche, digits = 1)),
    
    "ETHNICITY" =   
      
      list("Zulu/Pedi/Xhosa/Tswana/Swazi"  =  ~ n_perc0(ETH == 1, na_rm = T, digits = 1),
           "Sotho"       =  ~ n_perc0(ETH == 2, na_rm = T, digits = 1),
           "Venda/Tsonga"   =  ~ n_perc0(ETH == 3, na_rm = T, digits = 1),
           "Ndebele"        =  ~ n_perc0(ETH == 4, na_rm = T, digits = 1)),
    
    "BMI" =
      
      list("Mean (SD)" = ~ mean_sd(Bmi, digits = 1)),
    
    "BMI Level" =   
      
      list("Underweight"  =  ~ n_perc0(bmi_cat == 1, na_rm = T, digits = 1),
           "Normal"       =  ~ n_perc0(bmi_cat == 2, na_rm = T, digits = 1),
           "Overweight"   =  ~ n_perc0(bmi_cat == 3, na_rm = T, digits = 1),
           "Obese"        =  ~ n_perc0(bmi_cat == 4, na_rm = T, digits = 1)),
    
    "BMI Level 1" =   
      
      list("Underweight"  =  ~ n_perc0(bmi_cat_1 == 1, na_rm = T, digits = 1),
           "Normal"       =  ~ n_perc0(bmi_cat_1 == 2, na_rm = T, digits = 1),
           "Overweight"   =  ~ n_perc0(bmi_cat_1 == 3, na_rm = T, digits = 1),
           "Obese class I"        =  ~ n_perc0(bmi_cat_1 == 4, na_rm = T, digits = 1),
           "Obese class II or II"        =  ~ n_perc0(bmi_cat_1 == 5, na_rm = T, digits = 1)),
  
    "Height (cm)" =
      
      list("Mean (SD)" = ~ mean_sd(Height, digits = 1)),
    
    "WC (cm)" =
      
      list("Mean (SD)" = ~ mean_sd(WC, digits = 1)),
    
     "HC (cm)" =
      
      list("Mean (SD)" = ~ mean_sd(HC, digits = 1)),
    
    "Smoking status" = 
      
      list("Never" = ~ n_perc0(smokingEver== 0, na_rm = T, digits = 1),
           "Past"  = ~ n_perc0(smokingEver == 1, na_rm = T, digits = 1),
          "Current" = ~ n_perc0(smokingEver== 2, na_rm = T, digits = 1)),
    
    "EDUCATIONAL LEVELS" = 
      
      list("None/primary" = ~ n_perc0(eduCat_1 == 1, na_rm = T, digits = 1),
           "High school"  = ~ n_perc0(eduCat_1  == 2, na_rm = T, digits = 1),
           "College/University/postgraduate" = ~ n_perc0(eduCat_1  == 3, na_rm = T, digits = 1)),
    
    "PA categories" =   
      
      list(" Sedentary"  =  ~ n_perc0(PA_cat == 1, na_rm = T, digits = 1),
           "Lightly active"       =  ~ n_perc0(PA_cat == 2, na_rm = T, digits = 1),
           "Moderately active"   =  ~ n_perc0(PA_cat == 3, na_rm = T, digits = 1),
           "Very active"        =  ~ n_perc0(PA_cat == 4, na_rm = T, digits = 1)),
    
    "PA (METs/week)" =
      
      list("Mean (SD)" = ~ mean_sd(PA_total, digits = 1)),
    
   # "Alcohol intake (g/day)" =
      
    #  list("Mean (SD)" = ~ mean_sd(ALCOHOL, digits = 1)),
    
    "Ever breast fed in parous women " =
      
      list("No"  =  ~ n_perc0(breastfed == 0, na_rm = T, digits = 1),
           "Yes"    =  ~ n_perc0(breastfed == 1, na_rm = T, digits = 1)),
    
    "Duration of breast feeding (months)" =
      
      list("Mean (SD)" = ~ mean_sd(breastMos, digits = 1)),
    
    "Previous oral contraceptive use" = 
      
      list("Never" = ~ n_perc0(OCEver_1 == 1, na_rm = T, digits = 1),
           "Ever" = ~ n_perc0(OCEver_1 == 2, na_rm = T, digits = 1)),
    
    " Previous hormonal replacement therapy (HRT) use" =
      
      list("Never" =  ~ n_perc0(HRTEver == 0, na_rm = T, digits = 1),
           "Ever"  =  ~ n_perc0(HRTEver == 1, na_rm = T, digits = 1)),
    
    "Number of full-term pregnancies" =
      
      list("Nulliparous"  = ~ n_perc0(parous_1 == 1, na_rm = T, digits = 1),
           "Parous" = ~ n_perc0(parous_1 == 2, na_rm = T, digits = 1)),
   
   "Number of full-term pregnancies 2" =
     
     list("Nulliparous"  = ~ n_perc0(Parity_1 == 0, na_rm = T, digits = 1),
          "1-2 children" = ~ n_perc0(Parity_1 == 1, na_rm = T, digits = 1), 
          "3-5 children" = ~ n_perc0(Parity_1 == 2, na_rm = T, digits = 1),
          ">6 children" = ~ n_perc0(Parity_1 == 3, na_rm = T, digits = 1)),
   
   " Number of children in parous women" =
     
     list("Mean (SD)" = ~ mean_sd(as.numeric(Parity), digits = 1)),
   
    " Age at end of first full-term pregnancy" =
      
      list("Mean (SD)" = ~ mean_sd(AgeFFTP, digits = 1)),
   
    "Family history of breast cancer" =
      
      list("No"  = ~ n_perc0(FamHist == 0, na_rm = T, digits = 1),
           "Yes" = ~ n_perc0(FamHist == 1, na_rm = T, digits = 1)),
   
   "HIV positivity" =
     
     list("No"  = ~ n_perc0(VIH == 0, na_rm = T, digits = 1),
          "Yes" = ~ n_perc0(VIH == 1, na_rm = T, digits = 1)),
   
   "DIABETES" =
     
     list("No"  = ~ n_perc0(DMT == 0, na_rm = T, digits = 1),
          "Yes" = ~ n_perc0(DMT == 1, na_rm = T, digits = 1)),
   
   "HYPERTENSION" =
     
     list("No"  = ~ n_perc0(HP_diag == 0, na_rm = T, digits = 1),
          "Yes" = ~ n_perc0(HP_diag == 1, na_rm = T, digits = 1))
   )
# ---------------------------------------------------------------

# By case-control status
st <- data_pos_all_p %>% group_by(Status) %>% summary_table(sum)
print(st)
print(st, cnames = c("Controls (N=396)", "Cases (N=396)"))
class(st)
st1 <- data.frame(st)
st1 <- setDT(st, keep.rownames=TRUE)[]
write_csv(st, path = "XXXX/Table 1.csv")

# ---------------------------------------------------------------
sum2 <-
  list(#"Total subjects"                  = list("N"   =   ~ n()),
    
    "Age at diagnosis" =
      
      list("Mean (SD)" = ~ mean_sd(AgeDiagIndex, digits = 1)),
    
    "ER+" =
      
      list("No"  = ~ n_perc0(ER_statusIndex == 0, digits = 1),
           "Yes" = ~ n_perc0(ER_statusIndex == 1, digits = 1)),
    "PR+" =
      
      list("No"  = ~ n_perc0(PR_statusIndex == 0, digits = 1),
           "Yes" = ~ n_perc0(PR_statusIndex == 1, digits = 1)),
    "HER2+" =
      
      list("No"  = ~ n_perc0(HER2_statusIndex == 0, digits = 1),
           "Yes" = ~ n_perc0(HER2_statusIndex == 1, digits = 1),
           "Unknown" = ~ n_perc0(HER2_statusIndex == 888, digits = 1)),
  
    "Grade" =
      
      list("Well differentiated"  = ~ n_perc0(GradeIndex == 1, digits = 1),
           "Moderately differentiated" = ~ n_perc0(GradeIndex == 2, digits = 1),
           "Poorly/un-differentiated" = ~ n_perc0(GradeIndex == 3, digits = 1),
           "Unknown" = ~ n_perc0(GradeIndex == 888, digits = 1)))

# ---------------------------------------------------------------

st2 <- data_pos_all_ %>% group_by(Status) %>% 
  summary_table(summ2)
print(st2)

print(st2, cnames = c("Controls (N=396)", "Cases (N=396)"))
class(st2)

st2 <- data.frame(st2)
st2 <- setDT(st2, keep.rownames=TRUE)[]
write_csv(st2, path = "XXXX/Table_1.csv")


st3 <- data_pos_all_ %>% filter(Status == 1) %>% summary_table(summ2)
print(st3)
print(st3, cnames = c("Cases (N=396)"))