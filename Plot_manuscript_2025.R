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


#Selection of the annotated metabolites
load("XXX/Tables/scale_meta_neg.RData")

#dataset contains metabolomics data (negative mode), which were previously scaled and log2-transformed.
scale_meta_neg <- setDT(scale_meta_neg, keep.rownames=TRUE)[]
scale_meta_neg <- scale_meta_neg %>% rename (PERSONID=rn)

scale_meta <- scale_meta_pos %>% left_join(scale_meta_neg, by="PERSONID")
dim(scale_meta)

data_pos_all_CT <- data_pos_all %>% select(PERSONID, Status) %>% filter(Status==0)
scale_meta_0 <- scale_meta %>% right_join(data_pos_all_CT, by="PERSONID")

names(scale_meta_0)
dim(scale_meta_0)
scale_meta_ <- scale_meta_0 %>% select(2:1733)
dim(scale_meta_)


scale_meta_all <- scale_meta_ %>% rename("Kynurenine (related ion)"="93.0579@1.8251286",
                                         "Kynurenine (ammonia loss)"="191.0584@1.8263695",
                                         "Kynurenine (main ion)"="208.0849@1.8271979",
                                         "Octenoylcarnitine"="285.1948@3.955888",
                                         "Cortisol"="362.2101@5.142576",
                                         "Unknown_1_a"="269.1109@0.62916785",
                                         "Unknown_2_a"="356.2929@7.2220035",
                                         "Unknown_2_b"="373.3193@7.22367",
                                         "Unknown_3"="78.0141@0.85349643",
                                         "Unknown_4"="151.0633@1.8031964",
                                         "Unknown_5"="376.2596@7.071281",
                                         "Cortisol (related ion)"="494.1618@5.1330748",
                                         "Cortisol (main ion, FA adduct)"="408.2156@5.1332083",
                                         "Unknown_1_b"="269.1122@0.62074476",
                                         "Unknown_6_a"="108.0431@0.6281819",
                                         "Unknown_6_b"="162.053@0.6281759",
                                         "Unknown_6_c"="72.0228@0.6281736",
                                         "Unknown_7"="305.0789@0.6263058",
                                         "Unknown_8"="402.2985@7.2131867",
                                         "Unknown_9"="564.1505@0.6536153",
                                         "Unknown_11"="60.0222@0.6277415",
                                         "Unknown_10"="591.3543@6.978287",
                                         "Unknown_12"="662.4191@7.3528366") %>% 
                                    select("Kynurenine (related ion)",
                                           "Kynurenine (ammonia loss)",
                                           "Kynurenine (main ion)",
                                           "Octenoylcarnitine",
                                           "Cortisol",
                                           "Cortisol (related ion)",
                                           "Cortisol (main ion, FA adduct)",
                                           "Unknown_1_a",
                                           "Unknown_1_b",
                                           "Unknown_2_a",
                                           "Unknown_2_b",
                                           "Unknown_3",
                                           "Unknown_4",
                                           "Unknown_5",
                                           "Unknown_6_a",
                                           "Unknown_6_b",
                                           "Unknown_6_c",
                                           "Unknown_7",
                                           "Unknown_8",
                                           "Unknown_9",
                                           "Unknown_11",
                                           "Unknown_10",
                                           "Unknown_12")


cormat_feat_ <- cor(scale_meta_all, method = c("pearson"), use='pairwise.complete.obs')
library(ggcorrplot)
library(corrplot)

#CORRELATION BETWEEN THE ANNOTATED FEATURES (figure 2)

cordf_feat_ <- as_tibble(cormat_feat_)
plot_feat_A  <- ggcorrplot(cordf_feat_, hc.order = T, hc.method = "ward", insig = "blank", legend.title = "Pearson\ncorrelation") + 
  theme_minimal() + scale_x_continuous(expand = c(-1,1), position="top") + ggtitle("") +
  theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab("")

plot_feat_A_  <- ggcorrplot(cordf_feat_, hc.order = F, hc.method = "ward", legend.title = "Pearson\ncorrelation") + 
  theme_minimal() + scale_x_continuous(expand = c(0,0)) + ggtitle("") +
  theme(axis.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab("")

plot_feat_A_0  <- ggcorrplot(cordf_feat_, hc.order = T, hc.method = "ward", legend.title = "Pearson\ncorrelation") + 
  theme_minimal() + scale_x_continuous(expand = c(0,0)) + ggtitle("") +
  theme(axis.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab("")

corrplot(cormat_feat_, method = "shade", tl.col = "black", tl.cex = 1,  tl.srt = 90,
         hclust.method = "ward.D2", order = "original", type = "upper") 



scale_meta_all_ <- scale_meta_ %>% rename("Kynurenine (related ion)"="93.0579@1.8251286",
                                         "Kynurenine (ammonia loss)"="191.0584@1.8263695",
                                         "Kynurenine (main ion)"="208.0849@1.8271979",
                                         "Octenoylcarnitine"="285.1948@3.955888",
                                         "Cortisol"="362.2101@5.142576",
                                         "Unknown_1_a"="269.1109@0.62916785",
                                         "Unknown_2_a"="356.2929@7.2220035",
                                         "Unknown_2_b"="373.3193@7.22367",
                                         "Unknown_3"="78.0141@0.85349643",
                                         "Unknown_4"="151.0633@1.8031964",
                                         "Unknown_5"="376.2596@7.071281",
                                         "Cortisol (related ion)"="494.1618@5.1330748",
                                         "Cortisol (main ion, FA adduct)"="408.2156@5.1332083",
                                         "Unknown_1_b"="269.1122@0.62074476",
                                         "Unknown_6_a"="108.0431@0.6281819",
                                         "Unknown_6_b"="162.053@0.6281759",
                                         "Unknown_6_c"="72.0228@0.6281736",
                                         "Unknown_7"="305.0789@0.6263058",
                                         "Unknown_8"="402.2985@7.2131867",
                                         "Unknown_9"="564.1505@0.6536153",
                                         "Unknown_11"="60.0222@0.6277415",
                                         "Unknown_10"="591.3543@6.978287",
                                         "Unknown_12"="662.4191@7.3528366")


library(reticulate)
library(janitor)

meta_names_all <- scale_meta_all_[1,]
meta_names_all_ <- t(meta_names_all)
meta_names_all_0 <- as.data.frame(meta_names_all_)
meta_names_all_1 <- setDT(meta_names_all_0, keep.rownames=TRUE)[]
meta_names_all_1 <- meta_names_all_1 %>% rename (Compound=rn)
meta_names_all_1 <- meta_names_all_1 %>% select(Compound)

meta_names_all_1 <- meta_names_all_1 %>% mutate(display_name=Compound, description=Compound)


dim(meta_names_all_1)
names(meta_names_all_1)


scale_meta_f <- scale_meta %>% select(2:1733)
dim(scale_meta_f)


scale_meta_all_f <- scale_meta_f %>% rename("Kynurenine (related ion)"="93.0579@1.8251286",
                                         "Kynurenine (ammonia loss)"="191.0584@1.8263695",
                                         "Kynurenine (main ion)"="208.0849@1.8271979",
                                         "Octenoylcarnitine"="285.1948@3.955888",
                                         "Cortisol"="362.2101@5.142576",
                                         "Unknown_1_a"="269.1109@0.62916785",
                                         "Unknown_2_a"="356.2929@7.2220035",
                                         "Unknown_2_b"="373.3193@7.22367",
                                         "Unknown_3"="78.0141@0.85349643",
                                         "Unknown_4"="151.0633@1.8031964",
                                         "Unknown_5"="376.2596@7.071281",
                                         "Cortisol (related ion)"="494.1618@5.1330748",
                                         "Cortisol (main ion, FA adduct)"="408.2156@5.1332083",
                                         "Unknown_1_b"="269.1122@0.62074476",
                                         "Unknown_6_a"="108.0431@0.6281819",
                                         "Unknown_6_b"="162.053@0.6281759",
                                         "Unknown_6_c"="72.0228@0.6281736",
                                         "Unknown_7"="305.0789@0.6263058",
                                         "Unknown_8"="402.2985@7.2131867",
                                         "Unknown_9"="564.1505@0.6536153",
                                         "Unknown_11"="60.0222@0.6277415",
                                         "Unknown_10"="591.3543@6.978287",
                                         "Unknown_12"="662.4191@7.3528366")


meta = c("Kynurenine (related ion)",
       "Kynurenine (ammonia loss)",
       "Kynurenine (main ion)",
       "Octenoylcarnitine",
       "Cortisol",
       "Cortisol (related ion)",
       "Cortisol (main ion, FA adduct)",
       "Unknown_1_a",
       "Unknown_1_b",
       "Unknown_2_a",
       "Unknown_2_b",
       "Unknown_3",
       "Unknown_4",
       "Unknown_5",
       "Unknown_6_a",
       "Unknown_6_b",
       "Unknown_6_c",
       "Unknown_7",
       "Unknown_8",
       "Unknown_9",
       "Unknown_11",
       "Unknown_10",
       "Unknown_12")

#fits2 <- apply(scale_meta_pos_, 2, function(x) clogit(as.numeric(Status) ~ x + ETH + WC + V443_Menopausal + Bmi + OCEver_1 + breastfed + 
                                                        #ever_alcohol + parous_1 + VIH + strata(MATCH), data=covariate_all_))

mclr1_0_all <- function(x, dat) {clogit(as.numeric(Status) ~ x + ETH + WC + V443_Menopausal + Bmi + OCEver_1 + breastfed + 
                                          ever_alcohol + parous_1 + VIH + strata(MATCH), data = dat)}


#all_0_all_f <- meta_names_all_1 %>% filter(description %in% meta)

names(meta_names_all_1)

m1_0_all <- apply(scale_meta_all_f, 2, mclr1_0_all, dat = covariate_all_) %>% map_df(tidy, exponentiate = T) %>%
  filter(term == "x") %>% mutate(p.adj = p.adjust(p.value, method = "fdr")) %>% bind_cols(meta_names_all_1)

all_0_all <- bind_rows("All Women" = m1_0_all, .id = "analysis") %>% 
  filter(term == "x") %>% mutate(p.adj = p.adjust(p.value, method = "fdr"))

all_0_all$sign <- ifelse (all_0_all$p.adj < 0.05,"red","grey")
all_0_all$sign[all_0_all$p.adj < 0.05 & all_0_all$estimate < 1] <- "blue"


all_0_all_f <- all_0_all %>% filter(description %in% meta)

#Figure 1 manuscript

ggplot(all_0_all, aes(x = estimate, y = -log10(p.value), color = sign, label = `description`)) +
        geom_point(aes(alpha = ifelse(p.value < 0.05, 1, 0.01) ),size = 4, show.legend = F) +
        xlab(expression("Odds ratio per SD increase")) +
        ylab(expression("-Log"[10]*" P value")) +
        labs(color = "Threshold\n") +
        geom_hline(
        yintercept = c(-log10(0.05), -log10(0.000316)), #yintercept = log10(0.001), linetype = "dotted"
        col = "black",
        linetype = c("dotted","dotdash"),
        size = 1)+
        theme_bw() +
        theme(text =  element_text(size=16), axis.text=element_text(size=16))+
        scale_size(range = c(1.75, 2.50))+
        scale_alpha(range = c(0.1, 1.00))+ 
        #theme(legend.position = "none")+
        scale_colour_manual(values = c("red", "blue", "red")) +
        geom_label_repel(data=subset(all_0_all_f, p.value < 0.05),
        aes(label = `description`, 
        fontface = "bold"),size = 4.3, color="steelblue") + 
        labs(title = "", subtitle =  "All women") + 
        annotate("text", x = 0.3, y=c(1.45, 3.65), size = 4.5, 
           label = c(expression(paste("Raw ", italic(P),"-threshold" )), 
                     expression(paste("FDR ", italic(P),"-threshold" )))) + xlim(0.1, 2)
        
# Plot faceting by analysis
library(ggplot2)
ggplot(all_0_all, aes(x = estimate, y = -log10(p.value))) + geom_point(shape = 1) +
  facet_wrap(fct_inorder(analysis) ~ ., scales = "free_x") + theme_bw() +
  geom_vline(xintercept = 1, size = 0.2, colour = "grey60") + 
  geom_hline(yintercept = -log10(0.05), size = 0.2, colour = "grey60")


# Create base plot to cut down code
# All participants
base1_0_all <- ggplot(m1_0_all, aes((estimate), log10(p.value))) + geom_point(shape = 1) + 
  theme_bw(base_size = 10) +
  xlab("Odds ratio per SD increase concentration") + 
  ylab(expression(paste(italic(P), "-value"))) +
  geom_vline(xintercept = 1, size = 0.2, colour = "grey60") + 
  geom_hline(yintercept = log10(0.05), size = 0.2, colour = "grey60") +
  theme(panel.grid.minor = element_blank() , panel.grid.major = element_blank(),
        axis.title.x = element_text(size = 9))

p1_all <- base1_0_all %+% xlim(0.5, 2) +
  scale_y_reverse(limits = c(0, -8), breaks = c(0:-5), labels = function(x) 10^x) +
  geom_text_repel(aes(label = description), size = 3, data = m1_0_all[m1_0_all$p.value < 0.05, ]) + 
  geom_hline(yintercept = log10(0.001), linetype = "dotted") +
  labs(title = "", subtitle =  "All participants") +
  annotate("text", x = 0.85, y=-1.37, size = 3, 
           label = expression(paste("Raw ", italic(P),"-threshold" )))

p1_all_ <- base1_0_all %+% xlim(0.2, 2.0) +
  scale_y_reverse(limits = c(0, -8), breaks = c(0:-5), labels = function(x) 10^x) +
  geom_text_repel(aes(label = description), size = 3.5, data = m1_0_all[m1_0_all$p.adj < 0.05, ]) + 
  geom_hline(yintercept = log10(0.001), linetype = "dotted") +
  labs(title = "", subtitle =  "All participants") +
  annotate("text", x = 0.5, y=c(-1.4, -3.1), size = 3, 
           label = c(expression(paste("Raw ", italic(P),"-threshold" )), 
           expression(paste("FDR ", italic(P),"-threshold" ))))