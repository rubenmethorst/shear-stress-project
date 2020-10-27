#setwd("~/Documents/Studie/Major/Artikel Scriptie/DATA/v19 mQTLs/") #mac
setwd("~/Genetics/Artikel Scriptie/DATA/v19 mQTLs") #windows

library(readr)
library(TwoSampleMR)
library(MRInstruments)
library(ggplot2)

#read mQTL data
mQTLs_shear_stress <- read_csv("~/Genetics/Artikel Scriptie/DATA/v19 mQTLs/mQTLs_shear_stress_v19_filtered_Illumina_v2.csv") #windows

#unify RS-identifiers for mQTL data

mQTLs_shear_stress$VARIANT <- as.character(mQTLs_shear_stress$VARIANT)

has.rsid = grep(mQTLs_shear_stress$VARIANT, pattern="rs")
hasnt.rsid = setdiff(1:nrow(mQTLs_shear_stress), has.rsid)


mQTLs_shear_stress[has.rsid, "VARIANT"] <-  gsub(x = unlist(mQTLs_shear_stress[has.rsid, "VARIANT"]), pattern = "(rs\\d+):.*$", replacement = "\\1")
mQTLs_shear_stress[hasnt.rsid, "VARIANT"] <- gsub(x = unlist(mQTLs_shear_stress[hasnt.rsid, "VARIANT"]), pattern = "(:[ACTG].*$)", replacement = "")

uniqueSNPs <- list(unique(mQTLs_shear_stress$VARIANT))

# List available GWASs
ao <- available_outcomes()

data(gwas_catalog)
head(gwas_catalog)

#read mQTL data as MR

write.csv(mQTLs_shear_stress, file = "mQTLs_shear_stress_v19_input.csv", row.names = F)

mQTL_exposure_data <- read_exposure_data(
  filename = "mQTLs_shear_stress_v19_input.csv",
  sep = ",",
  snp_col = "VARIANT",
  beta_col = "Beta",
  se_col = "SE",
  effect_allele_col = "CodedAlleleA",
  other_allele_col = "OtherAlleleB",
  eaf_col = "MAF",
  pval_col = "Perm_P",
  gene_col = "GeneName",
  samplesize_col = "N"
)

mQTL_exposure_data$exposure <- "DNAm"

##MR on Coronary Artery Disease

#clumping

mQTL_exposure_data_clumped <- clump_data(mQTL_exposure_data)

#extract outcome data

out_dat_7 <- extract_outcome_data(
  snps = mQTL_exposure_data_clumped$SNP,                  
  outcomes = 7,                                           #enter your MR Base outcomes here
  proxies = TRUE
)

#harmonize data

mQTL_vs_out_7 <- harmonise_data(
  exposure_dat = mQTL_exposure_data_clumped,
  outcome_dat = out_dat_7
)


#MR analysis

mQTL_res_out_7 <- mr(mQTL_vs_out_7)

res_single_out_7 <- mr_singlesnp(mQTL_vs_out_7)

res_pleiotropy_7 <- mr_pleiotropy_test(mQTL_vs_out_7)

write.csv(mQTL_res_out_7, file = "mQTL_res_out_7_20200210.csv", row.names = F)
write.csv(res_single_out_7, file = "res_single_out_7_20200210.csv", row.names = F)
write.csv(res_pleiotropy_7, file = "res_pleiotropy_7_20200210.csv", row.names = F)

#PLOTS

PLOT_7 <- mr_scatter_plot(mQTL_res_out_7, mQTL_vs_out_7)

PLOT_7 <- PLOT_7[[paste0(mQTL_res_out_7$id.exposure[1], ".", mQTL_res_out_7$id.outcome[1])]] + #COPY THE ID FROM THE MR TABLE!
  theme_minimal() + 
  theme(legend.position="bottom",
        panel.background = element_rect(fill = "#FFFFFF",
                                        colour = "#FFFFFF"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  labs(x = "SNP effect on DNA methylation") + 
  labs(y = "SNP effect on Coronary Heart Disease") +
  scale_color_manual(values=c("#5EB17F","#F59D10","#E55738","#8D5B9A","#E35493"))

plot(PLOT_7)

#forestplot

PLOT_7_forest <- mr_forest_plot(res_single_out_7)
PLOT_7_forest[[1]]


##MR on Ischemic stroke

#clumping

mQTL_exposure_data_clumped <- clump_data(mQTL_exposure_data)

#extract outcome data

out_dat_1108 <- extract_outcome_data(
  snps = mQTL_exposure_data_clumped$SNP,                  
  outcomes = 1108,                                      #enter your MR Base outcomes here
  proxies = TRUE
)

#harmonize data

mQTL_vs_out_1108 <- harmonise_data(
  exposure_dat = mQTL_exposure_data_clumped,
  outcome_dat = out_dat_1108
)


#MR analysis

mQTL_res_out_1108 <- mr(mQTL_vs_out_1108)

res_single_out_1108 <- mr_singlesnp(mQTL_vs_out_1108)

res_pleiotropy_1108 <- mr_pleiotropy_test(mQTL_vs_out_1108)

write.csv(mQTL_res_out_1108, file = "mQTL_res_out_1108_20200210.csv", row.names = F)
write.csv(res_single_out_1108, file = "res_single_out_1108_20200210.csv", row.names = F)
write.csv(res_pleiotropy_1108, file = "res_pleiotropy_1108.csv", row.names = F)

#PLOTS

PLOT_1108 <- mr_scatter_plot(mQTL_res_out_1108, mQTL_vs_out_1108)


PLOT_1108 <- PLOT_1108[["7GOYOu.1108"]] + #COPY THE ID FROM THE MR TABLE!
  theme_minimal() + 
  theme(legend.position="bottom",
        panel.background = element_rect(fill = "#FFFFFF",
                                        colour = "#FFFFFF"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  labs(x = "SNP effect on DNA methylation") + 
  labs(y = "SNP effect on Stroke") +
  scale_color_manual(values=c("#5EB17F","#F59D10","#E55738","#8D5B9A","#E35493"))

PLOT_1108

#forestplot

PLOT_1108_forest <- mr_forest_plot(res_single_out_1108)
PLOT_1108_forest[[1]]
