# Apparently the manifest file of Illumina is not correct, 
# therefore I used a R library to determine the CpGs of interest.
# First I did mQTL analysis with (approx.) 11 genes +/- 250kb and permutation.
# This results in 184 SNPs, but also SNPs outside of the promotor site of the genes
# (-2000TSS, 1st exon). 
# 
# I now made a list of -2000TSS -- end 1st exon of each gene using UCSC manually.
# Followed by this R library and HM450K locations of CpGs and filtered on all the CpGs within
# these regions. 
# 
# Followed by filtering out all mQTLs not related to one of those CpGs_mQTLs
# 
# This results in a mere 14 SNPs. 
# 
# For MR this results in *too few* SNPs available for analysis (from 18 & 14 SNPs to 1 & 0 SNP)
# 
# So did I use the correct manifest file, because this makes even less sense/genome build problems?
# I used hg19 in UCSC and hg19 is the manifest file based on, 
# but on what build is the mQTL analysis based?

#set directories and load packages
#setwd("~/Documents/Studie/Major/Artikel Scriptie/DATA/v19 mQTLs")

setwd("~/Genetics/Artikel Scriptie/DATA/v19 mQTLs") #windows

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

#load data
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)

data(Locations)

as.data.frame(Locations)

data(Islands.UCSC)

as.data.frame(Islands.UCSC)

Illumina_data <- as.data.frame(Locations)

#read the table I made (by hand...) using the UCSC browser Human Feb. 2009 (GRCh37/hg19)
library(readr)
Gene_regions <- read_delim("Gene_regions.csv", 
                           ";", escape_double = FALSE, trim_ws = TRUE)

#extract CpG sites (quick and neat :), not so efficient)

library(tidyverse)

HOXA5 <- subset(Illumina_data, pos > 27182665 & pos < 27185287 & chr == "chr7")
TMEM184B <- subset(Illumina_data, pos > 38668890 & pos < 38671040 & chr == "chr16")
ADAMTSL5 <- subset(Illumina_data, pos > 1512962 & pos < 1515188 & chr == "chr19")
KLF4 <- subset(Illumina_data, pos > 110251449 & pos < 110253927 & chr == "chr9")
KLF3 <- subset(Illumina_data, pos > 38663817 & pos < 38666082 & chr == "chr4")
CMKLR1 <- subset(Illumina_data, pos > 108732804 & pos < 108735094 & chr == "chr12")
PKP4 <- subset(Illumina_data, pos > 159311611 & pos < 159313730 & chr == "chr2")
ACVRL1 <- subset(Illumina_data, pos > 52298202 & pos < 52301479 & chr == "chr12")
DOK4 <- subset(Illumina_data, pos > 57520217 & pos < 57522407 & chr == "chr16")
SPRY2 <- subset(Illumina_data, pos > 80914757 & pos < 80917086 & chr == "chr13")
ENOSF1 <- subset(Illumina_data, pos > 712504 & pos < 715676 & chr == "chr18")

CpGs_mQTLs <- rbind(HOXA5, TMEM184B, ADAMTSL5, KLF4, KLF3, CMKLR1, PKP4, ACVRL1, DOK4, SPRY2, ENOSF1)

#order the data
CpGs_mQTLs$CpG <- rownames(CpGs_mQTLs)
rownames(CpGs_mQTLs) <- c()

CpGs_mQTLs <- CpGs_mQTLs[order(CpGs_mQTLs$chr), ]

CpGs_mQTLs <- CpGs_mQTLs[, c(4, 3, 1, 2)]

#list of CpGs that are in my regions
CpGs_region <- list(unique(CpGs_mQTLs$CpG))

#read mQTL data
library(plyr)

mydir = "mQTLs_v20/Nominal results"
myfiles = list.files(path=mydir, pattern="*.txt", full.names=TRUE)

mQTLs_shear_stress_nominal_v20 = ldply(myfiles, read_csv)

#filter out shear stress mQTLs based on CpG sites

shear_stress_mQTLs <- mQTLs_shear_stress_nominal_v20[mQTLs_shear_stress_nominal_v20$ProbeID %in% CpGs_region[[1]], ]

write.csv(shear_stress_mQTLs, file = "mQTLs_shear_stress_nominal_v20_filtered.csv", row.names = F)


#Manual input to verify results: filter based on CpG location

HOXA5 <- subset(mQTLs_shear_stress_v19, BP_CpG > 27182665 & BP_CpG < 27185287 & Chr_CpG == "7")
TMEM184B <- subset(mQTLs_shear_stress_v19, BP_CpG > 38668890 & BP_CpG < 38671040 & Chr_CpG == "16")
ADAMTSL5 <- subset(mQTLs_shear_stress_v19, BP_CpG > 1512962 & BP_CpG < 1515188 & Chr_CpG == "19")
KLF4 <- subset(mQTLs_shear_stress_v19, BP_CpG > 110251449 & BP_CpG < 110253927 & Chr_CpG == "9")
KLF3 <- subset(mQTLs_shear_stress_v19, BP_CpG > 38663817 & BP_CpG < 38666082 & Chr_CpG == "4")
CMKLR1 <- subset(mQTLs_shear_stress_v19, BP_CpG > 108732804 & BP_CpG < 108735094 & Chr_CpG == "12")
PKP4 <- subset(mQTLs_shear_stress_v19, BP_CpG > 159311611 & BP_CpG < 159313730 & Chr_CpG == "2")
ACVRL1 <- subset(mQTLs_shear_stress_v19, BP_CpG > 52298202 & BP_CpG < 52301479 & Chr_CpG == "12")
DOK4 <- subset(mQTLs_shear_stress_v19, BP_CpG > 57520217 & BP_CpG < 57522407 & Chr_CpG == "16")
SPRY2 <- subset(mQTLs_shear_stress_v19, BP_CpG > 80914757 & BP_CpG < 80917086 & Chr_CpG == "13")
ENOSF1 <- subset(mQTLs_shear_stress_v19, BP_CpG > 712504 & BP_CpG < 715676 & Chr_CpG == "18")

TEST <- rbind(HOXA5, TMEM184B, ADAMTSL5, KLF4, KLF3, CMKLR1, PKP4, ACVRL1, DOK4, SPRY2, ENOSF1)

#same results
