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
# For MR this results in *no* SNPs available for analysis.
# 
# So did I use the correct manifest file, because this makes even less sense/genome build problems?
# I used hg19 in UCSC and hg19 is the manifest file based on, 
# but on what build is the mQTL analysis based?

#set directories and load packages
setwd("~/Documents/Studie/Major/Artikel Scriptie/DATA/v19 mQTLs")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

library(plyr)
library(readr)

#read mQTL data

mydir = "mQTLs_v20"
myfiles = list.files(path=mydir, pattern="*.txt", full.names=TRUE)

mQTLs_shear_stress_v20 = ldply(myfiles, read_csv)


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

#TEST SOMETHING: filter out shear stress mQTLs based on CpG sites

#mQTLs_shear_stress_v19 <- read_delim("mQTLs_shear_stress_v19.csv", 
#                                     ";", escape_double = FALSE, trim_ws = TRUE)

TEST <- mQTLs_shear_stress_v20[mQTLs_shear_stress_v20$ProbeID %in% CpGs_region[[1]], ]

write.csv(TEST, file = "mQTLs_shear_stress_v19_filtered_Illumina.csv", row.names = F)


#TEST SOMETHING: filter based on CpG location

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

TEST2 <- rbind(HOXA5, TMEM184B, ADAMTSL5, KLF4, KLF3, CMKLR1, PKP4, ACVRL1, DOK4, SPRY2, ENOSF1)

#same results








#inlezen data

CPG <- read.delim("CPG.txt", header = FALSE)


Illumina450K <-  read_csv("~/Documents/Studie/Major/Artikel Scriptie/Scripts/HumanMethylation450_15017482_v1-2.csv")

#mijn 11 genen eruit halen

IlluminaTEST <- Illumina450K[grep("KLF4|HOXA5|TMEM184B|ADAMTSL5|KLF3|CMKLR1|PKP4|ACVRL1|DOK4|SPRY2|ENOSF1", 
                                  Illumina450K$MAPINFO), ]


Illumina_old <- IlluminaTEST[-grep("Body|3'UTR", IlluminaTEST$UCSC_RefGene_Group), ]

ROI <- data.frame(Illumina_old$IlmnID, Illumina_old$UCSC_RefGene_Name, Illumina_old$Chromosome_36, Illumina_old$Coordinate_36, Illumina_old$Coordinate_36, Illumina_old$Strand)





#Range, lead, trait column erbij

range <- numeric(157)

ROI["range"] <- range

ROI$range <- 250000


lead <- numeric(157)
ROI["lead"] <- lead

ROI$lead <- "new_lead"

trait <- numeric(157)
ROI["trait"] <- trait

ROI$trait <- "SS"

#eerste regel goed maken

ROI[1, 8] <- "previous_&_new_lead"

#to and from column erbij

ROI$from <- ROI$Illumina_data.Coordinate_36

ROI$to <- ROI$Illumina_data.Coordinate_36

#volgorde goed

ROI2 <- ROI[c(1, 2, 3, 4, 10, 11, 7, 8, 9)]

#from to change values to TSS and -/+ 250kb

TEST <- ROI2

TEST$from[TEST$from %like% "^1590"] <- "158771722"
TEST$from[TEST$from %like% "^383"] <- "38092185"
TEST$from[TEST$from %like% "^271"] <- "27149812"
TEST$from[TEST$from %like% "^109"] <- "109291576"
TEST$from[TEST$from %like% "^107"] <- "107257216" 
TEST$from[TEST$from %like% "^505"] <- "50342380"
TEST$from[TEST$from %like% "^798"] <- "79813087"
TEST$from[TEST$from %like% "^560"] <- "56077886"
TEST$from[TEST$from %like% "^702"] <- "702676" 
TEST$from[TEST$from %like% "^703"] <- "702676" 
TEST$from[TEST$from %like% "^146"] <- "1464019" 
TEST$from[TEST$from %like% "^369"] <- "36998962"


TEST$to[TEST$to %like% "^1590"] <- "159021722" 
TEST$to[TEST$to %like% "^383"] <- "38342185" 
TEST$to[TEST$to %like% "^271"] <- "27399812" 
TEST$to[TEST$to %like% "^109"] <- "109541576" 
TEST$to[TEST$to %like% "^107"] <- "107507216" 
TEST$to[TEST$to %like% "^505"] <- "50592380" 
TEST$to[TEST$to %like% "^798"] <- "80063087" 
TEST$to[TEST$to %like% "^560"] <- "56327886" 
TEST$to[TEST$to %like% "^702"] <- "952676" 
TEST$to[TEST$to %like% "^703"] <- "952676" 
TEST$to[TEST$to %like% "^146"] <- "1714019" 
TEST$to[TEST$to %like% "^369"] <- "37248962" 

ROI2 <- TEST


#export

write.table(ROI2, file = "ROI", sep = " ", row.names = FALSE, col.names = FALSE)



#stuff

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("FDb.InfiniumMethylation.hg19")
#BiocManager::install("DNAmArray")

library(FDb.InfiniumMethylation.hg19)
#library(DNAmArray) #does not work in my version of R


# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# install.packages("reticulate")
# 
# BiocManager::install("IlluminaHumanMethylation450kmanifest")

library(reticulate)
library(IlluminaHumanMethylation450kmanifest)

data(IlluminaHumanMethylation450kmanifest)

head(IlluminaHumanMethylation450kmanifest)


