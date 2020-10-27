#read data

CPG <- read.delim("CPG.txt", header = FALSE)


Illumina450K <-  read_csv("~/Documents/Studie/Major/Artikel Scriptie/Scripts/HumanMethylation450_15017482_v1-2.csv")

#extract 11 genes of interest

Illumina450 <- Illumina450K[grep("KLF4|HOXA5|TMEM184B|ADAMTSL5|KLF3|CMKLR1|PKP4|ACVRL1|DOK4|SPRY2|ENOSF1", 
                                  Illumina450K$MAPINFO), ]


Illumina <- Illumina450[-grep("Body|3'UTR", Illumina450$UCSC_RefGene_Group), ]

ROI <- data.frame(Illumina$IlmnID, Illumina$UCSC_RefGene_Name, Illumina$Chromosome_36, Illumina$Coordinate_36, Illumina$Coordinate_36, Illumina$Strand)



#Range, lead, trait column added

range <- numeric(157)

ROI["range"] <- range

ROI$range <- 250000


lead <- numeric(157)
ROI["lead"] <- lead

ROI$lead <- "new_lead"

trait <- numeric(157)
ROI["trait"] <- trait

ROI$trait <- "SS"

#fix the first line

ROI[1, 8] <- "previous_&_new_lead"

#to and from column added

ROI$from <- ROI$Illumina_data.Coordinate_36

ROI$to <- ROI$Illumina_data.Coordinate_36

#order

ROI <- ROI[c(1, 2, 3, 4, 10, 11, 7, 8, 9)]

#from to change values to TSS and -/+ 250kb
#Manually, due to double-checked using Ensemble genome browser

TEST <- ROI

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

ROI <- TEST


#export

write.table(ROI, file = "ROI", sep = " ", row.names = FALSE, col.names = FALSE)
