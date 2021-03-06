setwd("/40/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/02-FASe-Achal")
load("FASe_Pheno_data.RData")
## FASe PHenotype from Aquilla
FASE_PHENO.Aquilla <- read.table("FASe_3894_Phenoscope.txt", header = T, stringsAsFactors = F)
dim(FASE_PHENO.Aquilla)
# 3892
View(FASE_PHENO.Aquilla)
as.data.frame(table(FASE_PHENO.Aquilla$PR))
# Var1 Freq
# 1  201907_USUHS_gDNA_SHERMAN   43
# 2                  Broad_WGS  142
# 3              Genentech_WES   92
# 4              Genentech_WGS   47
# 5                   LOAD_WES   33
# 6               Macrogen_WGS   20
# 7                 MAPT_A152T   21
# 8        MGI_FASeEOAD_201605  423
# 9            Otogenetics_WES  834
# 10          phs000572_201508  117
# 11          phs000572_201612  190
# 12          phs000572_201707   94
# 13          phs000572_201802 1540
# 14                   TGI_WES  298


# BLOOMFIELD <- read.table("/40/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/02-FASe-Achal/Bloomfield-WES-WGS-project.csv", header = T, sep = ",", stringsAsFactors = T)
# dim(BLOOMFIELD)
# # 9810
# 
# ## Check what other projects fall under FASe
# ADSPFAM <- read.delim("/40/AD/AD_Seq_Data/03.-phenotype/202105-ADSP_Umbrella_NG00067.v5/ADSPFamilyBasedPhenotypes_DS_2021.02.19_ALL.txt", header = T, sep = "\t", stringsAsFactors = F)
# dim(ADSPFAM)
# ADSPCACO <- read.delim("/40/AD/AD_Seq_Data/03.-phenotype/202105-ADSP_Umbrella_NG00067.v5/ADSPCaseControlPhenotypes_DS_2021.02.19_ALL.txt", header = T, sep = "\t", stringsAsFactors = F)
# dim(ADSPCACO)
# # 27547
# BLOOMFIELD$FASE <- BLOOMFIELD$Bloomfield.gvcf.id..SM...9810. %in% c(ADSPFAM$SUBJID, ADSPCACO$SUBJID)
# as.data.frame(table(BLOOMFIELD$Seq_project[BLOOMFIELD$FASE == TRUE])[table(BLOOMFIELD$Seq_project[BLOOMFIELD$FASE == TRUE]) > 0])

## Get project list for FASe and then extract samples from those projects
FASE_PR <- c(unique(FASE_PHENO.Aquilla$PR),"202103_ADSP_FUS-familial_WGS_UNNAMED", "202104_ADSP_site27-sync-n303_WGS_UNNAMED")
PHENO.Bloomfield <- read.delim("/100/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/01-Bloomfield-preQC/03-PLINK-QC-files/Bloomfield_8751_metadata_20211203.csv", header = T, sep = ",", stringsAsFactors = T)
dim(PHENO.Bloomfield)
# # 8751
# head(PHENO.Bloomfield)

sum(PHENO.Bloomfield$Seq_project %in% FASE_PR)
# 3931
FASE_PHENO.BLOOMFIELD <- PHENO.Bloomfield[PHENO.Bloomfield$Seq_project %in% FASE_PR,]

## GET FID and IID from FAM
FASE_FAM <- read.table("Bloomfield_9810-hwe-geno0.05-mind0.1-WXSm-missing-projects-include-good-IDS-V2_no_chr_FASE.fam", header = F)
FASE_FAM$V2 <- as.character(FASE_FAM$V2)
dim(FASE_PHENO.BLOOMFIELD)
# [1] 3931   26
sum(FASE_PHENO.BLOOMFIELD$Bloomfield.gvcf.id..SM...9810. %in% FASE_FAM$V2)
# 3931
FASE_PHENO.BLOOMFIELD <- cbind.data.frame(setNames(FASE_FAM[match(FASE_PHENO.BLOOMFIELD$Bloomfield.gvcf.id..SM...9810., FASE_FAM$V2),1:2], c("FID", "IID")), FASE_PHENO.BLOOMFIELD)
sum(FASE_PHENO.BLOOMFIELD$IID == FASE_PHENO.BLOOMFIELD$Bloomfield.gvcf.id..SM...9810.)
# 3931

## FIX # N/A characters
FASE_PHENO.BLOOMFIELD[FASE_PHENO.BLOOMFIELD=="#N/A"] <- NA
FASE_PHENO.BLOOMFIELD <- FASE_PHENO.BLOOMFIELD[-grep("^X$|X.", colnames(FASE_PHENO.BLOOMFIELD))]

write.table(FASE_PHENO.BLOOMFIELD[1:2], "/40/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/02-FASe-Achal/FASE_sample_list.txt", sep ="\t", col.names = T, quote = F, row.names = FALSE)

# #########
# ## IBD ##
# #########
# ## IBD
# library(data.table)
# IBD <- fread("Bloomfield_9810-hwe-geno0.05-mind0.1-WXSm-missing-projects-include-good-IDS-V2_no_chr_FASE-IBD.genome")
# 
# # ## Subset IBD
# # smallIBD <- IBD[IBD$Z0 < 0.5,]
# # 
# # ## Select samples of interest
# # if (exists("SELECTED_IBD")) {
# # rm(SELECT_IBD)
# # }
# # if (file.exists("selected_points.csv")) {
# #   file.remove(selected_points.csv)
# # }
# # ## Select interactively
# # source("https://raw.githubusercontent.com/achalneupane/rcodes/master/IBD_selection.r")
# # SELECT_IBD(smallIBD)
# # 
# # SELECTED_SAMPLES <- read.table("selected_points.csv", header = F)
# # IBD$color <- ifelse(IBD$key %in% SELECTED_SAMPLES$V1, "red", "black")
# # 
# # library(ggplot2)
# # ggplot(IBD, aes(x=Z0, y=Z1, color = as.factor(color)))+ geom_point() + ggtitle("Bloomfield-FASe - 3931") +
# #   scale_color_identity()
# # ggsave("Bloomfield-FASE_3931_IBD_selected.jpg", plot = last_plot(), device = NULL, scale = 1, width = 16, height = 9, dpi = 300, limitsize = TRUE)
# 
# 
# library(ggplot2)
# ggplot(IBD, aes(x=Z0, y=Z1, color = as.factor(color)))+ geom_point() + ggtitle("Bloomfield-FASe - 3931") +
#   scale_color_identity()
# ggsave("Bloomfield-FASE_3931_IBD.jpg", plot = last_plot(), device = NULL, scale = 1, width = 16, height = 9, dpi = 300, limitsize = TRUE)
# 
# 
# parent_offsping <- IBD[(IBD$Z0 < 0.25 & IBD$Z1 > 0.75), ]
# 
# parent_offsping$Relationship <- "parent-offspring"
# dim(parent_offsping)
# # 132
# 
# sibPairs <- IBD[(IBD$Z0 < 0.5 &
#                    IBD$Z0 > 0.10 &
#                    IBD$Z1 < 0.75 &
#                    IBD$Z1 > 0.25),] 
# sibPairs$Relationship <- "sib-pairs"
# dim(sibPairs)
# # 2191
# duplicates <- IBD[(IBD$Z0 < 0.25 &
#                      IBD$Z1 < 0.25), ]
# duplicates$Relationship <- "duplicates"
# dim(duplicates)
# # 0
# relatives.ALL <- rbind(parent_offsping, sibPairs, duplicates)
# # check how many times each sample is related
# relatives.ALL$nIID1 <- with(transform(relatives.ALL, n = 1),  ave(n, IID1, FUN = length))
# relatives.ALL$nIID2 <- with(transform(relatives.ALL, n = 1),  ave(n, IID2, FUN = length))
# write.table(relatives.ALL, "relatives_ALL.csv", sep ="\t", col.names = T, quote = F)
# 
# 
# 
# ## PCA
# # Merge HAPMAP ethnicity
# PCA <- read.table("Bloomfield_9810-hwe-geno0.05-mind0.1-WXSm-missing-projects-include-good-IDS-V2_no_chr_FASE-HAPMAP-MERGED3-for_PCA.eigenvec", header =T, stringsAsFactors=FALSE)
# dim(PCA)
# HAPMAP.ethnicty <- read.table("relationships_w_pops_121708.txt", header = T )
# head(HAPMAP.ethnicty)
# 
# 
# PCA$COHORT <- "FASe"
# PCA$COHORT <- HAPMAP.ethnicty$population[match(PCA$IID, HAPMAP.ethnicty$IID)]
# PCA <- PCA[c("FID", "IID", c(paste0("PC", 1:10), "COHORT"))]
# PCA$COHORT <- as.character(PCA$COHORT)
# PCA$COHORT[is.na(PCA$COHORT)] <- "FASe"
# write.table(PCA, "Bloomfield_FASe_3931-round1.txt", sep ="\t", col.names = T, quote = F)
# 
# 
# 
# ## Generate a new file that has IID, PC1,PC2, and a new column COHORT 
# library(ggplot2)
# p <- ggplot(PCA, aes(x=PC1, y=PC2, color=COHORT)) + geom_point() + xlab("PC1") + ylab("PC2") + ggtitle("Bloomfield-FASe-3931") +
#   scale_color_manual(values = c('green', 'black', 'red', "blue"))
# 
# 
# ## Select NHW samples only
# p <- p + annotate("rect", xmin=-0.006, xmax=0.02, ymin=-0.02, ymax= 0.03, 
#              fill=NA, colour="red") +
#   annotate("text", x=0.012, y=0.025, label="NHW", size=4, color = "red")
# p 
# ggsave("/100/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/02-FASe-Achal/Bloomfield-FASe-3931-PCs-COHORT.jpg", plot = p, device = NULL, scale = 1, width = 12, height = 8, dpi = 600, limitsize = TRUE)


# After running this PCA, I noticed that FID and IIDs are the same, so fixing that
# FIX MAP
FASE_PHENO.Aquilla$IID[(grepl("^MAP_", FASE_PHENO.Aquilla$FID))] <- paste0("MAP_", gsub("MAP_", "", as.character(FASE_PHENO.Aquilla$IID[(grepl("^MAP_", FASE_PHENO.Aquilla$FID))])))

## CLEAN .WGS and .WES strings in Bloomfield phenotype
FASE_PHENO.BLOOMFIELD$CLEANED_ID <- as.character(FASE_PHENO.BLOOMFIELD$Bloomfield.gvcf.id..SM...9810.)
FASE_PHENO.BLOOMFIELD$CLEANED_ID[(grepl(".WGS|.WES", FASE_PHENO.BLOOMFIELD$CLEANED_ID))] <- gsub("\\..*","",FASE_PHENO.BLOOMFIELD$CLEANED_ID[(grepl(".WGS|.WES", FASE_PHENO.BLOOMFIELD$CLEANED_ID))])

FASE_PHENO.BLOOMFIELD$CLEANED_ID[grepl ("SHERMAN", FASE_PHENO.BLOOMFIELD$Seq_project)] <- sapply(strsplit(as.character(FASE_PHENO.BLOOMFIELD$FullSample.FSM)[grepl ("SHERMAN", FASE_PHENO.BLOOMFIELD$Seq_project)], "\\^"), `[`, 1)



## MERGE Aquilla Phenotype with Bloomfield Phenotype
colnames(FASE_PHENO.Aquilla) <- paste0("AQ_", colnames(FASE_PHENO.Aquilla))
sum(FASE_PHENO.BLOOMFIELD$CLEANED_ID %in% FASE_PHENO.Aquilla$AQ_IID)
# 3648

FASE_PHENO.BLOOMFIELD <- cbind(FASE_PHENO.BLOOMFIELD, FASE_PHENO.Aquilla[match(FASE_PHENO.BLOOMFIELD$CLEANED_ID, FASE_PHENO.Aquilla$AQ_IID),])



FASE_PHENO.BLOOMFIELD$FOUND  <- ifelse(FASE_PHENO.BLOOMFIELD$CLEANED_ID %in% FASE_PHENO.BLOOMFIELD$AQ_IID, "YES", "NO")
table(FASE_PHENO.BLOOMFIELD$FOUND )
# NO  YES 
# 283 3648 

FASE_PHENO_NOT_FOUND <- FASE_PHENO.BLOOMFIELD[grepl("NO", FASE_PHENO.BLOOMFIELD$FOUND),]
## These are the samples not found in older release Aquilla
table(as.character(FASE_PHENO_NOT_FOUND$Seq_project))
# 201907_USUHS_gDNA_SHERMAN     202103_ADSP_FUS-familial_WGS_UNNAMED 202104_ADSP_site27-sync-n303_WGS_UNNAMED                                Broad_WGS 
# 2                                      164                                       86                                       31



## Merge old (Aquilla) and new (Bloomfield) Phenotype




# For any samples that were in previous release, I will update the FID from that release
FASE_PHENO.BLOOMFIELD$NEW_FID <- as.character(FASE_PHENO.BLOOMFIELD$AQ_FID)
FASE_PHENO.BLOOMFIELD$FID <- as.character(FASE_PHENO.BLOOMFIELD$FID)

# If NEW_FID is not empty, update FID from older release
FASE_PHENO.BLOOMFIELD$FID[!is.na(FASE_PHENO.BLOOMFIELD$NEW_FID)] <-  FASE_PHENO.BLOOMFIELD$NEW_FID[!is.na(FASE_PHENO.BLOOMFIELD$NEW_FID)]


# Check FID for BROAD_WGS
AURORA <- read.table("/40/AD/AD_Seq_Data/05.-Analyses/05-Aurora_201904/03b-Aurora-FASeME/201904-Aurora-phenotye.ped", header = T, stringsAsFactors = FALSE)
head(AURORA)
sum(FASE_PHENO_NOT_FOUND$FID [(FASE_PHENO_NOT_FOUND$CLEANED_ID %in% AURORA$IID)] %in% AURORA$FID)
# 31
# BROAD FID Looks ok!

FASE_PHENO_NOT_FOUND <- FASE_PHENO_NOT_FOUND[!grepl("Broad_WGS", FASE_PHENO_NOT_FOUND$Seq_project),]

## ??
## Get the missing phenotype from Aquilla
## STATUS
FASE_PHENO.BLOOMFIELD$STATUS_ADDED_FROM_AQUILLA[is.na(FASE_PHENO.BLOOMFIELD$STATUS..CC.) & !is.na(FASE_PHENO.BLOOMFIELD$AQ_STATUS)] <- "YES"
FASE_PHENO.BLOOMFIELD$STATUS..CC.[is.na(FASE_PHENO.BLOOMFIELD$STATUS..CC.) & !is.na(FASE_PHENO.BLOOMFIELD$AQ_STATUS)] <- FASE_PHENO.BLOOMFIELD$AQ_STATUS[is.na(FASE_PHENO.BLOOMFIELD$STATUS..CC.) & !is.na(FASE_PHENO.BLOOMFIELD$AQ_STATUS)]

## APOE
FASE_PHENO.BLOOMFIELD$APOE_ADDED_FROM_AQUILLA[is.na(FASE_PHENO.BLOOMFIELD$APOE) & !is.na(FASE_PHENO.BLOOMFIELD$AQ_APOE)] <- "YES"
FASE_PHENO.BLOOMFIELD$APOE[is.na(FASE_PHENO.BLOOMFIELD$APOE) & !is.na(FASE_PHENO.BLOOMFIELD$AQ_APOE)]  <- FASE_PHENO.BLOOMFIELD$AQ_APOE[is.na(FASE_PHENO.BLOOMFIELD$APOE) & !is.na(FASE_PHENO.BLOOMFIELD$AQ_APOE)] 

## AAO
FASE_PHENO.BLOOMFIELD$AAO <- as.numeric(as.character(FASE_PHENO.BLOOMFIELD$AAO))
FASE_PHENO.BLOOMFIELD$AAO_ADDED_FROM_AQUILLA[is.na(FASE_PHENO.BLOOMFIELD$AAO) & !is.na(FASE_PHENO.BLOOMFIELD$AQ_AAO)] <- "YES"
FASE_PHENO.BLOOMFIELD$AAO[is.na(FASE_PHENO.BLOOMFIELD$AAO) & !is.na(FASE_PHENO.BLOOMFIELD$AQ_AAO)]  <- FASE_PHENO.BLOOMFIELD$AQ_AAO[is.na(FASE_PHENO.BLOOMFIELD$AAO) & !is.na(FASE_PHENO.BLOOMFIELD$AQ_AAO)]

## ALA
FASE_PHENO.BLOOMFIELD$ALA <- as.numeric(as.character(FASE_PHENO.BLOOMFIELD$ALA))
FASE_PHENO.BLOOMFIELD$ALA_ADDED_FROM_AQUILLA[is.na(FASE_PHENO.BLOOMFIELD$ALA) & !is.na(FASE_PHENO.BLOOMFIELD$AQ_ALA)] <- "YES"
FASE_PHENO.BLOOMFIELD$ALA[is.na(FASE_PHENO.BLOOMFIELD$ALA) & !is.na(FASE_PHENO.BLOOMFIELD$AQ_ALA)] <- FASE_PHENO.BLOOMFIELD$AQ_ALA[is.na(FASE_PHENO.BLOOMFIELD$ALA) & !is.na(FASE_PHENO.BLOOMFIELD$AQ_ALA)]


FASE_PHENO.BLOOMFIELD$FOUND  <- ifelse(FASE_PHENO.BLOOMFIELD$CLEANED_ID %in% FASE_PHENO.BLOOMFIELD$AQ_IID, "YES", "NO")
FASE_PHENO.BLOOMFIELD$FOUND [grepl("Broad_WGS", FASE_PHENO.BLOOMFIELD$Seq_project)] <- "YES"
table(FASE_PHENO.BLOOMFIELD$FOUND)
# NO  YES 
# 252 3679 

# 1=CO, 2=CA.
## Get the Phenotype for 202103_ADSP_FUS-familial_WGS_UNNAMED and 202104_ADSP_site27-sync-n303_WGS_UNNAMED 
# Read ADSP FAM file
ADSPFambased <- read.delim("/40/AD/AD_Seq_Data/03.-phenotype/202105-ADSP_Umbrella_NG00067.v5/ADSPFamilyBasedPhenotypes_DS_2021.02.19_ALL.txt", header = T, stringsAsFactors = F, sep = "\t")
# ADSP_CACO <- read.delim("/40/AD/AD_Seq_Data/03.-phenotype/202105-ADSP_Umbrella_NG00067.v5/ADSPCaseControlPhenotypes_DS_2021.02.19_ALL.txt", header = T, stringsAsFactors = F, sep = "\t")
colnames(ADSPFambased) <- paste0("ADSPFAM_", colnames(ADSPFambased))
ADSPFambased$ADSPFAM_STATUS <- ADSPFambased$ADSPFAM_AD

ADSPFambased$ADSPFAM_STATUS[grepl("1|2|3|4", ADSPFambased$ADSPFAM_STATUS)] <- 2
ADSPFambased$ADSPFAM_STATUS[grepl("0", ADSPFambased$ADSPFAM_STATUS)] <- 1
ADSPFambased$ADSPFAM_STATUS[grepl("9|5", ADSPFambased$ADSPFAM_STATUS)] <- -9
sum(FASE_PHENO.BLOOMFIELD$CLEANED_ID %in% ADSPFambased$ADSPFAM_SUBJID)
# 250

FASE_PHENO.BLOOMFIELD <- cbind(FASE_PHENO.BLOOMFIELD, ADSPFambased[match(FASE_PHENO.BLOOMFIELD$CLEANED_ID, ADSPFambased$ADSPFAM_SUBJID),])

## First, fix FID from ADSPFAM
FASE_PHENO.BLOOMFIELD$FID[grepl ("202103_ADSP_FUS-familial_WGS_UNNAMED", FASE_PHENO.BLOOMFIELD$Seq_project) & !is.na(FASE_PHENO.BLOOMFIELD$ADSPFAM_FamID)] <- FASE_PHENO.BLOOMFIELD$ADSPFAM_FamID[grepl ("202103_ADSP_FUS-familial_WGS_UNNAMED", FASE_PHENO.BLOOMFIELD$Seq_project) & !is.na(FASE_PHENO.BLOOMFIELD$ADSPFAM_FamID)]

## AGE
FASE_PHENO.BLOOMFIELD$AGE <- as.numeric(gsub("\\+", "",as.character(FASE_PHENO.BLOOMFIELD$AGE)))
FASE_PHENO.BLOOMFIELD$AGE_ADDED_FROM_ADSPFAM[grepl ("202103_ADSP_FUS-familial_WGS_UNNAMED", FASE_PHENO.BLOOMFIELD$Seq_project) & !is.na(FASE_PHENO.BLOOMFIELD$ADSPFAM_FamID)] <- "YES"
FASE_PHENO.BLOOMFIELD$AGE[grepl ("202103_ADSP_FUS-familial_WGS_UNNAMED", FASE_PHENO.BLOOMFIELD$Seq_project) & !is.na(FASE_PHENO.BLOOMFIELD$ADSPFAM_FamID)] <- FASE_PHENO.BLOOMFIELD$ADSPFAM_Age[grepl ("202103_ADSP_FUS-familial_WGS_UNNAMED", FASE_PHENO.BLOOMFIELD$Seq_project) & !is.na(FASE_PHENO.BLOOMFIELD$ADSPFAM_FamID)]
# 
FASE_PHENO.BLOOMFIELD$STATUS..CC. <- as.numeric(as.character(FASE_PHENO.BLOOMFIELD$STATUS..CC.))
FASE_PHENO.BLOOMFIELD$STATUS_ADDED_FROM_ADSPFAM[grepl ("202103_ADSP_FUS-familial_WGS_UNNAMED", FASE_PHENO.BLOOMFIELD$Seq_project) & !is.na(FASE_PHENO.BLOOMFIELD$ADSPFAM_FamID)] <- "YES"
FASE_PHENO.BLOOMFIELD$STATUS..CC.[grepl ("202103_ADSP_FUS-familial_WGS_UNNAMED", FASE_PHENO.BLOOMFIELD$Seq_project) & !is.na(FASE_PHENO.BLOOMFIELD$ADSPFAM_FamID)] <- FASE_PHENO.BLOOMFIELD$ADSPFAM_STATUS[grepl ("202103_ADSP_FUS-familial_WGS_UNNAMED", FASE_PHENO.BLOOMFIELD$Seq_project) & !is.na(FASE_PHENO.BLOOMFIELD$ADSPFAM_FamID)]


# APOE
FASE_PHENO.BLOOMFIELD$APOE <- as.character(FASE_PHENO.BLOOMFIELD$APOE)
FASE_PHENO.BLOOMFIELD$APOE_ADDED_FROM_ADSPFAM[grepl ("202103_ADSP_FUS-familial_WGS_UNNAMED", FASE_PHENO.BLOOMFIELD$Seq_project)] <- "YES"
FASE_PHENO.BLOOMFIELD$APOE[grepl ("202103_ADSP_FUS-familial_WGS_UNNAMED", FASE_PHENO.BLOOMFIELD$Seq_project)] <- ADSPFambased$ADSPFAM_APOE[match(FASE_PHENO.BLOOMFIELD$IID[grepl ("202103_ADSP_FUS-familial_WGS_UNNAMED", FASE_PHENO.BLOOMFIELD$Seq_project)], ADSPFambased$ADSPFAM_SUBJID)]



table(FASE_PHENO.BLOOMFIELD$STATUS..CC.)
# -9    1    2     
# 271 1146 2194   


## Check Missing FAM
# STATUS
library(stringr)
## FIDs that needs to be updated
FID_MISSING <- as.data.frame(FASE_PHENO.BLOOMFIELD$IID[grepl ("202103_ADSP_FUS-familial_WGS_UNNAMED", FASE_PHENO.BLOOMFIELD$Seq_project)][is.na(FASE_PHENO.BLOOMFIELD$FID[grepl ("202103_ADSP_FUS-familial_WGS_UNNAMED", FASE_PHENO.BLOOMFIELD$Seq_project)])])
FID_MISSING$PR <- "202103_ADSP_FUS-familial_WGS_UNNAMED"
colnames(FID_MISSING) <- c("IID", "PR")

## NO IIDs match, so there is nothing to merge from ADSP FAM file (ADSPFamilyBasedPhenotypes_DS_2021.02.19_ALL.txt). ALL FIDs from 202104_ADSP_site27-sync-n303_WGS_UNNAMED needs to be updated
ADSPFambased$IID[match(FASE_PHENO.BLOOMFIELD$IID[grepl ("202104_ADSP_site27-sync-n303_WGS_UNNAMED", FASE_PHENO.BLOOMFIELD$Seq_project)], ADSPFambased$SUBJID)]

table(FASE_PHENO.BLOOMFIELD$STATUS..CC.)
# -9    1    2     
# 271 1146 2194   

## Remove FID_MISSING from the NHW list; Also remove 202104_ADSP_site27-sync-n303_WGS_UNNAMED samples because none of them are in ADSPFambased
FID_MISSING
FID_MISSING2 <- as.data.frame(FASE_PHENO.BLOOMFIELD$FID[grepl ("202104_ADSP_site27-sync-n303_WGS_UNNAMED", FASE_PHENO.BLOOMFIELD$Seq_project)]) 
FID_MISSING2$PR <- "202104_ADSP_site27-sync-n303_WGS_UNNAMED"
colnames(FID_MISSING2) <- c("IID", "PR")

FID_MISSING <- rbind.data.frame(FID_MISSING, FID_MISSING2)
FID_MISSING$FID <- FID_MISSING$IID

write.table(FID_MISSING, "Samples_with_FID_problem.txt", sep ="\t", col.names = T, quote = F)

FASE_PHENO.BLOOMFIELD <- FASE_PHENO.BLOOMFIELD [!FASE_PHENO.BLOOMFIELD$IID %in% FID_MISSING$IID,]

table(FASE_PHENO.BLOOMFIELD$STATUS..CC.)
# -9    1    2 
# 283 1174 2335

# Fix the age for 202103_ADSP_FUS-familial_WGS_UNNAMED and 202104_ADSP_site27-sync-n303_WGS_UNNAMED, 
# 1=CO, 2=CA
FASE_PHENO.BLOOMFIELD$ADDED_AAO_FROM_ADSPFAM[grepl ("202103_ADSP_FUS-familial_WGS_UNNAMED", FASE_PHENO.BLOOMFIELD$Seq_project) & FASE_PHENO.BLOOMFIELD$STATUS..CC. == 2 & !is.na(FASE_PHENO.BLOOMFIELD$AGE)] <- "YES"
FASE_PHENO.BLOOMFIELD$AAO[grepl ("202103_ADSP_FUS-familial_WGS_UNNAMED", FASE_PHENO.BLOOMFIELD$Seq_project) & FASE_PHENO.BLOOMFIELD$STATUS..CC. == 2 & !is.na(FASE_PHENO.BLOOMFIELD$AGE)] <- as.numeric(as.character(FASE_PHENO.BLOOMFIELD$AGE[grepl ("202103_ADSP_FUS-familial_WGS_UNNAMED", FASE_PHENO.BLOOMFIELD$Seq_project) & FASE_PHENO.BLOOMFIELD$STATUS..CC. == 2 & !is.na(FASE_PHENO.BLOOMFIELD$AGE)]))

FASE_PHENO.BLOOMFIELD$ADDED_ALA_FROM_ADSPFAM[grepl ("202103_ADSP_FUS-familial_WGS_UNNAMED", FASE_PHENO.BLOOMFIELD$Seq_project) & FASE_PHENO.BLOOMFIELD$STATUS..CC. == 1 & !is.na(FASE_PHENO.BLOOMFIELD$AGE)] <- "YES"
FASE_PHENO.BLOOMFIELD$ALA[grepl ("202103_ADSP_FUS-familial_WGS_UNNAMED", FASE_PHENO.BLOOMFIELD$Seq_project) & FASE_PHENO.BLOOMFIELD$STATUS..CC. == 1 & !is.na(FASE_PHENO.BLOOMFIELD$AGE)] <- as.numeric(as.character(FASE_PHENO.BLOOMFIELD$AGE[grepl ("202103_ADSP_FUS-familial_WGS_UNNAMED", FASE_PHENO.BLOOMFIELD$Seq_project) & FASE_PHENO.BLOOMFIELD$STATUS..CC. == 1 & !is.na(FASE_PHENO.BLOOMFIELD$AGE)]))

# # Also fix for 202104_ADSP_site27-sync-n303_WGS_UNNAMED
# FASE_PHENO.BLOOMFIELD$AAO[grepl ("202104_ADSP_site27-sync-n303_WGS_UNNAMED", FASE_PHENO.BLOOMFIELD$Seq_project) & FASE_PHENO.BLOOMFIELD$STATUS..CC. == 2 & !is.na(FASE_PHENO.BLOOMFIELD$AGE)] <- as.numeric(as.character(FASE_PHENO.BLOOMFIELD$AGE[grepl ("202104_ADSP_site27-sync-n303_WGS_UNNAMED", FASE_PHENO.BLOOMFIELD$Seq_project) & FASE_PHENO.BLOOMFIELD$STATUS..CC. == 2 & !is.na(FASE_PHENO.BLOOMFIELD$AGE)]))
# FASE_PHENO.BLOOMFIELD$ALA[grepl ("202104_ADSP_site27-sync-n303_WGS_UNNAMED", FASE_PHENO.BLOOMFIELD$Seq_project) & FASE_PHENO.BLOOMFIELD$STATUS..CC. == 1 & !is.na(FASE_PHENO.BLOOMFIELD$AGE)] <- as.numeric(as.character(FASE_PHENO.BLOOMFIELD$AGE[grepl ("202104_ADSP_site27-sync-n303_WGS_UNNAMED", FASE_PHENO.BLOOMFIELD$Seq_project) & FASE_PHENO.BLOOMFIELD$STATUS..CC. == 1 & !is.na(FASE_PHENO.BLOOMFIELD$AGE)]))


# Recode APOE4ANY
FASE_PHENO.BLOOMFIELD$APOE4ANY <- FASE_PHENO.BLOOMFIELD$APOE
FASE_PHENO.BLOOMFIELD$APOE4ANY[is.na(FASE_PHENO.BLOOMFIELD$APOE)] <- -9
FASE_PHENO.BLOOMFIELD$APOE4ANY[grepl("22", FASE_PHENO.BLOOMFIELD$APOE)] <- 0
FASE_PHENO.BLOOMFIELD$APOE4ANY[grepl("23|32", FASE_PHENO.BLOOMFIELD$APOE)] <- 0
FASE_PHENO.BLOOMFIELD$APOE4ANY[grepl("33", FASE_PHENO.BLOOMFIELD$APOE)] <- 0
FASE_PHENO.BLOOMFIELD$APOE4ANY[grepl("24|42", FASE_PHENO.BLOOMFIELD$APOE)] <- 1
FASE_PHENO.BLOOMFIELD$APOE4ANY[grepl("34|43", FASE_PHENO.BLOOMFIELD$APOE)] <- 1
FASE_PHENO.BLOOMFIELD$APOE4ANY[grepl("44", FASE_PHENO.BLOOMFIELD$APOE)] <- 1

table(FASE_PHENO.BLOOMFIELD$APOE4ANY)
# -9    0    1 
# 114 1688 2043 


BLOOMFIELD.AQUILLA.ADSPFAMbased <- FASE_PHENO.BLOOMFIELD

colnames(BLOOMFIELD.AQUILLA.ADSPFAMbased)[grepl("AQ_", colnames(BLOOMFIELD.AQUILLA.ADSPFAMbased))] <- paste0("AQUILLA_", gsub("AQ_", "", colnames(BLOOMFIELD.AQUILLA.ADSPFAMbased)[grepl("AQ_", colnames(BLOOMFIELD.AQUILLA.ADSPFAMbased))]))

write.table(BLOOMFIELD.AQUILLA.ADSPFAMbased, "/40/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/02-FASe-Achal/BLOOMFIELD.AQUILLA.ADSPFAMbased_phenotype_with_added_status.csv", sep =",", col.names = T, quote = F, row.names = FALSE)


# Select NHW and FASe samples, which we call NHW_SAMPLES here
# NHW_SAMPLES <- PCA[PCA$PC1 < 0.02 & PCA$PC1 > -0.006 & PCA$PC2 < 0.03 & PCA$PC2 > -0.02,]
# dim(NHW_SAMPLES)
# # 3952 13

NHW_SAMPLES <- PCA[grepl("CEU|FASe", PCA$COHORT),]


# Any samples with FID_MISSING
# https://www.notion.so/Comments-on-samples-with-missing-FID-c71d62de52a54c0696f6566da61ee70f
sum(NHW_SAMPLES$IID %in% FID_MISSING$IID)
# 137

# Any samples with FID_MISSING; but are related ? I will confirm with IBD as well
sum(FID_MISSING$IID %in% c(relatives.ALL$IID1, relatives.ALL$IID2))
# 4 of these are seemingly related, but only one of them is NHW (A-ADC-AD010911)
# [1] "A-ADC-AD010430" "A-ADC-AD001490" "A-ADC-AD010527" "A-ADC-AD010911"
RL1 <- relatives.ALL[relatives.ALL$IID1 %in% FID_MISSING$IID,]
RL2 <- relatives.ALL[relatives.ALL$IID2 %in% FID_MISSING$IID,]
tt <- rbind(RL1, RL2)
View(tt)


# > dim(NHW_SAMPLES)
# [1] 3952   14

## select samples within 3sd
# # NHW samples
# threeSD <- NHW_SAMPLES %>% 
#   filter(PC1 < 0.02 & PC1 > -0.006, PC2 < 0.03 & PC2 > -0.02) %>% 
#   filter(PC1 > mean(PC1)-(3*sd(PC1)) & PC1 < mean(PC1)+(3*sd(PC1)) & 
#            PC2 > mean(PC2)-(3*sd(PC2)) & PC2 < mean(PC2)+(3*sd(PC2))) %>% 
#   mutate(FID=FID,IID=IID) %>% dplyr::select(FID, IID)

## First add reported NHW samples from ADSP family to plot in PCA
## Add ADSP Family (ADSPFambased) ethnicity 
# recode Race and Ethnicity
ADSPFambased$Race_recoded[ADSPFambased$Race == 1] <- "AI"
ADSPFambased$Race_recoded[ADSPFambased$Race == 2] <- "Asian"
ADSPFambased$Race_recoded[ADSPFambased$Race == 3] <- "N_Hawaiian"
ADSPFambased$Race_recoded[ADSPFambased$Race == 4] <- "AA"
ADSPFambased$Race_recoded[ADSPFambased$Race == 5] <- "White"
ADSPFambased$Race_recoded[ADSPFambased$Race == 6] <- "-9"

ADSPFambased$Ethnicity_recoded[ADSPFambased$Ethnicity == 0] <- "NHL"
ADSPFambased$Ethnicity_recoded[ADSPFambased$Ethnicity == 1] <- "HL"

ADSPFambased$Race_Ethnicity <- paste(ADSPFambased$Race_recoded, ADSPFambased$Ethnicity_recoded, sep = ":")

## ADD this to PCA
table(ADSPFambased$Race_Ethnicity)
ADSPFambased_NHW <- ADSPFambased[grepl("White:NHL", ADSPFambased$Race_Ethnicity),]

sum(PCA$IID %in% ADSPFambased_NHW$SUBJID)
# 186
PCA$ADSPFambased_Ethnicity <- ""
PCA$ADSPFambased_Ethnicity <- ADSPFambased_NHW$SUBJID[match(PCA$IID, ADSPFambased_NHW$SUBJID)]


# Samples within SD cutoff in reference to CEU HAPMAP samples
NHW_SAMPLES_CEU <- NHW_SAMPLES[NHW_SAMPLES$COHORT == "CEU",]


## (HAPMAP)
## Select NHW samples only
p.sd <- ggplot(PCA, aes(x=PC1, y=PC2, color=COHORT)) + geom_point() + xlab("PC1") + ylab("PC2") + ggtitle("Bloomfield-FASe-3931") +
  scale_color_manual(values = c('green', 'black', 'red', "blue")) +
  annotate("text", x=0.008, y=0.025, label="NHW", size=4, color = "red")

library(tidyverse)


SDSelection.Table <- list()
SD.cutoff.all <- 3:5
for (i in 1:length(SD.cutoff.all)){
SD.cutoff <- SD.cutoff.all[i]  
PC1min <- (mean(NHW_SAMPLES_CEU$PC1) - (SD.cutoff*sd(NHW_SAMPLES_CEU$PC1)))
PC1max <- (mean(NHW_SAMPLES_CEU$PC1) + (SD.cutoff*sd(NHW_SAMPLES_CEU$PC1)))
PC2min <- (mean(NHW_SAMPLES_CEU$PC2) - (SD.cutoff*sd(NHW_SAMPLES_CEU$PC2)))
PC2max <- (mean(NHW_SAMPLES_CEU$PC2) + (SD.cutoff*sd(NHW_SAMPLES_CEU$PC2)))

SDSelection <- NHW_SAMPLES[NHW_SAMPLES$PC1 > PC1min & 
  NHW_SAMPLES$PC1 < PC1max &
  NHW_SAMPLES$PC2 > PC2min &
  NHW_SAMPLES$PC2 < PC2max,]

SDSelection.Table[[i]] <- as.vector(SDSelection$IID)

## Select NHW samples only
p.sd <- p.sd + annotate("rect", xmin=PC1min, xmax=PC1max, ymin=PC2min, ymax=PC2max, 
                  fill=NA, colour="red") +
  annotate("text", x=PC1max, y=PC2max, label=paste0("sd: ",SD.cutoff.all[i]), size=4, color = "black")
}


ggsave("/100/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/02-FASe-Achal/Bloomfield-FASe-3931-PCs-COHORT-different-sd-HapMap.jpg", plot = p.sd, device = NULL, scale = 1, width = 12, height = 8, dpi = 600, limitsize = TRUE)

#########################
## Get table of samples shared between three different SD cutoffs
as.character(lapply(SDSelection.Table, length))
# [1] "3826" "3923" "3939"

names(SDSelection.Table) <- c(paste("sd", SD.cutoff.all, sep = "."))

SD3 <- setNames(as.data.frame(SDSelection.Table$sd.3), "sd.3")
SD4 <- setNames(as.data.frame(SDSelection.Table$sd.4), "sd.4")
SD5 <- setNames(as.data.frame(SDSelection.Table$sd.5), "sd.5")


SD3_SD4.merged <- merge(x = SD4, y = transform(SD3, sd.3t = sd.3), by.x = "sd.4", by.y = "sd.3t", all = T)
SD3_SD4_SD5.merged <- merge(x = SD5, y = transform(SD3_SD4.merged, sd.4t = sd.4), by.x = "sd.5", by.y = "sd.4t", all = T)
# reorder cols
SD3_SD4_SD5.merged <- SD3_SD4_SD5.merged[, c("sd.3", "sd.4", "sd.5")]


# Exclude HapMap
sum(SD3_SD4_SD5.merged$sd.5 %in% FASE_PHENO.BLOOMFIELD$IID)
# 3690
SD3_SD4_SD5.merged <- SD3_SD4_SD5.merged[SD3_SD4_SD5.merged$sd.5 %in% FASE_PHENO.BLOOMFIELD$IID,]


# SD3_SD4_SD5.merged is the table of samples at sd 3-5; I will use this below; see heading "Finding samples within different SD and the related samples outside the SD margin"
#########################

## Add ethnicity from ADSP Family
df.ethnicity <- setNames(cbind.data.frame(PCA$PC1[!is.na(PCA$ADSPFambased_Ethnicity)], PCA$PC2[!is.na(PCA$ADSPFambased_Ethnicity)]), c("PC1", "PC2"))

p.sd.reportedNHW <- p.sd + geom_point(data = df.ethnicity, aes(col="Reported_NHW")) +
  scale_color_manual(values = c('green', 'black', 'red','yellow', "blue")) 
p.sd.reportedNHW

ggsave("/100/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/02-FASe-Achal/Bloomfield-FASe-3931-PCs-COHORT-different-sd-HapMap-with-ADSP-FAM-reportedNHW.jpg", plot = p.sd.reportedNHW, device = NULL, scale = 1, width = 12, height = 8, dpi = 600, limitsize = TRUE)

## (HAPMAP)
# #######################################################################################
# ## Finding samples within different SD and the related samples outside the SD margin ##
# #######################################################################################

# Nos of samples withing SD3, SD4 and SD5
SD3_SD4_SD5.merged %>% summarise_all(~ sum(!is.na(.)))
# sd.3 sd.4 sd.5
# 3578 3674 3690

# Nos of extra samples between SD3, SD4 and SD5
SD3_SD4_SD5.merged %>% summarise_all(~ sum(is.na(.)))
# sd.3 sd.4 sd.5
# 112   16    0

# create a list to store dataframe from different SDs
FASE_NHW_ALL_MEM_SD <- list()

for (i in 1:length(SD.cutoff.all)){
## Get phenot for NHW samples
FASE_PHENO_NHW <- FASE_PHENO.BLOOMFIELD[FASE_PHENO.BLOOMFIELD$IID %in% SD3_SD4_SD5.merged[,paste0("sd.", SD.cutoff.all[i])],]
dim(FASE_PHENO_NHW)
# [1] 3690   32



# Just to make sure I am not breaking any family apart, I will get all samples from the same family for three SDs I checked above
FASE_NHW_ALL_MEM <- NULL

FIDUNIQUE <- unique(FASE_PHENO_NHW$FID)
length(FIDUNIQUE)
# 2438
for (j in 1:length(FIDUNIQUE)){
print(paste0("Doing_SM_", j, "_from_SD:", SD.cutoff.all[i])) 
FASE_NHW_ALL_MEM.tmp <- FASE_PHENO.BLOOMFIELD [grepl(paste0("^", FIDUNIQUE[j], "$"), FASE_PHENO.BLOOMFIELD$FID),]
FASE_NHW_ALL_MEM <- rbind.data.frame(FASE_NHW_ALL_MEM, FASE_NHW_ALL_MEM.tmp)
print(nrow(FASE_NHW_ALL_MEM))
}

FASE_NHW_ALL_MEM_SD[[i]] <- FASE_NHW_ALL_MEM
names(FASE_NHW_ALL_MEM_SD)[[i]] <- paste0("sd.", SD.cutoff.all[i])

}

as.character(lapply(FASE_NHW_ALL_MEM_SD, nrow))
# "3628" "3684" "3692" # Nos of samples at SD3, SD4 and SD5 including any relatives outside

sum(FASE_NHW_ALL_MEM_SD$sd.5$IID %in% SD3_SD4_SD5.merged$sd.5)


## Additional related samples outside each SD group
extra.sd.3 <- FASE_NHW_ALL_MEM_SD$sd.3$IID[!FASE_NHW_ALL_MEM_SD$sd.3$IID %in% SD3_SD4_SD5.merged$sd.3]
extra.sd.4 <- FASE_NHW_ALL_MEM_SD$sd.4$IID[!FASE_NHW_ALL_MEM_SD$sd.4$IID %in% SD3_SD4_SD5.merged$sd.4]
extra.sd.5 <- FASE_NHW_ALL_MEM_SD$sd.5$IID[!FASE_NHW_ALL_MEM_SD$sd.5$IID %in% SD3_SD4_SD5.merged$sd.5]
sq <- seq(max(length(extra.sd.5), length(extra.sd.3)))
Relatives_found_outside_each_SD <- data.frame(extra.sd.3[sq], extra.sd.4[sq], extra.sd.5[sq])


# Add FID
SD3_SD4_SD5.merged$FID <- FASE_PHENO.BLOOMFIELD$FID[match(SD3_SD4_SD5.merged$sd.5, FASE_PHENO.BLOOMFIELD$IID)]

# Samples within each SD
write.table(SD3_SD4_SD5.merged, "/40/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/02-FASe-Achal/Samples_witihin_each_SD_Category.csv", sep =",", col.names = T, quote = F, row.names = FALSE)


## color layers to related samples outside SD
p.sd.reportedNHW <- p.sd + geom_point(data = df.ethnicity, aes(col="Reported_NHW")) +
  scale_color_manual(values = c('green', 'black', 'red','yellow', "blue")) 
p.sd.reportedNHW


Relatives_found_outside_each_SD[ , sapply(Relatives_found_outside_each_SD, function(x) all(is.na(x)) ) ] <- NULL
Relatives_found_outside_each_SD

library(tidyr)
library(dplyr)

COLOR_CODE <- Relatives_found_outside_each_SD %>% 
  pivot_longer(
    everything(),
    names_to = "SD",
    values_to = "Samples"
  ) %>% 
  na.omit() %>% 
  arrange(Samples)


COLOR_CODE <- cbind(COLOR_CODE,PCA[c(3:4)][match(COLOR_CODE$Samples, PCA$IID),])

## SD3
SD.cutoff <- 3  
PC1min <- (mean(NHW_SAMPLES_CEU$PC1) - (SD.cutoff*sd(NHW_SAMPLES_CEU$PC1)))
PC1max <- (mean(NHW_SAMPLES_CEU$PC1) + (SD.cutoff*sd(NHW_SAMPLES_CEU$PC1)))
PC2min <- (mean(NHW_SAMPLES_CEU$PC2) - (SD.cutoff*sd(NHW_SAMPLES_CEU$PC2)))
PC2max <- (mean(NHW_SAMPLES_CEU$PC2) + (SD.cutoff*sd(NHW_SAMPLES_CEU$PC2)))
ggplot(PCA, aes(x=PC1, y=PC2, color=COHORT)) + geom_point() + xlab("PC1") + ylab("PC2") + ggtitle("Bloomfield-FASe-3931") +
  scale_color_manual(values = c('green', 'black', 'red', "blue")) +
  annotate("text", x=0.008, y=0.025, label="NHW", size=4, color = "red") +
    annotate("rect", xmin=PC1min, xmax=PC1max, ymin=PC2min, ymax=PC2max, 
             fill=NA, colour="red") +
    geom_point(data = COLOR_CODE[3:4][COLOR_CODE$SD == paste0("extra.sd.", SD.cutoff, ".sq."),], aes(col="related")) +
    scale_color_manual(values = c('green', 'black', 'red','yellow', "blue")) +
    annotate("text", x=PC1max, y=PC2max, label=paste0("sd: ",SD.cutoff), size=4, color = "black")



## SD4
SD.cutoff <- 4
PC1min <- (mean(NHW_SAMPLES_CEU$PC1) - (SD.cutoff*sd(NHW_SAMPLES_CEU$PC1)))
PC1max <- (mean(NHW_SAMPLES_CEU$PC1) + (SD.cutoff*sd(NHW_SAMPLES_CEU$PC1)))
PC2min <- (mean(NHW_SAMPLES_CEU$PC2) - (SD.cutoff*sd(NHW_SAMPLES_CEU$PC2)))
PC2max <- (mean(NHW_SAMPLES_CEU$PC2) + (SD.cutoff*sd(NHW_SAMPLES_CEU$PC2)))
ggplot(PCA, aes(x=PC1, y=PC2, color=COHORT)) + geom_point() + xlab("PC1") + ylab("PC2") + ggtitle("Bloomfield-FASe-3931") +
  scale_color_manual(values = c('green', 'black', 'red', "blue")) +
  annotate("text", x=0.008, y=0.025, label="NHW", size=4, color = "red") +
  annotate("rect", xmin=PC1min, xmax=PC1max, ymin=PC2min, ymax=PC2max, 
           fill=NA, colour="red") +
  geom_point(data = COLOR_CODE[3:4][COLOR_CODE$SD == paste0("extra.sd.", SD.cutoff, ".sq."),], aes(col="related")) +
  scale_color_manual(values = c('green', 'black', 'red','yellow', "blue")) +
  annotate("text", x=PC1max, y=PC2max, label=paste0("sd: ",SD.cutoff), size=4, color = "black")


## SD5
SD.cutoff <- 5
PC1min <- (mean(NHW_SAMPLES_CEU$PC1) - (SD.cutoff*sd(NHW_SAMPLES_CEU$PC1)))
PC1max <- (mean(NHW_SAMPLES_CEU$PC1) + (SD.cutoff*sd(NHW_SAMPLES_CEU$PC1)))
PC2min <- (mean(NHW_SAMPLES_CEU$PC2) - (SD.cutoff*sd(NHW_SAMPLES_CEU$PC2)))
PC2max <- (mean(NHW_SAMPLES_CEU$PC2) + (SD.cutoff*sd(NHW_SAMPLES_CEU$PC2)))
ggplot(PCA, aes(x=PC1, y=PC2, color=COHORT)) + geom_point() + xlab("PC1") + ylab("PC2") + ggtitle("Bloomfield-FASe-3931") +
  scale_color_manual(values = c('green', 'black', 'red', "blue")) +
  annotate("text", x=0.008, y=0.025, label="NHW", size=4, color = "red") +
  annotate("rect", xmin=PC1min, xmax=PC1max, ymin=PC2min, ymax=PC2max, 
           fill=NA, colour="red") +
  geom_point(data = COLOR_CODE[3:4][COLOR_CODE$SD == paste0("extra.sd.", SD.cutoff, ".sq."),], aes(col="related")) +
  scale_color_manual(values = c('green', 'black', 'red','yellow', "blue")) +
  annotate("text", x=PC1max, y=PC2max, label=paste0("sd: ",SD.cutoff), size=4, color = "black")



## (1000 genome)
####################################
## PLOT PCA with 1000 genome data ##
####################################

setwd("/100/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/02-FASe-Achal/1000genome_PCA")
eigen<- read.table(file="PCA_FASe_1KG_merged.eigenvec",header=F) 
colnames(eigen)<-c("FID","IID", "C1","C2","C3","C4","C5","C6","C7","C8","C9","C10")

#Using super population as population code to plot 
race<- read.table(file="racefile_with_superpopcodes.txt",header=TRUE)
race$FID <- NULL 

datafile<- merge(eigen,race,by=c("IID"))
table(datafile$race)
# AFR  AMR  ASN  EUR FASe  SAN 
# 671  348  515  522 3931  492
##Save datafile 
write.table(datafile, file= "FASe_PCA.csv", col.names = TRUE, row.names = FALSE, sep = '\t', quote = F) 







# PCA plot with super populations
# AFR
# AMR
# ASN
# EUR
# SAN


p.sd.1KG.nosd <- ggplot(datafile, aes(x=C1, y=C2, color=race)) + geom_point() + xlab("PC1") + ylab("PC2") + ggtitle("Bloomfield-FASe-3931-With-1000genome") +
  scale_color_manual(values = c(AFR='green',
                                FASe='black',
                                EUR='red',
                                AMR="blue",
                                ASN="purple",
                                SAN="orange"))

NHW_SAMPLES_CEU_1KG <- datafile[grepl("EUR",datafile$race),]

SDSelection.Table <- list()
SD.cutoff.all <- 3:5
for (i in 1:length(SD.cutoff.all)){
  SD.cutoff <- SD.cutoff.all[i]  
  PC1min <- (mean(NHW_SAMPLES_CEU_1KG$C1) - (SD.cutoff*sd(NHW_SAMPLES_CEU_1KG$C1)))
  PC1max <- (mean(NHW_SAMPLES_CEU_1KG$C1) + (SD.cutoff*sd(NHW_SAMPLES_CEU_1KG$C1)))
  PC2min <- (mean(NHW_SAMPLES_CEU_1KG$C2) - (SD.cutoff*sd(NHW_SAMPLES_CEU_1KG$C2)))
  PC2max <- (mean(NHW_SAMPLES_CEU_1KG$C2) + (SD.cutoff*sd(NHW_SAMPLES_CEU_1KG$C2)))
  
  SDSelection <- datafile[datafile$C1 > PC1min & 
                            datafile$C1 < PC1max &
                            datafile$C2 > PC2min &
                            datafile$C2 < PC2max,]
  
  SDSelection.Table[[i]] <- as.vector(SDSelection$IID)
  
  ## Select NHW samples only
  p.sd.1KG.nosd <- p.sd.1KG.nosd + annotate("rect", xmin=PC1min, xmax=PC1max, ymin=PC2min, ymax=PC2max, 
                                  fill=NA, colour="red") +
    annotate("text", x=PC1max, y=PC2max, label=paste0("sd: ",SD.cutoff.all[i]), size=4, color = "black")
}

p.sd.1KG.nosd
ggsave("Bloomfield_9810-hwe-geno0.05-mind0.1-WXSm-missing-projects-include-good-IDS-V2_no_chr_FASE_UPDATED_FID_with_1000Genome.jpg", plot = p.sd.1KG.nosd, device = NULL, scale = 1, width = 12, height = 8, dpi = 600, limitsize = TRUE)


## (1000 genome)
#########################
## Get table of samples shared between three different SD cutoffs
as.character(lapply(SDSelection.Table, length))
# [1] "4239" "4320" "4367"

names(SDSelection.Table) <- c(paste("sd", SD.cutoff.all, sep = "."))

SD3 <- setNames(as.data.frame(SDSelection.Table$sd.3), "sd.3")
SD4 <- setNames(as.data.frame(SDSelection.Table$sd.4), "sd.4")
SD5 <- setNames(as.data.frame(SDSelection.Table$sd.5), "sd.5")


SD3_SD4.merged.1000G <- merge(x = SD4, y = transform(SD3, sd.3t = sd.3), by.x = "sd.4", by.y = "sd.3t", all = T)
SD3_SD4_SD5.merged.1000G <- merge(x = SD5, y = transform(SD3_SD4.merged.1000G, sd.4t = sd.4), by.x = "sd.5", by.y = "sd.4t", all = T)
# reorder cols
SD3_SD4_SD5.merged.1000G <- SD3_SD4_SD5.merged.1000G[, c("sd.3", "sd.4", "sd.5")]


# Exclude HapMap
sum(SD3_SD4_SD5.merged.1000G$sd.5 %in% FASE_PHENO.BLOOMFIELD$IID)
# 3692
SD3_SD4_SD5.merged.1000G <- SD3_SD4_SD5.merged.1000G[SD3_SD4_SD5.merged.1000G$sd.5 %in% FASE_PHENO.BLOOMFIELD$IID,]
# Two related samples excluded with HAPMAP analysis at 5 SD are included in this analysis

## (1000 genome)
# #######################################################################################
# ## Finding samples within different SD and the related samples outside the SD margin ##
# #######################################################################################

# Nos of samples withing SD3, SD4 and SD5; Object SD3_SD4_SD5.merged is for HAPMAP
SD3_SD4_SD5.merged.1000G %>% summarise_all(~ sum(!is.na(.)))
## HAPMAP
# sd.3 sd.4 sd.5
# 3578 3674 3690
## 1000 Genome
# sd.3 sd.4 sd.5
# 3614 3682 3692


# Nos of extra samples between SD3, SD4 and SD5
SD3_SD4_SD5.merged.1000G %>% summarise_all(~ sum(is.na(.)))
## HAPMAP
# sd.3 sd.4 sd.5
# 112   16    0
## 1000 Genome
# sd.3 sd.4 sd.5
# 78   10    0


# create a list to store dataframe from different SDs
FASE_NHW_ALL_MEM_SD.1000G <- list()

for (i in 1:length(SD.cutoff.all)){
  ## Get phenot for NHW samples
  FASE_PHENO_NHW <- FASE_PHENO.BLOOMFIELD[FASE_PHENO.BLOOMFIELD$IID %in% SD3_SD4_SD5.merged.1000G[,paste0("sd.", SD.cutoff.all[i])],]
  dim(FASE_PHENO_NHW)
  # [1] 3690   32
  
  
  
  # Just to make sure I am not breaking any family apart, I will get all samples from the same family for three SDs I checked above
  FASE_NHW_ALL_MEM <- NULL
  
  FIDUNIQUE <- unique(FASE_PHENO_NHW$FID)
  length(FIDUNIQUE)
  # 2438
  for (j in 1:length(FIDUNIQUE)){
    print(paste0("Doing_SM_", j, "_from_SD:", SD.cutoff.all[i])) 
    FASE_NHW_ALL_MEM.tmp <- FASE_PHENO.BLOOMFIELD [grepl(paste0("^", FIDUNIQUE[j], "$"), FASE_PHENO.BLOOMFIELD$FID),]
    FASE_NHW_ALL_MEM <- rbind.data.frame(FASE_NHW_ALL_MEM, FASE_NHW_ALL_MEM.tmp)
    print(nrow(FASE_NHW_ALL_MEM))
  }
  
  FASE_NHW_ALL_MEM_SD.1000G[[i]] <- FASE_NHW_ALL_MEM
  names(FASE_NHW_ALL_MEM_SD.1000G)[[i]] <- paste0("sd.", SD.cutoff.all[i])
  
}

as.character(lapply(FASE_NHW_ALL_MEM_SD.1000G, nrow))
## HAPMAP (all samples including relatives found outside SD3, SD4 and SD5 from HAPMAP analysis)
# "3628" "3684" "3692" # Nos of samples at SD3, SD4 and SD5 including any relatives outside
## 1000 Genome (all samples including relatives found outside SD3, SD4 and SD5 from 1000G analysis)
# "3646" "3688" "3699"


sum(FASE_NHW_ALL_MEM_SD.1000G$sd.5$IID %in% SD3_SD4_SD5.merged.1000G$sd.5)


## Additional related samples outside each SD group
extra.sd.3 <- FASE_NHW_ALL_MEM_SD.1000G$sd.3$IID[!FASE_NHW_ALL_MEM_SD.1000G$sd.3$IID %in% SD3_SD4_SD5.merged.1000G$sd.3]
extra.sd.4 <- FASE_NHW_ALL_MEM_SD.1000G$sd.4$IID[!FASE_NHW_ALL_MEM_SD.1000G$sd.4$IID %in% SD3_SD4_SD5.merged.1000G$sd.4]
extra.sd.5 <- FASE_NHW_ALL_MEM_SD.1000G$sd.5$IID[!FASE_NHW_ALL_MEM_SD.1000G$sd.5$IID %in% SD3_SD4_SD5.merged.1000G$sd.5]
sq <- seq(max(length(extra.sd.5), length(extra.sd.3)))

Relatives_found_outside_each_SD <- data.frame(extra.sd.3[sq], extra.sd.4[sq], extra.sd.5[sq])


# Add FID
SD3_SD4_SD5.merged.1000G$FID <- FASE_PHENO.BLOOMFIELD$FID[match(SD3_SD4_SD5.merged.1000G$sd.5, FASE_PHENO.BLOOMFIELD$IID)]

# Samples within each SD
write.table(SD3_SD4_SD5.merged.1000G, "/40/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/02-FASe-Achal/Samples_witihin_each_SD_Category_1000G.csv", sep =",", col.names = T, quote = F, row.names = FALSE)


## color layers to related samples outside SD
datafile$ADSPFambased_Ethnicity <- ADSPFambased_NHW$SUBJID[match(datafile$IID, ADSPFambased_NHW$SUBJID)]
df.ethnicity.C <- setNames(cbind.data.frame(datafile$C1[!is.na(datafile$ADSPFambased_Ethnicity)], datafile$C2[!is.na(datafile$ADSPFambased_Ethnicity)]), c("C1", "C2"))

p.sd.1KG.nosd.reportedNHW <- p.sd.1KG.nosd + geom_point(data = df.ethnicity.C, aes(col="Reported_NHW")) +
  scale_color_manual(
    values = c(
      AFR = 'green',
      FASe = 'black',
      EUR = 'red',
      AMR = 'blue',
      ASN = 'purple',
      SAN = 'orange',
      Reported_NHW = "yellow"
    )
  ) 
p.sd.1KG.nosd.reportedNHW
ggsave("Bloomfield_9810-hwe-geno0.05-mind0.1-WXSm-missing-projects-include-good-IDS-V2_no_chr_FASE_UPDATED_FID_with_1000Genome_with_reported_NHW.jpg", plot = p.sd.1KG.nosd.reportedNHW, device = NULL, scale = 1, width = 12, height = 8, dpi = 600, limitsize = TRUE)



Relatives_found_outside_each_SD[ , sapply(Relatives_found_outside_each_SD, function(x) all(is.na(x)) ) ] <- NULL
Relatives_found_outside_each_SD

library(tidyr)
library(dplyr)

COLOR_CODE <- Relatives_found_outside_each_SD %>% 
  pivot_longer(
    everything(),
    names_to = "SD",
    values_to = "Samples"
  ) %>% 
  na.omit() %>% 
  arrange(Samples)


COLOR_CODE <- cbind(COLOR_CODE,datafile[c(3:4)][match(COLOR_CODE$Samples, datafile$IID),])

## SD3
SD.cutoff <- 3  
PC1min <- (mean(NHW_SAMPLES_CEU_1KG$C1) - (SD.cutoff*sd(NHW_SAMPLES_CEU_1KG$C1)))
PC1max <- (mean(NHW_SAMPLES_CEU_1KG$C1) + (SD.cutoff*sd(NHW_SAMPLES_CEU_1KG$C1)))
PC2min <- (mean(NHW_SAMPLES_CEU_1KG$C2) - (SD.cutoff*sd(NHW_SAMPLES_CEU_1KG$C2)))
PC2max <- (mean(NHW_SAMPLES_CEU_1KG$C2) + (SD.cutoff*sd(NHW_SAMPLES_CEU_1KG$C2)))


ggplot(datafile, aes(x=C1, y=C2, color=race)) + geom_point() + xlab("PC1") + ylab("PC2") + ggtitle("Bloomfield-FASe-3931-1000G") +
  annotate("text", x=0.008, y=0.025, label="NHW", size=4, color = "red") +
  annotate("rect", xmin=PC1min, xmax=PC1max, ymin=PC2min, ymax=PC2max, 
           fill=NA, colour="red") +
  geom_point(data = COLOR_CODE[3:4][COLOR_CODE$SD == paste0("extra.sd.", SD.cutoff, ".sq."),], aes(col="related")) +
  scale_color_manual(values = c(AFR='green',
                                FASe='black',
                                EUR='red',
                                AMR="blue",
                                ASN="purple",
                                SAN="orange", 
                                related = "white")) +
  annotate("text", x=PC1max, y=PC2max, label=paste0("sd: ",SD.cutoff), size=4, color = "black")



## SD4
SD.cutoff <- 4  
PC1min <- (mean(NHW_SAMPLES_CEU_1KG$C1) - (SD.cutoff*sd(NHW_SAMPLES_CEU_1KG$C1)))
PC1max <- (mean(NHW_SAMPLES_CEU_1KG$C1) + (SD.cutoff*sd(NHW_SAMPLES_CEU_1KG$C1)))
PC2min <- (mean(NHW_SAMPLES_CEU_1KG$C2) - (SD.cutoff*sd(NHW_SAMPLES_CEU_1KG$C2)))
PC2max <- (mean(NHW_SAMPLES_CEU_1KG$C2) + (SD.cutoff*sd(NHW_SAMPLES_CEU_1KG$C2)))



ggplot(datafile, aes(x=C1, y=C2, color=race)) + geom_point() + xlab("PC1") + ylab("PC2") + ggtitle("Bloomfield-FASe-3931-1000G") +
  annotate("text", x=0.008, y=0.025, label="NHW", size=4, color = "red") +
  annotate("rect", xmin=PC1min, xmax=PC1max, ymin=PC2min, ymax=PC2max, 
           fill=NA, colour="red") +
  geom_point(data = COLOR_CODE[3:4][COLOR_CODE$SD == paste0("extra.sd.", SD.cutoff, ".sq."),], aes(col="related")) +
  scale_color_manual(values = c(AFR='green',
                                FASe='black',
                                EUR='red',
                                AMR="blue",
                                ASN="purple",
                                SAN="orange", 
                                related = "white")) +
  annotate("text", x=PC1max, y=PC2max, label=paste0("sd: ",SD.cutoff), size=4, color = "black")


## SD5
SD.cutoff <- 5
PC1min <- (mean(NHW_SAMPLES_CEU_1KG$C1) - (SD.cutoff*sd(NHW_SAMPLES_CEU_1KG$C1)))
PC1max <- (mean(NHW_SAMPLES_CEU_1KG$C1) + (SD.cutoff*sd(NHW_SAMPLES_CEU_1KG$C1)))
PC2min <- (mean(NHW_SAMPLES_CEU_1KG$C2) - (SD.cutoff*sd(NHW_SAMPLES_CEU_1KG$C2)))
PC2max <- (mean(NHW_SAMPLES_CEU_1KG$C2) + (SD.cutoff*sd(NHW_SAMPLES_CEU_1KG$C2)))



ggplot(datafile, aes(x=C1, y=C2, color=race)) + geom_point() + xlab("PC1") + ylab("PC2") + ggtitle("Bloomfield-FASe-3931-1000G") +
  annotate("text", x=0.008, y=0.025, label="NHW", size=4, color = "red") +
  annotate("rect", xmin=PC1min, xmax=PC1max, ymin=PC2min, ymax=PC2max, 
           fill=NA, colour="red") +
  geom_point(data = COLOR_CODE[3:4][COLOR_CODE$SD == paste0("extra.sd.", SD.cutoff, ".sq."),], aes(col="related")) +
  scale_color_manual(values = c(AFR='green',
                                FASe='black',
                                EUR='red',
                                AMR="blue",
                                ASN="purple",
                                SAN="orange", 
                                related = "white")) +
  annotate("text", x=PC1max, y=PC2max, label=paste0("sd: ",SD.cutoff), size=4, color = "black")



## Plot all the population in 1000 genome data
#Read population file
race<- read.table(file="racefile_with_popcodes.txt",header=T)
race$FID <- NULL 
datafile<- merge(eigen,race,by=c("IID"))
table(datafile$race)
# ACB  ASW  BEB  CDX  CEU  CHB  CHS  CLM  ESN FASe  FIN  GBR  GIH  GWD  IBS  ITU  JPT  KHV  LWK  MSL  MXL  PEL  PJL  PUR  STU 
# 97   61   86  100   99  106  105   95  100 3931  105  100  106  113  107  102  105   99  103   90   64   85   96  104  102 
# TSI  YRI 
# 111  107 

myCol = c("pink1", "violet", "mediumpurple1", "slateblue1", "purple", "purple3",
          "turquoise2", "skyblue", "steelblue", "black", "blue2", "navyblue",
          "orange", "tomato", "coral2", "palevioletred", "violetred", "red2",
          "springgreen2", "yellowgreen", "palegreen4",
          "wheat2", "tan", "tan2", "tan3",
          "grey70", "grey50")
length(myCol)
p <- ggplot(datafile, aes(x=C1, y=C2)) + geom_point(size=2, aes(colour =race)) + scale_color_manual(values= myCol)
p
ggsave("Bloomfield_9810-hwe-geno0.05-mind0.1-WXSm-missing-projects-include-good-IDS-V2_no_chr_FASE_UPDATED_FID_with_1000Genome_with_reported_NHW_all_population.jpg", plot = p, device = NULL, scale = 1, width = 12, height = 8, dpi = 600, limitsize = TRUE)

# Rest of the code is based on HapMap analysis
##################
## Demographics ##
##################

########################################################################################################
####################### Age-of-onset stratified by STATUS, ethnicity and sex ###########################
########################################################################################################
covars <- FASE_NHW_ALL_MEM_SD$sd.5

## Start demographic table
# FASE_PHENO.BLOOMFIELD$IID <- as.character(FASE_PHENO.BLOOMFIELD$IID)
table(table(covars$FID))
# 1    2    3    4    5    6    7    8    9   10   12   14   24 
# 2061   20  101  128   64   37   14    5    4    1    1    1    1 



colnames(covars)[colnames(covars) == "AAO"] <- "AGE_AT_ONSET"
colnames(covars)[colnames(covars) == "ALA"] <- "AGE_LAST_VISIT"
covars$ETHNICITY <- "NHW"
colnames(covars)[colnames(covars) == "STATUS..CC."] <- "STATUS"
colnames(covars)[colnames(covars) == "Pheno_SEX"] <- "SEX"

## Set binary values to APOE 
covars$APOE4ANY <- covars$APOE
sum(grepl("22|23|33", covars$APOE4ANY))
covars$APOE4ANY[grepl("22|23|33|32", covars$APOE4ANY)] <- 0
covars$APOE4ANY[grepl("24|34|44|42|43", covars$APOE4ANY)] <- 1

## Fix age for missing age
covars$AGE_AT_ONSET[covars$STATUS == 2 & is.na(covars$AGE_AT_ONSET) & !is.na(covars$AGE)] <- as.numeric(gsub("\\+", "", as.character(covars$AGE[covars$STATUS == 2 & is.na(covars$AGE_AT_ONSET) & !is.na(covars$AGE)])))
covars$AGE_LAST_VISIT [covars$STATUS == 1 & is.na(covars$AGE_LAST_VISIT) & !is.na(covars$AGE)] <- as.numeric(gsub("\\+", "", as.character(covars$AGE[covars$STATUS == 1 & is.na(covars$AGE_LAST_VISIT) & !is.na(covars$AGE)])))


covars$AGE_AT_ONSET <- as.numeric(as.character(covars$AGE_AT_ONSET))
covars$AGE_LAST_VISIT <- as.numeric(as.character(covars$AGE_LAST_VISIT))



Get_STATs_FASE <- function(covars){
  N.controls <- sum(covars$STATUS==1, na.rm = T)
  N.cases <- sum(covars$STATUS==2, na.rm = T)
  N.OTHER <- sum(covars$STATUS==-9, na.rm = T)
  TOTAL= nrow(covars)
  
  ## Split groups
  CONTROLS <- covars[covars$STATUS==1, ]
  CASES <- covars[covars$STATUS==2, ]
  UNKNOWNS <- covars[covars$STATUS==-9, ]
  
  ##############
  ## CONTROLS ##
  ##############
  ## Percent FEMALE
  FE.CO <- sum(CONTROLS$SEX == 2, na.rm = T)
  MA.CO <- sum(CONTROLS$SEX == 1, na.rm = T)
  PERC.FEMALE.CO <- (FE.CO/(MA.CO + FE.CO))*100
  ## Percent APOE4
  POS <- sum(CONTROLS[, c("APOE4ANY")] == 1, na.rm = T)
  NEG <- sum(CONTROLS[, c("APOE4ANY")] == 0, na.rm = T)
  PERC.APOE.CO <- (POS/ (POS + NEG)) *100
  # MISSING AGES CO
  N.CO.missing.age <- sum(is.na(CONTROLS$AGE_LAST_VISIT))
  Avg.AGE.SD.CO <- paste0(round(mean(CONTROLS$AGE_LAST_VISIT, na.rm = T), 2), " (", round(sd(CONTROLS$AGE_LAST_VISIT, na.rm = T), 2), ")")
  Age.range.CO <- paste(range(CONTROLS$AGE_LAST_VISIT[!CONTROLS$AGE_LAST_VISIT== -9], na.rm = T), collapse = "-")
  
  ###########
  ## CASES ##
  ###########
  ## Percent FEMALE
  FE.CA <- sum(CASES$SEX == 2, na.rm = T)
  MA.CA <- sum(CASES$SEX == 1, na.rm = T)
  PERC.FEMALE.CA <- (FE.CA/(MA.CA + FE.CA))*100
  ## Percent APOE4
  POS <- sum(CASES[, c("APOE4ANY")] == 1, na.rm = T)
  NEG <- sum(CASES[, c("APOE4ANY")] == 0, na.rm = T)
  PERC.APOE.CA <- round((POS/ (POS + NEG)) *100, 2)
  # MISSING AGES CA
  N.CA.missing.age <- sum(is.na(CASES$AGE_AT_ONSET))
  Avg.AGE.SD.CA <- paste0(round(mean(CASES$AGE_AT_ONSET, na.rm = T), 2), " (", round(sd(CASES$AGE_AT_ONSET, na.rm = T), 2), ")")
  Age.range.CA <- paste(range(CASES$AGE_AT_ONSET[!CASES$AGE_AT_ONSET== -9], na.rm = T), collapse = "-")
  
  
  ##########
  ## FASe ##
  ##########
  ## Percent FEMALE
  FE <- sum(covars$SEX == 2, na.rm = T)
  MA <- sum(covars$SEX == 1, na.rm = T)
  PERC.FEMALE <- (FE/(MA + FE))*100
  ## Percent APOE4
  POS <- sum(covars[, c("APOE4ANY")] == 1, na.rm = T)
  NEG <- sum(covars[, c("APOE4ANY")] == 0, na.rm = T)
  PERC.APOE <- (POS/ (POS + NEG)) *100
  
  
  STATS <- cbind.data.frame(N = rbind('CASES (2)' = N.cases, 'CONTROLS (1)' = N.controls, 'FASe (All)' = TOTAL), 
                            '%Female' = round(rbind(PERC.FEMALE.CA, PERC.FEMALE.CO, PERC.FEMALE), 2), 
                            '%APOE4' = round(rbind(PERC.APOE.CA, PERC.APOE.CO, PERC.APOE),2), 
                            'MISSING_AGE' = rbind(N.CA.missing.age, N.CO.missing.age, ""), 
                            'AvgerageAge (SD)' = rbind(Avg.AGE.SD.CA, Avg.AGE.SD.CO, ""), 
                            'AgeRange' = rbind(Age.range.CA, Age.range.CO, ""))
  
  return(STATS)
}



Demographic <- Get_STATs_FASE (covars)
View(Demographic)



# Check if family size of 1 are related
FAM_SIZE1 <- setNames(as.data.frame(table(covars$FID)), c("FID", "Freq"))
FAM_SIZE1 <- FAM_SIZE1[FAM_SIZE1$Freq  ==1,]
FAM_SIZE1$FID <- as.character(FAM_SIZE1$FID)
FAM_OF_SIZE1 <- covars[covars$FID %in% FAM_SIZE1$FID, c("FID", "IID")]
write.table(FAM_OF_SIZE1, "/100/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/02-FASe-Achal/samples_in_fam-size1.txt", col.names = T, row.names = F, quote = F, sep = "\t")



IBD0.3 <- IBD[IBD$PI_HAT > 0.3,]
IBD0.3$key <- paste(IBD0.3$FID1, IBD0.3$FID2, sep = ":")


# Check if any of the 1-member family is related
related.1.FAM <- {}
for (i in 1:nrow(FAM_SIZE1)){
  print(paste0("Doing iteration: ", i, " for FID: ",  FAM_SIZE1$FID[i]))
  patterns <- paste0(paste0("^", FAM_SIZE1$FID[i], ":"),"|", paste0(":", FAM_SIZE1$FID[i], "$"))
  related.1.FAM.tmp <- IBD0.3[grepl(patterns, IBD0.3$key),]
  related.1.FAM <- rbind.data.frame(related.1.FAM, related.1.FAM.tmp)
  # print(paste0("***", dim(related.1.FAM), "***"))
}

## remove duplicate rows
library(dplyr)
related.1.FAM <- distinct(related.1.FAM)

# Keep origial FID in covar
covars$FID_original <- covars$FID


# Since these samples are related, I will replace the FID of these samples as provided by Vicky
covars$FID [grepl("^A-ADC-AD010911$|^MAP_10071$", covars$IID)] <- "MR113"
covars$FID [grepl("^MAP_12367$|^MAP_84628$", covars$IID)] <- "FD106"
covars$FID [grepl("^MAP_61952$|^MAP_64869$", covars$IID)] <- "MR181"
covars$FID [grepl("^MAP_61401.WES.phs000572_201508$|^MAP_68161$", covars$IID)] <- "MR045"
covars$FID [grepl("^MAP_86353$|^MAP_86361$", covars$IID)] <- "FD238"
covars$FID [grepl("^MAP_86495.WES.MGI_FASeEOAD_201605$|^MAP_86496$", covars$IID)] <- "FD257"

table(table(covars$FID))
# 1    2    3    4    5    6    7    8    9   10   12   14   24 
# 2051   25  101  128   64   37   14    5    4    1    1    1    1 

## Find cases and controls of age < 65 yo
covars_controls <- as.data.table(covars[covars$STATUS == 1,])
covars_cases <- as.data.table(covars[covars$STATUS == 2,])

# Table of sample counts at different age cutoffs; check young cases and controls
covars_controls$AGE_LAST_VISIT <- as.numeric(as.character(covars_controls$AGE_LAST_VISIT))
covars_cases$AGE_AT_ONSET <- as.numeric(as.character(covars_cases$AGE_AT_ONSET))



ageGroup.cases <- cut(covars_cases$AGE_AT_ONSET,
                      breaks=c(0, 40, 50, 60, 65, 70, 75, 80, 85, 90, Inf),
                      include.lowest=TRUE,
                      labels=c("<=40", ">40-50", ">50-60", ">60-65", ">65-70",
                               ">70-75", ">75-80", ">80-85", ">85-90", ">90"))

table(ageGroup.cases)
# ageGroup.cases
# <=40 >40-50 >50-60 >60-65 >65-70 >70-75 >75-80 >80-85 >85-90    >90 
# 3     24    292    793    237    299    253    110     75     20 


ageGroup.controls <- cut(covars_controls$AGE_LAST_VISIT,
                        breaks=c(0, 40, 50, 60, 65, 70, 75, 80, 85, 90, Inf),
                        include.lowest=TRUE,
                        labels=c("<=40", ">40-50", ">50-60", ">60-65", ">65-70",
                                 ">70-75", ">75-80", ">80-85", ">85-90", ">90"))

table(ageGroup.controls)
# ageGroup.controls
# <40 40-50 50-60 60-65 65-70 70-75 75-80 80-85 85-90   >90 
# 2     5    15    12    30    54    52    72   792    33 

# sum(covars_controls$AGE_LAST_VISIT > 0 & covars_controls$AGE_LAST_VISIT <= 40 , na.rm = T)
# sum(covars_controls$AGE_LAST_VISIT > 85 & covars_controls$AGE_LAST_VISIT <= 90 , na.rm = T)

AGE.Groups <- setNames(cbind.data.frame(table(ageGroup.cases), table(ageGroup.controls))[c(1,2,4)], c("Age", "CA", "CO"))
AGE.Groups
# Age  CA  CO
# 1    <=40   3   2
# 2  >40-50  24   5
# 3  >50-60 292  15
# 4  >60-65 793  12
# 5  >65-70 237  30
# 6  >70-75 299  54
# 7  >75-80 253  52
# 8  >80-85 110  72
# 9  >85-90  75 792
# 10    >90  20  33


###################
## Controls < 65 ##
###################
sum(covars_controls$AGE_LAST_VISIT < 65 & covars_controls$AGE_LAST_VISIT > 0, na.rm = T)
# 33
covars_controls.65.yo <- covars_controls[covars_controls$AGE_LAST_VISIT < 65 & covars_controls$AGE_LAST_VISIT > 0,]
# See if these younger controls are related to each other
table(covars_controls.65.yo$FID)
# 27_104   4_715 LD0179F LD0223F LD0949F LD1223F LD1315F LD1329F LD1778F UM0147F UM0152F UM0187F UM0196F UM0304F UM0453F UP0005F 
# 1       1       1       2       3       4       3       1       4       1       1       3       3       1       1       3

# See if these younger controls have any other relatives in phenotype file
table(covars$FID)[names(table(covars$FID)) %in% names(table(covars_controls.65.yo$FID))]
# 27_104   4_715 LD0179F LD0223F LD0949F LD1223F LD1315F LD1329F LD1778F UM0147F UM0152F UM0187F UM0196F UM0304F UM0453F UP0005F 
# 6       6       2       4       6       5       6      14       4       6       3       4       8       3       5      12

# Remove all younger controls (< 65); we remove all younger controls
to.remove.younger.controls.IIDs <- names(table(covars_controls.65.yo$IID))
length(to.remove.younger.controls.IIDs)
# 33
################
## Cases < 65 ##
################
sum(covars_cases$AGE_AT_ONSET < 65 & covars_cases$AGE_AT_ONSET > 0, na.rm = T)
# 917
covars_cases.65.yo <- covars_cases[covars_cases$AGE_AT_ONSET < 65 & covars_cases$AGE_AT_ONSET > 0,]
# See if these younger cases are related
table(covars_cases.65.yo$FID) [table(covars_cases.65.yo$FID) > 1]
# 10R_R30 15_2041   17_24   25_70  27_165   27_97   4_596 5_26202 5_26408 NC0075F UM0453F 
# 2       2       2       2       2       2       2       2       2       2       2 


# See if these younger cases have any more relatives in phenotype file
table(covars$FID)[names(table(covars$FID)) %in% names(table(covars_cases.65.yo$FID))] [table(covars$FID)[names(table(covars$FID)) %in% names(table(covars_cases.65.yo$FID))] >1]
# 0_3041   0_3060   0_3115   0_3159   0_3169   0_3220   0_3296   0_3308  10J_111   10J_91 10R_R100  10R_R15  10R_R30  13_LD02  14_1956  15_2041  15_8052 
# 6        4        3        5        5        4        3        5        3        3        4        5        4        4        4        3        3 
# 17_24 19_L0019 19_L0022 19_L0026  19_L005    22_72    25_18    25_19     25_2    25_26    25_34    25_40    25_41    25_48    25_53     25_6    25_64 
# 4        3        4        5        4        3        7        4        6        3        3        4        5        4        4        5        3 
# 25_70     25_9   26_BCR   26_HTB   26_MKG   26_SIP   26_VKC   27_141   27_165   27_176     27_3     27_7    27_90    27_97    4_550    4_553    4_596 
# 6        4        6        4        4        5        4        5        3        4        4        4        5        5        3        5        4 
# 4_609    4_610    4_649    4_650    4_653    4_654    4_658  5_26044  5_26170  5_26202  5_26372  5_26408    7_108  8_62818  8_64024  8_64039  8_64053 
# 4        4        8        4        5        3        6        4        6        4        4        4        5        5        5        6        4 
# AD-26  AUX-003  NC0075F  UM0147F  UM0196F  UM0453F  UM0458F  UM0464F 
# 2        3        2        6        8        5        9        6

# to remove younger cases with no relatives in FASe
to.remove.younger.cases <- names(table(covars$FID)[names(table(covars$FID)) %in% names(table(covars_cases.65.yo$FID))] [table(covars$FID)[names(table(covars$FID)) %in% names(table(covars_cases.65.yo$FID))] == 1])
to.remove.younger.cases.IIDs <- covars_cases.65.yo$IID[covars_cases.65.yo$FID %in% to.remove.younger.cases]
length(to.remove.younger.cases.IIDs)
# 828


write.table(to.remove.younger.controls.IIDs, "younger_controls.txt", col.names = F, row.names = F)
write.table(to.remove.younger.cases.IIDs, "younger_cases_without_relatives.txt", col.names = F, row.names = F)
younger.samples.to.remove <- c(to.remove.younger.controls.IIDs, to.remove.younger.cases.IIDs)


FINAL.Covars <- covars[!covars$IID %in% younger.samples.to.remove,]
# Also remove sample from household size 24
# Household 24 ???> FID 203

sum(FINAL.Covars$FID %in% "203")
# 24
FINAL.Covars <- FINAL.Covars[!FINAL.Covars$FID %in% "203",]

dim(FINAL.Covars)
# 2807 

# Update samples FID
write.table(FINAL.Covars[c("FID_original", "IID", "FID", "IID")], "/100/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/02-FASe-Achal/sample_list_postQC_2807_update_FID.txt", col.names = F, row.names = F, quote = F, sep = "\t")

write.table(FINAL.Covars[1:2], "/100/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/02-FASe-Achal/sample_list_postQC_2807.txt", col.names = F, row.names = F, quote = F, sep = "\t")


# check if there are samples in FAM file
setwd("/100/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/02-FASe-Achal")
samples_in_FAM <- read.table("Bloomfield_9810-hwe-geno0.05-mind0.1-WXSm-missing-projects-include-good-IDS-V2_no_chr_FASE_UPDATED_FID_V2_post_QC2_V2.fam")

FINAL.Covars$KEY <- paste(FINAL.Covars$FID, FINAL.Covars$IID, sep = ":")
samples_in_FAM$KEY <- paste(samples_in_FAM$V1, samples_in_FAM$V2, sep = ":")
sum(FINAL.Covars$KEY %in% samples_in_FAM$KEY)
# 2807


View(FINAL.Covars)

Get_STATs_FASE(FINAL.Covars)

table(table(FINAL.Covars$FID))
# 1    2    3    4    5    6    7    8    9   10   13 
# 1226   27  101  126   66   32   14    4    5    1    1


# write.table(FINAL.Covars, "/100/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/02-FASe-Achal/Phenotype_2807.txt", col.names = T, row.names = F, quote = F, sep = "\t")




sum(FINAL.Covars$IID %in% FASE_PHENO.Aquilla$IID)
# 2588



## ADD PID and MID
FINAL.Covars.PID.MID <- cbind.data.frame(FINAL.Covars, FASE_PHENO.Aquilla[match(FINAL.Covars$CLEANED_ID, FASE_PHENO.Aquilla$IID), c("PID", "MID")])
FINAL.Covars.PID.MID <- FINAL.Covars.PID.MID[-c(3:4)]
FINAL.Covars.PID.MID$PID[grepl("^\\.$", FINAL.Covars.PID.MID$PID)] <- NA
FINAL.Covars.PID.MID$MID[grepl("^\\.$", FINAL.Covars.PID.MID$MID)] <- NA

FINAL.Covars.PID.MID$PID <- as.character(FINAL.Covars.PID.MID$PID)
FINAL.Covars.PID.MID$MID <- as.character(FINAL.Covars.PID.MID$MID)

FINAL.Covars.PID.MID$PID[is.na(FINAL.Covars.PID.MID$MID)] <- ADSPFambased$Father[match(FINAL.Covars.PID.MID$CLEANED_ID[is.na(FINAL.Covars.PID.MID$MID)], ADSPFambased$SUBJID)]
FINAL.Covars.PID.MID$MID[is.na(FINAL.Covars.PID.MID$MID)] <- ADSPFambased$Mother[match(FINAL.Covars.PID.MID$CLEANED_ID[is.na(FINAL.Covars.PID.MID$MID)], ADSPFambased$SUBJID)]
sum(is.na(FINAL.Covars.PID.MID$MID))
# 170
FINAL.Covars.PID.MID <- FINAL.Covars.PID.MID[-grep("^X$|X.", colnames(FINAL.Covars.PID.MID))]


## POST QC HapMap PCA
PCA <- read.table("/40/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/02-FASe-Achal/Bloomfield_9810-FASe-2807-hwe-geno0.05-mind0.1-WXSm-missing-projects-include-good-IDS-V2_no_chr_FASE_UPDATED_FID_post_QC3-HAPMAP-MERGED3-for_PCA.eigenvec", header =T, stringsAsFactors=FALSE)

library(ggplot2)
HAPMAP.ethnicty <- read.table("relationships_w_pops_121708.txt", header = T )
head(HAPMAP.ethnicty)

dim(PCA)
PCA$COHORT <- "FASe"
PCA$COHORT <- HAPMAP.ethnicty$population[match(PCA$IID, HAPMAP.ethnicty$IID)]
PCA <- PCA[c("FID", "IID", c(paste0("PC", 1:10), "COHORT"))]
PCA$COHORT <- as.character(PCA$COHORT)
PCA$COHORT[is.na(PCA$COHORT)] <- "FASe"
write.table(PCA, "/40/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/02-FASe-Achal/Bloomfield_9810-FASe-2807-Hapmap.txt", sep ="\t", col.names = T, quote = F)

##plotting:
p <- ggplot(PCA, aes(x=PC1, y=PC2, color=COHORT)) + geom_point() + xlab("PC1") + ylab("PC2") + ggtitle("Bloomfield_9810-FASe-2807") +
  scale_color_manual(values = c('green', 'black', 'red', "blue"))  
p
ggsave("/40/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/02-FASe-Achal/Bloomfield_9810-FASe-2807_with_HAPMAP.jpg", plot = p, device = NULL, scale = 1, width = 8, height = 4, dpi = 600, limitsize = TRUE)


setwd("/40/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/02-FASe-Achal")
load("FASe_Pheno_data.RData")

## Remove samples with mendelian variants
MENDELIAN <- read.delim("20211124_Bloomfield_pheno.csv", header = T, sep = ",")
dim(MENDELIAN)
MENDELIAN <- MENDELIAN[!is.na(MENDELIAN$Mendlian_mutation_name),]
MENDELIAN <- MENDELIAN[grepl("APP|PSEN", MENDELIAN$Mendlian_mutation_name, ignore.case = T),]

MENDELIAN$Bloomfield.gvcf.id..SM...9810. <- as.character(MENDELIAN$Bloomfield.gvcf.id..SM...9810.)
sum(FINAL.Covars.PID.MID$IID %in% MENDELIAN$Bloomfield.gvcf.id..SM...9810.)
FINAL.Covars.PID.MID <- FINAL.Covars.PID.MID[!FINAL.Covars.PID.MID$IID %in% MENDELIAN$Bloomfield.gvcf.id..SM...9810.,]

# write.table(FINAL.Covars.PID.MID, "/100/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/02-FASe-Achal/FBAT/Phenotype_with_PID_and_MID_2807.txt", col.names = T, row.names = F, quote = F, sep = "\t")
# write.table(FINAL.Covars.PID.MID[36:37], "/100/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/02-FASe-Achal/FBAT/recode_PID_and_MID_2807.txt", col.names = T, row.names = F, quote = F, sep = "\t")
write.table(FINAL.Covars.PID.MID[1:2], "/100/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/02-FASe-Achal/FBAT/subset_2804.txt", col.names = T, row.names = F, quote = F, sep = "\t")



# Re- PLot PCA
PCAs<- read.table("/100/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/02-FASe-Achal/FBAT/Bloomfield_9810-hwe-geno0.05-mind0.1-WXSm-missing-projects-include-good-IDS-V2_no_chr_FASE_UPDATED_FID_V2_post_QC2_V2_post_QC3-No-MCI-no-mendelian-PCAS.eigenvec", header=T)
##plotting:
library(ggplot2)
P <- ggplot(PCAs, aes(x=PC1, y=PC2)) + geom_point() + xlab("PC1") + ylab("PC2") + ggtitle("Bloomfield_9810-FASe-2804")
ggsave("/100/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/02-FASe-Achal/FBAT/Bloomfield_9810-FASe-2804-PCs-COHORT_After_keeping_NHW_and_without_HAPMAP.jpg", plot = P, device = NULL, scale = 1, width = 16, height = 9, dpi = 300, limitsize = TRUE)


## ADD PCA
setwd("/100/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/02-FASe-Achal/FBAT")
PCA.NHW.FASE <- read.table("Bloomfield_9810-hwe-geno0.05-mind0.1-WXSm-missing-projects-include-good-IDS-V2_no_chr_FASE_UPDATED_FID_V2_post_QC2_V2_post_QC3-No-MCI-no-mendelian-PCAS.eigenvec", header = T)
sum(FINAL.Covars.PID.MID$IID %in% PCA.NHW.FASE$IID)
# 2804
FINAL.Covars.PID.MID <- cbind.data.frame(FINAL.Covars.PID.MID, PCA.NHW.FASE[match(FINAL.Covars.PID.MID$IID, PCA.NHW.FASE$IID), -c(1:2)])

# paste(colnames(FINAL.Covars.PID.MID), collapse = ",")  
FINAL.Covars.PID.MID <- FINAL.Covars.PID.MID[c("FID","IID","PID","MID","SEX","STATUS","Genetic_Sex","APOE","AGE_AT_ONSET","AGE_LAST_VISIT","AGE","ADCO","Seq_project","COHORT","CLEANED_ID","FOUND","NEW_FID","APOE4ANY","KEY","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")]
  
colnames(FINAL.Covars.PID.MID)

FINAL.Covars.PID.MID.SELECTED.COLS <- FINAL.Covars.PID.MID [c("FID", "IID", "PID", "MID", "SEX", "STATUS", "APOE", "AGE_AT_ONSET", "AGE_LAST_VISIT", "Seq_project",
  "COHORT", "APOE4ANY", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")]

colnames(FINAL.Covars.PID.MID.SELECTED.COLS) <- c("pid", "id", "fid", "mid", "SEX", "STATUS", "APOE", "AAO", "ALA", "Seq_project",
                                                  "COHORT", "APOE4ANY", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")


## recode SEX 
FINAL.Covars.PID.MID.SELECTED.COLS$SEX[FINAL.Covars.PID.MID.SELECTED.COLS$SEX == -9| is.na(FINAL.Covars.PID.MID.SELECTED.COLS$SEX)] <- 0
FINAL.Covars.PID.MID.SELECTED.COLS$STATUS[FINAL.Covars.PID.MID.SELECTED.COLS$STATUS == -9| is.na(FINAL.Covars.PID.MID.SELECTED.COLS$STATUS)] <- 0

write.table(FINAL.Covars.PID.MID.SELECTED.COLS, "/100/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/02-FASe-Achal/FBAT/Phenotype_with_PID_and_MID_2804.txt", col.names = T, row.names = F, quote = F, sep = "\t")

write.table(FINAL.Covars.PID.MID.SELECTED.COLS[1:4], "/100/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/02-FASe-Achal/FBAT/recode_PID_and_MID_2804.txt", col.names = T, row.names = F, quote = F, sep = "\t")
write.table(FINAL.Covars.PID.MID.SELECTED.COLS[c(1,2,5)], "/100/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/02-FASe-Achal/FBAT/update_sex.txt", col.names = T, row.names = F, quote = F, sep = "\t")
write.table(FINAL.Covars.PID.MID.SELECTED.COLS[c(1,2,6)], "/100/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/02-FASe-Achal/FBAT/update_status.txt", col.names = T, row.names = F, quote = F, sep = "\t")

tt <- FINAL.Covars.PID.MID.SELECTED.COLS
colnames(tt)[1:2] <- c("FID", "IID")
tt$AGE_AT_ONSET <- as.numeric(as.character(tt$AAO))
tt$AGE_LAST_VISIT <- as.numeric(as.character(tt$ALA))
tt$SEX
tt$STATUS
Get_STATs_FASE(tt)

save.image("/100/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/02-FASe-Achal/FASe_Pheno_data.RData")


FAM_PLINK <- read.table("Bloomfield_9810-hwe-geno0.05-mind0.1-WXSm-missing-projects-include-good-IDS-V2_no_chr_FASE_UPDATED_FID_V2_post_QC2_V2_post_QC3-No-MCI-no-mendelian-geno-0.02-maxmaf-0.01_with_sex.fam", header = F, check.names = F)
head(FAM_PLINK)
FAM_PLINK$KEY <- paste(FAM_PLINK$V1, FAM_PLINK$V2, sep = ":")

FINAL.Covars.PID.MID.SELECTED.COLS$KEY <- paste(FINAL.Covars.PID.MID.SELECTED.COLS$pid, FINAL.Covars.PID.MID.SELECTED.COLS$id, sep = ":")
sum(FAM_PLINK$KEY %in% FINAL.Covars.PID.MID.SELECTED.COLS$KEY)

FAM_PLINK$PHENO_SEX <- FINAL.Covars.PID.MID.SELECTED.COLS$SEX[match(FAM_PLINK$V2, FINAL.Covars.PID.MID.SELECTED.COLS$id)]
sum(FAM_PLINK$PHENO_SEX == FAM_PLINK$V5)

read.table("Bloomfield_9810-hwe-geno0.05-mind0.1-WXSm-missing-projects-include-good-IDS-V2_no_chr_FASE_UPDATED_FID_V2_post_QC2_V2_post_QC3-No-MCI-no-mendelian-geno-0.02-maxmaf-0.01_sex.sexcheck")

PLINK.sex <- read.delim("../Bloomfield_9810-hwe-geno0.05-mind0.1-WXSm-missing-projects-include-good-IDS-V2_no_chr_FASE_UPDATED_FID_V2_post_QC2_V2_post_QC3_sex.sexcheck_formatted", header = T, sep = "\t")
dim(PLINK.sex)
SEX <- read.table("update_sex.txt", header = T, sep = "\t")
SEX$KEY <- paste(SEX$pid, SEX$id, sep = ":")

sum(SEX$id %in% PLINK.sex$IID )
PLINK.sex$PHENO_SEX <- FINAL.Covars$Pheno_SEX[match(PLINK.sex$IID, FINAL.Covars$IID )]

PLINK.sex$MATCH <- ifelse(PLINK.sex$PHENO_SEX == PLINK.sex$SNPSEX, "YES", "NO")
table(PLINK.sex$MATCH)

sum(SEX$KEY %in% FINAL.Covars$KEY )

SEX$PHENO_SEX <- FINAL.Covars$Pheno_SEX[match(SEX$KEY, FINAL.Covars$KEY )]

table(SEX$PHENO_SEX == SEX$SEX)
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################

PHENO <- read.delim("/100/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/01-Bloomfield-preQC/03-PLINK-QC-files/Bloomfield_8751_metadata_20211203.csv", header = T, sep = ",", stringsAsFactors = T)
dim(PHENO)
# # 8751
# head(PHENO)

sum(PHENO$Seq_project %in% FASE_PR)
# 3931
FASE_PHENO.BLOOMFIELD <- PHENO[PHENO$Seq_project %in% FASE_PR,]

table(FASE_PHENO.BLOOMFIELD$STATUS..CC.)
sum(FASE_PHENO.BLOOMFIELD$Bloomfield.gvcf.id..SM...9810. %in% FASE_PHENO.Aquilla$IID)
# 3452
colnames(FASE_PHENO.Aquilla) <- paste0("Aquilla_", colnames(FASE_PHENO.Aquilla))

sum(FINAL.Covars.PID.MID.SELECTED.COLS$id %in% FASE_PHENO.BLOOMFIELD$Bloomfield.gvcf.id..SM...9810.)

tt <- FINAL.Covars.PID.MID.SELECTED.COLS
colnames(tt) <- paste0("NEW_", colnames(tt))

tt <- cbind.data.frame(tt, FASE_PHENO.BLOOMFIELD[match(tt$NEW_id, FASE_PHENO.BLOOMFIELD$Bloomfield.gvcf.id..SM...9810.),])
sum(tt$NEW_STATUS == tt$STATUS..CC.)
# 2295
sum(tt$NEW_SEX == tt$Pheno_SEX)
# 2610
sum(tt$NEW_SEX == tt$Genetic_Sex)
# 2525


###########################################
#### Fixing the issues with Bloomfield ####
###########################################
library(data.table)
AQUILLA.VARS <- fread("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-a/Aquilla_7983_WXS_SNPS-INDELS-geno0.05-hwe-mind0.1-WXSmissingCLEAN-sexupdate-annot-snpeff-dbnsfp.vcf.tsv", header = T, sep = "\t" )
BLOOMFIELD.VARS <- fread("/40/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/Bloomfield_9810-hwe-geno0.05-mind0.1-WXSm-missing-projects-include-good-IDS-V2_VCF-snpeff-dbnsfp-ExAC.0.3.GRCh38-annot.vcf.tsv", header = T, sep = "\t")
BLOOMFIELD.REMOVED.VARS <- fread("/40/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/01-Bloomfield-preQC/03-PLINK-QC-files/Bloomfield_9810-hwe-geno0.05-mind0.1-missing-WXS.list", header = F)
BLOOMFIELD.REMOVED.VARS$KEY <- gsub("chr", "", BLOOMFIELD.REMOVED.VARS$V1)

AQUILLA.VARS$KEY <- paste(AQUILLA.VARS$`#CHROM`, AQUILLA.VARS$POS, AQUILLA.VARS$REF, AQUILLA.VARS$ALT, sep = ":")
BLOOMFIELD.VARS$KEY <- paste(BLOOMFIELD.VARS$`#CHROM`, BLOOMFIELD.VARS$POS, BLOOMFIELD.VARS$REF, BLOOMFIELD.VARS$ALT, sep = ":")

## Variants in Aquilla that are missing in Bloomfield
sum(!AQUILLA.VARS$KEY %in% BLOOMFIELD.VARS$KEY)
# 477781
VAR.in.AQUILLA.missing.in.BLOOMFIELD <- AQUILLA.VARS[!AQUILLA.VARS$KEY %in% BLOOMFIELD.VARS$KEY,]


VAR.in.AQUILLA.missing.in.BLOOMFIELD$consequence <- sapply(strsplit(VAR.in.AQUILLA.missing.in.BLOOMFIELD$INFO,"\\|"), `[`, 2)
VAR.in.AQUILLA.missing.in.BLOOMFIELD$gene <- sapply(strsplit(VAR.in.AQUILLA.missing.in.BLOOMFIELD$INFO,"\\|"), `[`, 3)
VAR.in.AQUILLA.missing.in.BLOOMFIELD$type <- sapply(strsplit(VAR.in.AQUILLA.missing.in.BLOOMFIELD$INFO,"\\|"), `[`, 4)
VAR.in.AQUILLA.missing.in.BLOOMFIELD$region <- sapply(strsplit(VAR.in.AQUILLA.missing.in.BLOOMFIELD$INFO,"\\|"), `[`, 5)
## No drop INFO field
VAR.in.AQUILLA.missing.in.BLOOMFIELD <- as.data.frame(VAR.in.AQUILLA.missing.in.BLOOMFIELD)


table(VAR.in.AQUILLA.missing.in.BLOOMFIELD$gene)

###########################################################
## Variants in Aquilla that were discarded in Bloomfield ##
###########################################################
sum(VAR.in.AQUILLA.missing.in.BLOOMFIELD$KEY %in% BLOOMFIELD.REMOVED.VARS$KEY)
# 169692

## APP|PSEN variants missing in BLOOMFIELD
sum(grepl("PSEN|APP", VAR.in.AQUILLA.missing.in.BLOOMFIELD$gene))
# 894

## IN CHR21
VAR.in.AQUILLA.missing.in.BLOOMFIELD.chr21 <- VAR.in.AQUILLA.missing.in.BLOOMFIELD[(grepl("^21:", VAR.in.AQUILLA.missing.in.BLOOMFIELD$KEY)),]
length(unique(VAR.in.AQUILLA.missing.in.BLOOMFIELD.chr21$gene))
# 271


#######################################################################################################
## Variants in Aquilla that were present in BLOOMFIELD but discarded during differential missingness ##
#######################################################################################################
VARIANTS.in.Aquilla.that.were.discarded.in.BLOOMFIELD <- VAR.in.AQUILLA.missing.in.BLOOMFIELD[VAR.in.AQUILLA.missing.in.BLOOMFIELD$KEY %in% BLOOMFIELD.REMOVED.VARS$KEY,]

## APP|PSEN Variants discarded in Bloomfield
sum(grepl("PSEN|APP", VARIANTS.in.Aquilla.that.were.discarded.in.BLOOMFIELD$gene))
## 443


#########################################################################################################################
## For problem variants table in https://www.notion.so/2022-01-27-chr21-chr22-missing-51880d946a80481f92c4731c57544562 ##
#########################################################################################################################
library(data.table)
AQUILLA.VARS <- fread("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-a/Aquilla_7983_WXS_SNPS-INDELS-geno0.05-hwe-mind0.1-WXSmissingCLEAN-sexupdate-annot-snpeff-dbnsfp.vcf.tsv", header = T, sep = "\t" )
BLOOMFIELD.VARS <- fread("/40/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/Bloomfield_9810-hwe-geno0.05-mind0.1-WXSm-missing-projects-include-good-IDS-V2_VCF-snpeff-dbnsfp-ExAC.0.3.GRCh38-annot.vcf.tsv", header = T, sep = "\t")

AQUILLA.VARS$KEY <- paste(AQUILLA.VARS$`#CHROM`, AQUILLA.VARS$POS, AQUILLA.VARS$REF, AQUILLA.VARS$ALT, sep = ":")
BLOOMFIELD.VARS$KEY <- paste(BLOOMFIELD.VARS$`#CHROM`, BLOOMFIELD.VARS$POS, BLOOMFIELD.VARS$REF, BLOOMFIELD.VARS$ALT, sep = ":")


AQUILLA.VARS$consequence <- sapply(strsplit(AQUILLA.VARS$INFO,"\\|"), `[`, 2)
AQUILLA.VARS$gene <- sapply(strsplit(AQUILLA.VARS$INFO,"\\|"), `[`, 3)
AQUILLA.VARS$type <- sapply(strsplit(AQUILLA.VARS$INFO,"\\|"), `[`, 4)
AQUILLA.VARS$region <- sapply(strsplit(AQUILLA.VARS$INFO,"\\|"), `[`, 5)


BLOOMFIELD.VARS$consequence <- sapply(strsplit(BLOOMFIELD.VARS$INFO,"\\|"), `[`, 2)
BLOOMFIELD.VARS$gene <- sapply(strsplit(BLOOMFIELD.VARS$INFO,"\\|"), `[`, 3)
BLOOMFIELD.VARS$type <- sapply(strsplit(BLOOMFIELD.VARS$INFO,"\\|"), `[`, 4)
BLOOMFIELD.VARS$region <- sapply(strsplit(BLOOMFIELD.VARS$INFO,"\\|"), `[`, 5)


dim(AQUILLA.VARS)
# 2143066       11

dim(BLOOMFIELD.VARS)
# 1916054       11

## Aquilla
CHR=19:22
for(i in 1:length(CHR)){
  CHR[i]  
  VARs <- AQUILLA.VARS[(grepl(paste0("^",CHR[i], ":"), AQUILLA.VARS$KEY)),]
  nos.VARs <- nrow(VARs)
  print(paste0("NOS of Variants in CHR", CHR[i], " : ", nos.VARs))
  nos.GENES <- length(unique(VARs$gene))
  print(paste0("NOS of genes in CHR", CHR[i], " : ", nos.GENES))    
}

# [1] "NOS of Variants in CHR19 : 138429"
# [1] "NOS of genes in CHR19 : 1783"
# [1] "NOS of Variants in CHR20 : 57112"
# [1] "NOS of genes in CHR20 : 711"
# [1] "NOS of Variants in CHR21 : 14975"
# [1] "NOS of genes in CHR21 : 271"
# [1] "NOS of Variants in CHR22 : 32170"
# [1] "NOS of genes in CHR22 : 556"

## BLOOMFIELD
CHR=19:22
for(i in 1:length(CHR)){
  CHR[i]  
  VARs <- BLOOMFIELD.VARS[(grepl(paste0("^",CHR[i], ":"), BLOOMFIELD.VARS$KEY)),]
  nos.VARs <- nrow(VARs)
  print(paste0("NOS of Variants in CHR", CHR[i], " : ", nos.VARs))
  nos.GENES <- length(unique(VARs$gene))
  print(paste0("NOS of genes in CHR", CHR[i], " : ", nos.GENES))    
}
# [1] "NOS of Variants in CHR19 : 116327"
# [1] "NOS of genes in CHR19 : 1848"
# [1] "NOS of Variants in CHR20 : 52469"
# [1] "NOS of genes in CHR20 : 702"
# [1] "NOS of Variants in CHR21 : 288"
# [1] "NOS of genes in CHR21 : 71"
# [1] "NOS of Variants in CHR22 : 822"
# [1] "NOS of genes in CHR22 : 166"

################################################
## Correlation between Aquilla and Bloomfield ##
################################################

CHR=20:22
TYPE=c("DP", "QD", "MQ", "AC", "AN", "FS")

##############
## WGS SNPs ##
##############
for(i in 1:length(CHR)){
  Bloomfield.table <- read.table(paste0("/40/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/01-Bloomfield-preQC/01-VQSR-ExAC-tsSNP99.6-tsINDEL95/QC_plots/Bloomfield_9810_IDs_WGS_4870_", "chr", CHR[i], ".vcf.gz-norm-ref-annotate.vcf-SNPs.vcf_recalibrated_SNPs.vcf.table"), header = TRUE, stringsAsFactors = F)
  head(Bloomfield.table)
  
  Bloomfield.table$KEY <- paste(Bloomfield.table$CHROM, Bloomfield.table$POS, Bloomfield.table$REF, Bloomfield.table$ALT, sep = ":")
  Aquilla.table <- read.table(paste0("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-a/01-Aquilla-preQC/01-VQSR-ExAC-tsSNP99.6-tsINDEL95/QC_plots/Aquilla_7983_IDs_WGS_3044_chr", CHR[i],".vcf.gz-norm-ref-annotate.vcf-SNPs.vcf_recalibrated_SNPs.vcf.table"), header = TRUE, stringsAsFactors = F)
  head(Aquilla.table)
  Aquilla.table$KEY <- paste(Aquilla.table$CHROM, Aquilla.table$POS, Aquilla.table$REF, Aquilla.table$ALT, sep = ":")
  
  colnames(Aquilla.table) <- paste0("AQ_", colnames(Aquilla.table))
  
  MERGED.AQ.BL <- merge(Bloomfield.table, Aquilla.table, by.x="KEY", by.y="AQ_KEY", all = T)
  
  for(j in 1:length(TYPE)){
    BLOOMFIELD <- as.numeric(as.character(MERGED.AQ.BL[,TYPE[j]]))
    AQUILLA <- as.numeric(as.character(MERGED.AQ.BL[, paste0("AQ_", TYPE[j])]))
    
    print(paste0("DOING: CHR", CHR[i], " and VAR: ", TYPE[j]))
    
    jpeg(paste0("WGS_SNPs_CHR", CHR[i],"_",TYPE[j],".jpg")) 
    plot(BLOOMFIELD, AQUILLA, pch = 19, col = "lightblue")
    # Regression line
    abline(lm(AQUILLA ~ BLOOMFIELD), col = "red", lwd = 3)
    # Pearson correlation
    text(paste("r=", round(cor(BLOOMFIELD, AQUILLA, use = "pairwise.complete.obs"), 2)), x = max(BLOOMFIELD, na.rm = T)/2, y = max(AQUILLA, na.rm = T)/2)
    dev.off()
    Sys.sleep(3)
  }
}



################
## WGS INDELS ##
################
for(i in 1:length(CHR)){
  Bloomfield.table <- read.table(paste0("/40/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/01-Bloomfield-preQC/01-VQSR-ExAC-tsSNP99.6-tsINDEL95/QC_plots/Bloomfield_9810_IDs_WGS_4870_chr", CHR[i], ".vcf.gz-norm-ref-annotate.vcf-INDELS.vcf_recalibrated_INDELs.vcf.table"), header = TRUE, stringsAsFactors = F)
  head(Bloomfield.table)
  
  Bloomfield.table$KEY <- paste(Bloomfield.table$CHROM, Bloomfield.table$POS, Bloomfield.table$REF, Bloomfield.table$ALT, sep = ":")
  Aquilla.table <- read.table(paste0("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-a/01-Aquilla-preQC/01-VQSR-ExAC-tsSNP99.6-tsINDEL95/QC_plots/Aquilla_7983_IDs_WGS_3044_chr", CHR[i],".vcf.gz-norm-ref-annotate.vcf-INDELS.vcf_recalibrated_INDELs.vcf.table"), header = TRUE, stringsAsFactors = F)
  head(Aquilla.table)
  Aquilla.table$KEY <- paste(Aquilla.table$CHROM, Aquilla.table$POS, Aquilla.table$REF, Aquilla.table$ALT, sep = ":")
  
  colnames(Aquilla.table) <- paste0("AQ_", colnames(Aquilla.table))
  
  MERGED.AQ.BL <- merge(Bloomfield.table, Aquilla.table, by.x="KEY", by.y="AQ_KEY", all = T)
  
  for(j in 1:length(TYPE)){
    BLOOMFIELD <- as.numeric(as.character(MERGED.AQ.BL[,TYPE[j]]))
    AQUILLA <- as.numeric(as.character(MERGED.AQ.BL[, paste0("AQ_", TYPE[j])]))
    
    print(paste0("DOING: CHR", CHR[i], " and VAR: ", TYPE[j]))
    
    jpeg(paste0("WGS_INDELs_CHR", CHR[i],"_",TYPE[j],".jpg")) 
    plot(BLOOMFIELD, AQUILLA, pch = 19, col = "lightblue")
    # Regression line
    abline(lm(AQUILLA ~ BLOOMFIELD), col = "red", lwd = 3)
    # Pearson correlation
    text(paste("r=", round(cor(BLOOMFIELD, AQUILLA, use = "pairwise.complete.obs"), 2)), x = max(BLOOMFIELD, na.rm = T)/2, y = max(AQUILLA, na.rm = T)/2)
    dev.off()
    Sys.sleep(3)
  }
}


##############
## WES SNPs ##
##############
for(i in 1:length(CHR)){
  Bloomfield.table <- read.table(paste0("/40/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/01-Bloomfield-preQC/01-VQSR-ExAC-tsSNP99.6-tsINDEL95/QC_plots/Bloomfield_9810_IDs_WES_4940_chr", CHR[i], ".vcf.gz-norm-ref-annotate.vcf-SNPs.vcf_recalibrated_SNPs.vcf.table"), header = TRUE, stringsAsFactors = F)
  head(Bloomfield.table)
  
  Bloomfield.table$KEY <- paste(Bloomfield.table$CHROM, Bloomfield.table$POS, Bloomfield.table$REF, Bloomfield.table$ALT, sep = ":")
  Aquilla.table <- read.table(paste0("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-a/01-Aquilla-preQC/01-VQSR-ExAC-tsSNP99.6-tsINDEL95/QC_plots/Aquilla_7983_IDs_WES_4941_chr", CHR[i],".vcf.gz-norm-ref-annotate.vcf-SNPs.vcf_recalibrated_SNPs.vcf.table"), header = TRUE, stringsAsFactors = F)
  head(Aquilla.table)
  Aquilla.table$KEY <- paste(Aquilla.table$CHROM, Aquilla.table$POS, Aquilla.table$REF, Aquilla.table$ALT, sep = ":")
  
  colnames(Aquilla.table) <- paste0("AQ_", colnames(Aquilla.table))
  
  MERGED.AQ.BL <- merge(Bloomfield.table, Aquilla.table, by.x="KEY", by.y="AQ_KEY", all = T)
  
  for(j in 1:length(TYPE)){
    BLOOMFIELD <- as.numeric(as.character(MERGED.AQ.BL[,TYPE[j]]))
    AQUILLA <- as.numeric(as.character(MERGED.AQ.BL[, paste0("AQ_", TYPE[j])]))
    
    print(paste0("DOING: CHR", CHR[i], " and VAR: ", TYPE[j]))
    
    jpeg(paste0("WES_SNPs_CHR", CHR[i],"_",TYPE[j],".jpg")) 
    plot(BLOOMFIELD, AQUILLA, pch = 19, col = "lightblue")
    # Regression line
    abline(lm(AQUILLA ~ BLOOMFIELD), col = "red", lwd = 3)
    # Pearson correlation
    text(paste("r=", round(cor(BLOOMFIELD, AQUILLA, use = "pairwise.complete.obs"), 2)), x = max(BLOOMFIELD, na.rm = T)/2, y = max(AQUILLA, na.rm = T)/2)
    dev.off()
    Sys.sleep(3)
  }
}


################
## WES INDELS ##
################
for(i in 1:length(CHR)){
  Bloomfield.table <- read.table(paste0("/40/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/01-Bloomfield-preQC/01-VQSR-ExAC-tsSNP99.6-tsINDEL95/QC_plots/Bloomfield_9810_IDs_WES_4940_chr", CHR[i], ".vcf.gz-norm-ref-annotate.vcf-INDELS.vcf_recalibrated_INDELs.vcf.table"), header = TRUE, stringsAsFactors = F)
  head(Bloomfield.table)
  
  Bloomfield.table$KEY <- paste(Bloomfield.table$CHROM, Bloomfield.table$POS, Bloomfield.table$REF, Bloomfield.table$ALT, sep = ":")
  Aquilla.table <- read.table(paste0("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-a/01-Aquilla-preQC/01-VQSR-ExAC-tsSNP99.6-tsINDEL95/QC_plots/Aquilla_7983_IDs_WES_4941_chr", CHR[i],".vcf.gz-norm-ref-annotate.vcf-INDELS.vcf_recalibrated_INDELs.vcf.table"), header = TRUE, stringsAsFactors = F)
  head(Aquilla.table)
  Aquilla.table$KEY <- paste(Aquilla.table$CHROM, Aquilla.table$POS, Aquilla.table$REF, Aquilla.table$ALT, sep = ":")
  
  colnames(Aquilla.table) <- paste0("AQ_", colnames(Aquilla.table))
  
  MERGED.AQ.BL <- merge(Bloomfield.table, Aquilla.table, by.x="KEY", by.y="AQ_KEY", all = T)
  
  for(j in 1:length(TYPE)){
    BLOOMFIELD <- as.numeric(as.character(MERGED.AQ.BL[,TYPE[j]]))
    AQUILLA <- as.numeric(as.character(MERGED.AQ.BL[, paste0("AQ_", TYPE[j])]))
    
    print(paste0("DOING: CHR", CHR[i], " and VAR: ", TYPE[j]))
    
    jpeg(paste0("WES_INDELs_CHR", CHR[i],"_",TYPE[j],".jpg")) 
    plot(BLOOMFIELD, AQUILLA, pch = 19, col = "lightblue")
    # Regression line
    abline(lm(AQUILLA ~ BLOOMFIELD), col = "red", lwd = 3)
    # Pearson correlation
    text(paste("r=", round(cor(BLOOMFIELD, AQUILLA, use = "pairwise.complete.obs"), 2)), x = max(BLOOMFIELD, na.rm = T)/2, y = max(AQUILLA, na.rm = T)/2)
    dev.off()
    Sys.sleep(3)
  }
}


#####################################
## Correlation between WGS and WES ##
#####################################

# CHR=c(1:22, "X")
CHR=1
TYPE=c("QD", "MQ")

#####################
## SNPs Bloomfield ##
#####################
for(i in 1:length(CHR)){
  Bloomfield.table.WGS <- read.table(paste0("/40/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/01-Bloomfield-preQC/01-VQSR-ExAC-tsSNP99.6-tsINDEL95/QC_plots/Bloomfield_9810_IDs_WGS_4870_", "chr", CHR[i], ".vcf.gz-norm-ref-annotate.vcf-SNPs.vcf_recalibrated_SNPs.vcf.table"), header = TRUE, stringsAsFactors = F)
  head(Bloomfield.table.WGS)
  
  Bloomfield.table.WGS$KEY <- paste(Bloomfield.table.WGS$CHROM, Bloomfield.table.WGS$POS, Bloomfield.table.WGS$REF, Bloomfield.table.WGS$ALT, sep = ":")
  
  Bloomfield.table.WES <- read.table(paste0("/40/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/01-Bloomfield-preQC/01-VQSR-ExAC-tsSNP99.6-tsINDEL95/QC_plots/Bloomfield_9810_IDs_WES_4940_chr", CHR[i], ".vcf.gz-norm-ref-annotate.vcf-SNPs.vcf_recalibrated_SNPs.vcf.table"), header = TRUE, stringsAsFactors = F)
  head(Bloomfield.table.WES)
  
  Bloomfield.table.WES$KEY <- paste(Bloomfield.table.WES$CHROM, Bloomfield.table.WES$POS, Bloomfield.table.WES$REF, Bloomfield.table.WES$ALT, sep = ":")
  colnames(Bloomfield.table.WES) <- paste0("WES_", colnames(Bloomfield.table.WES))
  
  MERGED.AQ.BL <- merge(Bloomfield.table.WGS, Bloomfield.table.WES, by.x="KEY", by.y="WES_KEY", all = T)
  
  for(j in 1:length(TYPE)){
    WGS <- as.numeric(as.character(MERGED.AQ.BL[,TYPE[j]]))
    WES <- as.numeric(as.character(MERGED.AQ.BL[, paste0("WES_", TYPE[j])]))
    
    print(paste0("DOING: CHR", CHR[i], " and VAR: ", TYPE[j]))
    
    jpeg(paste0("BLOOMFIELD_WGS_Vs_WES_SNPs_CHR", CHR[i],"_",TYPE[j],".jpg")) 
    plot(WGS, WES, pch = 19, col = "lightblue")
    # Regression line
    abline(lm(WES ~ WGS), col = "red", lwd = 3)
    # Pearson correlation
    text(paste("r=", round(cor(WGS, WES, use = "pairwise.complete.obs"), 2)), x = max(WGS, na.rm = T)/2, y = max(WES, na.rm = T)/2)
    dev.off()
    Sys.sleep(3)
  }
}



#######################
## INDELs Bloomfield ##
#######################
for(i in 1:length(CHR)){
  Bloomfield.table.WGS <- read.table(paste0("/40/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/01-Bloomfield-preQC/01-VQSR-ExAC-tsSNP99.6-tsINDEL95/QC_plots/Bloomfield_9810_IDs_WGS_4870_chr", CHR[i], ".vcf.gz-norm-ref-annotate.vcf-INDELS.vcf_recalibrated_INDELs.vcf.table"), header = TRUE, stringsAsFactors = F)
  head(Bloomfield.table.WGS)
  
  Bloomfield.table.WGS$KEY <- paste(Bloomfield.table.WGS$CHROM, Bloomfield.table.WGS$POS, Bloomfield.table.WGS$REF, Bloomfield.table.WGS$ALT, sep = ":")
  
  Bloomfield.table.WES <- read.table(paste0("/40/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/01-Bloomfield-preQC/01-VQSR-ExAC-tsSNP99.6-tsINDEL95/QC_plots/Bloomfield_9810_IDs_WES_4940_chr", CHR[i], ".vcf.gz-norm-ref-annotate.vcf-INDELS.vcf_recalibrated_INDELs.vcf.table"), header = TRUE, stringsAsFactors = F)
  head(Bloomfield.table.WES)
  
  Bloomfield.table.WES$KEY <- paste(Bloomfield.table.WES$CHROM, Bloomfield.table.WES$POS, Bloomfield.table.WES$REF, Bloomfield.table.WES$ALT, sep = ":")
  colnames(Bloomfield.table.WES) <- paste0("WES_", colnames(Bloomfield.table.WES))
  
  MERGED.AQ.BL <- merge(Bloomfield.table.WGS, Bloomfield.table.WES, by.x="KEY", by.y="WES_KEY", all = T)
  
  for(j in 1:length(TYPE)){
    WGS <- as.numeric(as.character(MERGED.AQ.BL[,TYPE[j]]))
    WES <- as.numeric(as.character(MERGED.AQ.BL[, paste0("WES_", TYPE[j])]))
    
    print(paste0("DOING: CHR", CHR[i], " and VAR: ", TYPE[j]))
    
    jpeg(paste0("BLOOMFIELD_WGS_Vs_WES_INDELs_CHR", CHR[i],"_",TYPE[j],".jpg")) 
    plot(WGS, WES, pch = 19, col = "lightblue")
    # Regression line
    abline(lm(WES ~ WGS), col = "red", lwd = 3)
    # Pearson correlation
    text(paste("r=", round(cor(WGS, WES, use = "pairwise.complete.obs"), 2)), x = max(WGS, na.rm = T)/2, y = max(WES, na.rm = T)/2)
    dev.off()
    Sys.sleep(3)
  }
}


##################
## SNPs Aquilla ##
##################
for(i in 1:length(CHR)){
  Aquilla.table.WGS <- read.table(paste0("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-a/01-Aquilla-preQC/01-VQSR-ExAC-tsSNP99.6-tsINDEL95/QC_plots/Aquilla_7983_IDs_WGS_3044_chr", CHR[i],".vcf.gz-norm-ref-annotate.vcf-SNPs.vcf_recalibrated_SNPs.vcf.table"), header = TRUE, stringsAsFactors = F)
  head(Aquilla.table.WGS)
  
  Aquilla.table.WGS$KEY <- paste(Aquilla.table.WGS$CHROM, Aquilla.table.WGS$POS, Aquilla.table.WGS$REF, Aquilla.table.WGS$ALT, sep = ":")
  
  Aquilla.table.WES <- read.table(paste0("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-a/01-Aquilla-preQC/01-VQSR-ExAC-tsSNP99.6-tsINDEL95/QC_plots/Aquilla_7983_IDs_WES_4941_chr", CHR[i],".vcf.gz-norm-ref-annotate.vcf-SNPs.vcf_recalibrated_SNPs.vcf.table"), header = TRUE, stringsAsFactors = F)
  head(Aquilla.table.WES)
  
  Aquilla.table.WES$KEY <- paste(Aquilla.table.WES$CHROM, Aquilla.table.WES$POS, Aquilla.table.WES$REF, Aquilla.table.WES$ALT, sep = ":")
  colnames(Aquilla.table.WES) <- paste0("WES_", colnames(Aquilla.table.WES))
  
  MERGED.AQ.BL <- merge(Aquilla.table.WGS, Aquilla.table.WES, by.x="KEY", by.y="WES_KEY", all = T)
  
  for(j in 1:length(TYPE)){
    WGS <- as.numeric(as.character(MERGED.AQ.BL[,TYPE[j]]))
    WES <- as.numeric(as.character(MERGED.AQ.BL[, paste0("WES_", TYPE[j])]))
    
    print(paste0("DOING: CHR", CHR[i], " and VAR: ", TYPE[j]))
    
    jpeg(paste0("AQUILLA_WGS_Vs_WES_SNPs_CHR", CHR[i],"_",TYPE[j],".jpg")) 
    plot(WGS, WES, pch = 19, col = "lightblue")
    # Regression line
    abline(lm(WES ~ WGS), col = "red", lwd = 3)
    # Pearson correlation
    text(paste("r=", round(cor(WGS, WES, use = "pairwise.complete.obs"), 2)), x = max(WGS, na.rm = T)/2, y = max(WES, na.rm = T)/2)
    dev.off()
    Sys.sleep(3)
  }
}


####################
## INDELs Aquilla ##
####################
for(i in 1:length(CHR)){
  Aquilla.table.WGS <- read.table(paste0("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-a/01-Aquilla-preQC/01-VQSR-ExAC-tsSNP99.6-tsINDEL95/QC_plots/Aquilla_7983_IDs_WGS_3044_chr", CHR[i],".vcf.gz-norm-ref-annotate.vcf-INDELS.vcf_recalibrated_INDELs.vcf.table"), header = TRUE, stringsAsFactors = F)
  head(Aquilla.table.WGS)
  
  Aquilla.table.WGS$KEY <- paste(Aquilla.table.WGS$CHROM, Aquilla.table.WGS$POS, Aquilla.table.WGS$REF, Aquilla.table.WGS$ALT, sep = ":")
  
  Aquilla.table.WES <- read.table(paste0("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-a/01-Aquilla-preQC/01-VQSR-ExAC-tsSNP99.6-tsINDEL95/QC_plots/Aquilla_7983_IDs_WES_4941_chr", CHR[i],".vcf.gz-norm-ref-annotate.vcf-INDELS.vcf_recalibrated_INDELs.vcf.table"), header = TRUE, stringsAsFactors = F)
  head(Aquilla.table.WES)
  
  Aquilla.table.WES$KEY <- paste(Aquilla.table.WES$CHROM, Aquilla.table.WES$POS, Aquilla.table.WES$REF, Aquilla.table.WES$ALT, sep = ":")
  colnames(Aquilla.table.WES) <- paste0("WES_", colnames(Aquilla.table.WES))
  
  MERGED.AQ.BL <- merge(Aquilla.table.WGS, Aquilla.table.WES, by.x="KEY", by.y="WES_KEY", all = T)
  
  for(j in 1:length(TYPE)){
    WGS <- as.numeric(as.character(MERGED.AQ.BL[,TYPE[j]]))
    WES <- as.numeric(as.character(MERGED.AQ.BL[, paste0("WES_", TYPE[j])]))
    
    print(paste0("DOING: CHR", CHR[i], " and VAR: ", TYPE[j]))
    
    jpeg(paste0("AQUILLA_WGS_Vs_WES_INDELs_CHR", CHR[i],"_",TYPE[j],".jpg")) 
    plot(WGS, WES, pch = 19, col = "lightblue")
    # Regression line
    abline(lm(WES ~ WGS), col = "red", lwd = 3)
    # Pearson correlation
    text(paste("r=", round(cor(WGS, WES, use = "pairwise.complete.obs"), 2)), x = max(WGS, na.rm = T)/2, y = max(WES, na.rm = T)/2)
    dev.off()
    Sys.sleep(3)
  }
}




# ###################
# ## Piared t-test ##
# ###################
# 
# before <-c(200.1, 190.9, 192.7, 213, 241.4, 196.9, 172.2, 185.5, 205.2, 193.7)
# # Weight of the mice after treatment
# after <-c(392.9, 393.2, 345.1, 393, 434, 427.9, 422, 383.9, 392.3, 352.2)
# # Create a data frame
# my_data <- data.frame( 
#   group = rep(c("before", "after"), each = 10),
#   weight = c(before,  after)
# )
# 
# 
# before <- subset(my_data,  group == "before", weight,
#                  drop = TRUE)
# # subset weight data after treatment
# after <- subset(my_data,  group == "after", weight,
#                 drop = TRUE)
# # Plot paired data
# # install.packages("ggplot2")
# library(ggplot2)
# library(PairedData)
# pd <- paired(before, after)
# plot(pd, type = "profile") + theme_bw()
# 
# pd <- paired(MERGED.AQ.BL$QD, MERGED.AQ.BL$WES_QD)
# colnames(pd) <- c("QD", "WES_QD")
# 
# res <- t.test(pd$QD, pd$WES_QD, paired = TRUE)
# res
# #####################
# ## SNPs Bloomfield ##
# #####################
# CHR=22
# # TYPE=c("QD", "MQ")
# TYPE=QD
# i=1
# j=1
# # for(i in 1:length(CHR)){
#   Bloomfield.table.WGS <- read.table(paste0("/40/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/01-Bloomfield-preQC/01-VQSR-ExAC-tsSNP99.6-tsINDEL95/QC_plots/Bloomfield_9810_IDs_WGS_4870_", "chr", CHR[i], ".vcf.gz-norm-ref-annotate.vcf-SNPs.vcf_recalibrated_SNPs.vcf.table"), header = TRUE, stringsAsFactors = F)
#   ## Get PASS
#   Bloomfield.table.WGS.FAILED<- Bloomfield.table.WGS[!grepl("^PASS$", Bloomfield.table.WGS$FILTER),]
#   Bloomfield.table.WGS<- Bloomfield.table.WGS[grepl("^PASS$", Bloomfield.table.WGS$FILTER),]
#   head(Bloomfield.table.WGS)
#   
#   Bloomfield.table.WGS$KEY <- paste(Bloomfield.table.WGS$CHROM, Bloomfield.table.WGS$POS, Bloomfield.table.WGS$REF, Bloomfield.table.WGS$ALT, sep = ":")
#   
#   Bloomfield.table.WES <- read.table(paste0("/40/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/01-Bloomfield-preQC/01-VQSR-ExAC-tsSNP99.6-tsINDEL95/QC_plots/Bloomfield_9810_IDs_WES_4940_chr", CHR[i], ".vcf.gz-norm-ref-annotate.vcf-SNPs.vcf_recalibrated_SNPs.vcf.table"), header = TRUE, stringsAsFactors = F)
#   ## Get PASS
#   Bloomfield.table.WES.FAILED<- Bloomfield.table.WES[!grepl("^PASS$", Bloomfield.table.WES$FILTER),]
#   Bloomfield.table.WES<- Bloomfield.table.WES[grepl("^PASS$", Bloomfield.table.WES$FILTER),]
#   head(Bloomfield.table.WES)
#   head(Bloomfield.table.WES)
#   
#   Bloomfield.table.WES$KEY <- paste(Bloomfield.table.WES$CHROM, Bloomfield.table.WES$POS, Bloomfield.table.WES$REF, Bloomfield.table.WES$ALT, sep = ":")
#   colnames(Bloomfield.table.WES) <- paste0("WES_", colnames(Bloomfield.table.WES))
#   
#   MERGED.AQ.BL <- merge(Bloomfield.table.WGS, Bloomfield.table.WES, by.x="KEY", by.y="WES_KEY", all = T)
#   
# 
#   tt <- cbind.data.frame(MERGED.AQ.BL$KEY, MERGED.AQ.BL$QD, MERGED.AQ.BL$WES_QD)
#   
#   colnames(tt) <- c("KEY", "WGS_QD", "WES_QD")
#   head(tt)
#   # KEY WGS_QD WES_QD
#   # 1 chr21:10413594:T:C   6.03     NA
#   # 2 chr21:10413598:A:G   8.35   8.35
#   # 3 chr21:10413604:G:C   0.94   0.94
#   
#   library(dplyr)
#   tt <- tt %>% 
#     filter(if_all(everything(), ~ !is.na(.x)))
#   tt <- tt[1:100,]
#   
#   t.test(tt$WES_QD, tt$WGS_QD, paired = T)
# 
#   
#   my_data <- data.frame( 
#     TYPE = rep(c("WES_QD", "WGS_QD"), each = nrow(tt)),
#     QD = c(tt$WGS_QD,  tt$WES_QD)
#   )
#   
#   res <- t.test(QD ~ TYPE, data = my_data, paired = TRUE)
#   
#   
# #   
# #   
# #   for(j in 1:length(TYPE)){
# #     WGS <- as.numeric(as.character(MERGED.AQ.BL[,TYPE[j]]))
# #     WES <- as.numeric(as.character(MERGED.AQ.BL[, paste0("WES_", TYPE[j])]))
# #     
# #     print(paste0("DOING: CHR", CHR[i], " and VAR: ", TYPE[j]))
# #     
# #     jpeg(paste0("BLOOMFIELD_WGS_Vs_WES_SNPs_CHR", CHR[i],"_",TYPE[j],".jpg")) 
# #     plot(WGS, WES, pch = 19, col = "lightblue")
# #     # Regression line
# #     abline(lm(WES ~ WGS), col = "red", lwd = 3)
# #     # Pearson correlation
# #     text(paste("r=", round(cor(WGS, WES, use = "pairwise.complete.obs"), 2)), x = max(WGS, na.rm = T)/2, y = max(WES, na.rm = T)/2)
# #     dev.off()
# #     Sys.sleep(3)
# #   }
# # }
# # 
  
#################################  
## Plot Missingness Bloomfield ##
#################################
CHR=19:22
# TYPE=c("QD", "MQ")
TYPE="AN"
# i=1
# j=1
for(i in 1:length(CHR)){
Bloomfield.table.WGS <- read.table(paste0("/40/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/01-Bloomfield-preQC/01-VQSR-ExAC-tsSNP99.6-tsINDEL95/QC_plots/Bloomfield_9810_IDs_WGS_4870_", "chr", CHR[i], ".vcf.gz-norm-ref-annotate.vcf-SNPs.vcf_recalibrated_SNPs.vcf.table"), header = TRUE, stringsAsFactors = F)
## Get PASS
Bloomfield.table.WGS.FAILED<- Bloomfield.table.WGS[!grepl("^PASS$", Bloomfield.table.WGS$FILTER),]
Bloomfield.table.WGS<- Bloomfield.table.WGS[grepl("^PASS$", Bloomfield.table.WGS$FILTER),]
head(Bloomfield.table.WGS)

Bloomfield.table.WGS$KEY <- paste(Bloomfield.table.WGS$CHROM, Bloomfield.table.WGS$POS, Bloomfield.table.WGS$REF, Bloomfield.table.WGS$ALT, sep = ":")
  
Bloomfield.table.WES <- read.table(paste0("/40/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/01-Bloomfield-preQC/01-VQSR-ExAC-tsSNP99.6-tsINDEL95/QC_plots/Bloomfield_9810_IDs_WES_4940_chr", CHR[i], ".vcf.gz-norm-ref-annotate.vcf-SNPs.vcf_recalibrated_SNPs.vcf.table"), header = TRUE, stringsAsFactors = F)
## Get PASS
Bloomfield.table.WES.FAILED<- Bloomfield.table.WES[!grepl("^PASS$", Bloomfield.table.WES$FILTER),]
Bloomfield.table.WES<- Bloomfield.table.WES[grepl("^PASS$", Bloomfield.table.WES$FILTER),]
head(Bloomfield.table.WES)
head(Bloomfield.table.WES)
  
Bloomfield.table.WES$KEY <- paste(Bloomfield.table.WES$CHROM, Bloomfield.table.WES$POS, Bloomfield.table.WES$REF, Bloomfield.table.WES$ALT, sep = ":")
colnames(Bloomfield.table.WES) <- paste0("WES_", colnames(Bloomfield.table.WES))
  
MERGED.AQ.BL <- merge(Bloomfield.table.WGS, Bloomfield.table.WES, by.x="KEY", by.y="WES_KEY", all = T)
  
## Add missingness    
header <- read.table("/40/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/01-Bloomfield-preQC/03-PLINK-QC-files/Bloomfield_9810_SNPS_INDELS_picard_biallelic.vcf-plink-jan-28-2022-hwe-geno0.05-mind0.1-missing-WXS.missing", header = TRUE, nrow = 1)
MISSING.dat <- fread("/40/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/01-Bloomfield-preQC/03-PLINK-QC-files/Bloomfield_9810_SNPS_INDELS_picard_biallelic.vcf-plink-jan-28-2022-hwe-geno0.05-mind0.1-missing-WXS.missing", header = F, sep = "\t")
MISSING.dat <- setnames(MISSING.dat, colnames(header))
head(MISSING.dat)

MISSING.dat$PASS[MISSING.dat$P<=0.05e-8] <- "FAIL"
MISSING.dat$PASS[MISSING.dat$P>0.05e-8] <- "PASS"

MERGED.AQ.BL$KEY <- as.character(MERGED.AQ.BL$KEY)

sum(MERGED.AQ.BL$KEY %in% MISSING.dat$SNP)

MERGED.AQ.BL$PASS.var <- MISSING.dat$PASS[match(MERGED.AQ.BL$KEY, MISSING.dat$SNP)]

MERGED.AQ.BL$FAIL <- ifelse(grepl("FAIL", MERGED.AQ.BL$PASS.var), "FAIL", "PASS")

for(j in 1:length(TYPE)) {

print (paste0("Doing CHR_",CHR[i], "_QC_", TYPE[j] ))
MERGED.AQ.BL <- setNames(cbind.data.frame(MERGED.AQ.BL$KEY, MERGED.AQ.BL$POS, MERGED.AQ.BL[TYPE[j]], MERGED.AQ.BL[paste0("WES_", TYPE[j])], MERGED.AQ.BL$FAIL), c("KEY", "POS", paste0("WGS_", TYPE[j]), paste0("WES_",TYPE[j]), "FAIL"))

library(reshape)
tt <- reshape::melt(data = MERGED.AQ.BL, id = c("KEY", "POS", "FAIL"))
tt <- tt[order(tt$POS),]
tt$variable <- as.character(tt$variable)

# tt <- tt[1:30000,]
tt <- tt[!is.na(tt$value),]
tt$value <- as.numeric(tt$value)


P <- ggplot(tt, aes(x=POS, y=value, color = FAIL, shape = variable)) +
  geom_point(position=position_jitter(h=0.1, w=0.1), alpha = 0.5, size = 2) +
  scale_color_manual(values = c('blue', 'red')) +
  scale_shape_manual(values= c(23, 24)) +
  # facet_grid(cols = vars(variable)) +
  theme_bw()
  
ggsave(paste0("BLOOMFIELD_AN_PLOT_", CHR[i], "_", TYPE[j], ".jpg"), plot = P, device = NULL, scale = 1, width = 16, height = 9, dpi = 300, limitsize = TRUE)
}
}
  
##############################
## Plot Missingness Aquilla ##
##############################

CHR=19:22
# TYPE=c("QD", "MQ")
TYPE="AN"
i=1
j=1

for(i in 1:length(CHR)){
  Aquilla.table.WGS <- read.table(paste0("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-a/01-Aquilla-preQC/01-VQSR-ExAC-tsSNP99.6-tsINDEL95/QC_plots/Aquilla_7983_IDs_WGS_3044_chr", CHR[i],".vcf.gz-norm-ref-annotate.vcf-SNPs.vcf_recalibrated_SNPs.vcf.table"), header = TRUE, stringsAsFactors = F)
  ## Get PASS
  Aquilla.table.WGS.FAILED<- Aquilla.table.WGS[!grepl("^PASS$", Aquilla.table.WGS$FILTER),]
  Aquilla.table.WGS<- Aquilla.table.WGS[grepl("^PASS$", Aquilla.table.WGS$FILTER),]
  head(Aquilla.table.WGS)
  
  Aquilla.table.WGS$KEY <- paste(Aquilla.table.WGS$CHROM, Aquilla.table.WGS$POS, Aquilla.table.WGS$REF, Aquilla.table.WGS$ALT, sep = ":")
  
  Aquilla.table.WES <- read.table(paste0("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-a/01-Aquilla-preQC/01-VQSR-ExAC-tsSNP99.6-tsINDEL95/QC_plots/Aquilla_7983_IDs_WES_4941_chr", CHR[i],".vcf.gz-norm-ref-annotate.vcf-SNPs.vcf_recalibrated_SNPs.vcf.table"), header = TRUE, stringsAsFactors = F)
  ## Get PASS
  Aquilla.table.WES.FAILED<- Aquilla.table.WES[!grepl("^PASS$", Aquilla.table.WES$FILTER),]
  Aquilla.table.WES<- Aquilla.table.WES[grepl("^PASS$", Aquilla.table.WES$FILTER),]
  head(Aquilla.table.WES)
  head(Aquilla.table.WES)
  
  Aquilla.table.WES$KEY <- paste(Aquilla.table.WES$CHROM, Aquilla.table.WES$POS, Aquilla.table.WES$REF, Aquilla.table.WES$ALT, sep = ":")
  colnames(Aquilla.table.WES) <- paste0("WES_", colnames(Aquilla.table.WES))
  
  MERGED.AQ.BL <- merge(Aquilla.table.WGS, Aquilla.table.WES, by.x="KEY", by.y="WES_KEY", all = T)
  
  ## Add missingness    
  header <- read.table("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-a/01-Aquilla-preQC/03-PLINK-QC-files_all_Aquilla/Aquilla_7983_WXS_SNPS-INDELS-geno0.05-hwe-mind0.1-WXSmissing.missing", header = TRUE, nrow = 1)
  MISSING.dat <- fread("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-a/01-Aquilla-preQC/03-PLINK-QC-files_all_Aquilla/Aquilla_7983_WXS_SNPS-INDELS-geno0.05-hwe-mind0.1-WXSmissing.missing", header = F, sep = "\t")
  MISSING.dat <- setnames(MISSING.dat, colnames(header))
  head(MISSING.dat)
  
  MISSING.dat$PASS[MISSING.dat$P<=0.05e-8] <- "FAIL"
  MISSING.dat$PASS[MISSING.dat$P>0.05e-8] <- "PASS"
  
  MERGED.AQ.BL$KEY <- as.character(MERGED.AQ.BL$KEY)
  
  sum(MERGED.AQ.BL$KEY %in% MISSING.dat$SNP)
  
  MERGED.AQ.BL$PASS.var <- MISSING.dat$PASS[match(MERGED.AQ.BL$KEY, MISSING.dat$SNP)]
  
  MERGED.AQ.BL$FAIL <- ifelse(grepl("FAIL", MERGED.AQ.BL$PASS.var), "FAIL", "PASS")
  
  for(j in 1:length(TYPE)) {
    
    print (paste0("Doing CHR_",CHR[i], "_QC_", TYPE[j] ))
    MERGED.AQ.BL <- setNames(cbind.data.frame(MERGED.AQ.BL$KEY, MERGED.AQ.BL$POS, MERGED.AQ.BL[TYPE[j]], MERGED.AQ.BL[paste0("WES_", TYPE[j])], MERGED.AQ.BL$FAIL), c("KEY", "POS", paste0("WGS_", TYPE[j]), paste0("WES_",TYPE[j]), "FAIL"))
    
    library(reshape)
    tt <- reshape::melt(data = MERGED.AQ.BL, id = c("KEY", "POS", "FAIL"))
    tt <- tt[order(tt$POS),]
    tt$variable <- as.character(tt$variable)
    
    # tt <- tt[1:30000,]
    tt <- tt[!is.na(tt$value),]
    tt$value <- as.numeric(tt$value)
    
    
    P <- ggplot(tt, aes(x=POS, y=value, color = FAIL, shape = variable)) +
      geom_point(position=position_jitter(h=0.1, w=0.1), alpha = 0.5, size = 2) +
      scale_color_manual(values = c('blue', 'red')) +
      scale_shape_manual(values= c(23, 24)) +
      # facet_grid(cols = vars(variable)) +
      theme_bw()
    
    ggsave(paste0("AQUILLA_AN_PLOT_", CHR[i], "_", TYPE[j], ".jpg"), plot = P, device = NULL, scale = 1, width = 16, height = 9, dpi = 300, limitsize = TRUE)
  }
}
  
  