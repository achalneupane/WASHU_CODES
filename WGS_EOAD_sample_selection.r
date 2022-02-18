#########################
## Read phenotype file ##
#########################
setwd("/100/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/01-Bloomfield-preQC/03-PLINK-QC-files")
PHENO <- read.delim("/100/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/01-Bloomfield-preQC/03-PLINK-QC-files/Bloomfield_8751_metadata_20211203.csv", header = T, sep = ",", stringsAsFactors = T)
dim(PHENO)
# # 8751
# head(PHENO)

FAM <- read.table("Bloomfield_9810-hwe-geno0.05-mind0.1-WXSm.fam", header = F)
FAM$V2 <- as.character(FAM$V2)
dim(PHENO)
# [1] 3931   26
sum(PHENO$Bloomfield.gvcf.id..SM...9810. %in% FAM$V2)
# 8750 


## FIX # N/A characters
PHENO[PHENO=="#N/A"] <- NA

# View(PHENO)

##########################################
## Function to create demographic table ##
##########################################
covars <- PHENO

## Recode APOE 
covars$APOE4ANY <- as.character(covars$APOE)
covars$STATUS <- as.character(covars$STATUS..CC.)
covars$SEX <- as.character(covars$Pheno_SEX)
covars$AGE_LAST_VISIT <- as.numeric(as.character(covars$ALA))
covars$AGE_AT_ONSET <- as.numeric(as.character(covars$AAO))

covars$Bloomfield.gvcf.id..SM...9810. <- as.character(covars$Bloomfield.gvcf.id..SM...9810.)
sum(grepl ("MAP_|^27_", covars$Bloomfield.gvcf.id..SM...9810.))
# 1671

## Keep MAP only
covars <- covars[grepl ("MAP_|^27_", covars$Bloomfield.gvcf.id..SM...9810.),]

## Keep WGS only
covars <- covars[grepl("WGS", covars$WXS),]
dim(covars)
# 1039

## Let's recode APOE4
sum(grepl("22|23|33", covars$APOE4ANY))
# 253
covars$APOE4ANY[grepl("22|23|33|32", covars$APOE4ANY)] <- 0
covars$APOE4ANY[grepl("24|34|44|42|43", covars$APOE4ANY)] <- 1


Get_STATs <- function(covars){
  # covars$ETHNICITY <- as.character(covars$ETHNICITY)
  TOTAL = nrow(covars)
  # CA,CO
  N.controls <- sum(covars$STATUS==1, na.rm = T)
  N.cases <- sum(covars$STATUS==2, na.rm = T)
  
  # Number of CA <= 65 and 70 (Cases==2)
  CA.65 <- sum(as.vector(na.omit(as.numeric(as.character(covars[covars$STATUS == 2,"AGE_AT_ONSET"])))) <= 65 & 
                 as.vector(na.omit(as.numeric(as.character(covars[covars$STATUS == 2,"AGE_AT_ONSET"])))) > 0)
  
  CA.70 <- sum(as.vector(na.omit(as.numeric(as.character(covars[covars$STATUS == 2,"AGE_AT_ONSET"])))) <= 70 & 
                 as.vector(na.omit(as.numeric(as.character(covars[covars$STATUS == 2,"AGE_AT_ONSET"])))) > 0)
  
  # Number of CO > 70 (Cases==2)
  CO.70 <- sum(as.vector(na.omit(as.numeric(as.character(covars[covars$STATUS == 1,"AGE_LAST_VISIT"])))) > 70)
  
  CO.80 <- sum(as.vector(na.omit(as.numeric(as.character(covars[covars$STATUS == 1,"AGE_LAST_VISIT"])))) > 80)
  
  # Number of Others
  N.OTHERS <- sum(covars$STATUS == -9| is.na(covars$STATUS) )
  
  # Percent Female
  PERC.FEMALE <- (sum(covars$SEX == 2, na.rm = T)/(sum(covars$SEX == 1, na.rm = T)+ sum(covars$SEX == 2, na.rm = T))) *100
  
  # MISSING AGES CA and CO
  N.CA.missing.age <- sum(is.na(covars [covars$STATUS == 2, "AGE_AT_ONSET"]))
  N.CO.missing.age <- sum(is.na(covars [covars$STATUS == 1, "AGE_LAST_VISIT"]))
  
  
  # Percent APOE4
  POS <- sum(covars[, c("APOE4ANY")] == 1, na.rm = T)
  NEG <- sum(covars[, c("APOE4ANY")] == 0, na.rm = T)
  UNK <- sum(covars[, c("APOE4ANY")] == -9| is.na(covars$APOE4ANY), na.rm = T)
  PERC.APOE <- (POS/ (POS + NEG)) *100
  
  CA.APOE.PERC <- sum(covars[ covars$STATUS == 2 , c("APOE4ANY")] == 1, na.rm = T)/ (N.cases) *100
  
  CO.APOE.PERC <- sum(covars[ covars$STATUS == 1 , c("APOE4ANY")] == 1, na.rm = T)/ (N.controls) *100
  
  # For Ethnicity
  STATS <-
    cbind(
      ETHNICITY = unique(covars$ETHNICITY),
      TOTAL = TOTAL,
      '% FEMALE' = round(PERC.FEMALE, 2),
      '% APOE' = round(PERC.APOE, 2),
      'N CONTROLS (1)' = N.controls,
      'N CASES (2)' = N.cases,
      'N CONTROLS > 70 yo' = CO.70,
      'N CONTROLS > 80 yo' = CO.80,
      'N CONTROLS missing age' = N.CO.missing.age,
      'N CASES ≤ 65 yo' = CA.65,
      'N CASES ≤ 70 yo' = CA.70,
      'N CASES missing age' = N.CA.missing.age,
      'N OTHERS (-9)' = N.OTHERS,
      '% CONTROLS (1) APOE4+' = round(CO.APOE.PERC, 2),
      '% CASES (2) APOE4+' = round(CA.APOE.PERC, 2)
    )
  
  return(STATS)
}


# Vicky: "/40/PhenotypeData/WU/MAP_Base_phenotype_file/Current_MAP_Base_phenotype/MAP_unique_Pheno_base_20191118.csv this file has the most updated data-freeze for MAP that I know of"
# I will try to find AGE and Status of those samples that are missing in Bloomfield phenotype
MAP_datafreeze <- read.delim("/40/PhenotypeData/WU/MAP_Base_phenotype_file/Current_MAP_Base_phenotype/MAP_unique_Pheno_base_20191118.csv", header = T, sep = ",")
head(MAP_datafreeze)

covars$CleanedID <- gsub("\\..*","",covars$Bloomfield.gvcf.id..SM...9810. )
sum(covars$CleanedID %in% MAP_datafreeze$Primary_ID)
# 1032 of 1039 WGS samples are in this MAP datafreeze

covars$ORIGINAL_STATUS <- covars$STATUS
covars$STATUS[is.na(covars$STATUS)| covars$STATUS == -9]  <- as.character(MAP_datafreeze$Final_CC_status[match(covars$CleanedID[is.na(covars$STATUS)| covars$STATUS == -9],  MAP_datafreeze$Primary_ID)])
# Samples with missing status in Bloomfield which can be filled from MAP datafreeze
# table(as.character(MAP_datafreeze$Final_CC_status[match(covars$CleanedID[is.na(covars$ORIGINAL_STATUS)| covars$ORIGINAL_STATUS == -9],  MAP_datafreeze$Primary_ID)]))
# Var1 Freq
# 1                    29
# 2     ADAD_carrier    2
# 3      ADAD_family    9
# 4         C9ORF72+    1
# 5               CA  103
# 6         Clin_DLB    1
# 7               CO  176
# 8              DLB   10
# 9              FTD   20
# 10     MAPT_family    1
# 11        Neuro_AD   38
# 12    Neuro_AD_DLB    2
# 13    Neuro_AD_FTD    1
# 14     Neuro_AD_PD    1
# 15        Neuro_CO    7
# 16       Neuro_FTD    4
# 17        Neuro_OT    3
# 18 Neuro_Presyntom    2
# 19              OT   65
# 20          OT(CO)   72
# 21              PD    3


################################
## fill in the missing Status ##
################################
covars$STATUS[is.na(covars$STATUS)| covars$STATUS == -9]  <- as.character(MAP_datafreeze$Final_CC_status[match(covars$CleanedID[is.na(covars$STATUS)| covars$STATUS == -9],  MAP_datafreeze$Primary_ID)])

###########################
## Finally recode STATUS ##
###########################
# Recoding CA/CO status
covars$STATUS[grepl("^CO$|^Neuro_CO$|OT\\(CO\\)", covars$STATUS)] <- 1
covars$STATUS[grepl("Clin_AD|^CA$|^Neuro_AD$|^Neuro_AD_DLB$|^Neuro_PreSymptomatic_AD$", covars$STATUS)] <- 2
covars$STATUS[!grepl("^1$|^2$", covars$STATUS)] <- -9

#############################
## fill in the missing age ##
#############################
covars$AAO <- as.numeric(as.character(covars$AAO))
covars$ALA <- as.numeric(as.character(covars$ALA))
covars$AAO[is.na(covars$AAO)]  <- as.numeric(as.character(MAP_datafreeze$age_onset[match(covars$CleanedID[is.na(covars$AAO)],  MAP_datafreeze$Primary_ID)]))
covars$ALA[is.na(covars$ALA)]  <- as.numeric(as.character(MAP_datafreeze$Age_at_last[match(covars$CleanedID[is.na(covars$ALA)],  MAP_datafreeze$Primary_ID)]))

covars$AGE_LAST_VISIT <- as.numeric(as.character(covars$ALA))
covars$AGE_AT_ONSET <- as.numeric(as.character(covars$AAO))
##############################
## fill in the missing APOE ##
##############################
covars$APOE[is.na(covars$APOE)]  <- as.numeric(as.character(MAP_datafreeze$apoe[match(covars$CleanedID[is.na(covars$APOE)],  MAP_datafreeze$Primary_ID)]))
covars$APOE4ANY[grepl("22|23|33|32", covars$APOE)] <- 0
covars$APOE4ANY[grepl("24|34|44|42|43", covars$APOE)] <- 1
##############################
## Now get the demographics ##
##############################
demographics <- t(Get_STATs(covars))
demographics


## Find cases and controls of age groups
library(data.table)
covars_controls <- as.data.table(covars[covars$STATUS == 1,])
covars_cases <- as.data.table(covars[covars$STATUS == 2,])

ageGroup.cases <- cut(covars_cases$AGE_AT_ONSET,
                      breaks=c(0, 40, 50, 60, 65, 70, 75, 80, 85, 90, Inf),
                      include.lowest=TRUE,
                      labels=c("<=40", ">40-50", ">50-60", ">60-65", ">65-70",
                               ">70-75", ">75-80", ">80-85", ">85-90", ">90"))

table(ageGroup.cases)

ageGroup.controls <- cut(covars_controls$AGE_LAST_VISIT,
                         breaks=c(0, 40, 50, 60, 65, 70, 75, 80, 85, 90, Inf),
                         include.lowest=TRUE,
                         labels=c("<=40", ">40-50", ">50-60", ">60-65", ">65-70",
                                  ">70-75", ">75-80", ">80-85", ">85-90", ">90"))

table(ageGroup.controls)

AGE.Groups <- setNames(cbind.data.frame(table(ageGroup.cases), table(ageGroup.controls))[c(1,2,4)], c("Age", "CA", "CO"))
AGE.Groups




write.table(covars, "/100/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/01-Bloomfield-preQC/03-PLINK-QC-files/WGS_MAP_samples.csv", sep ="\t", col.names = T, quote = F, row.names = FALSE)



## CO
CO <- covars[covars$STATUS == 1,]
CO <- CO[CO$AGE_LAST_VISIT > 80,]
CO <- CO[!is.na(CO$AGE_LAST_VISIT),]
dim(CO)

## CA
CA <- covars[covars$STATUS == 2,]
CA <- CA[CA$AGE_AT_ONSET <= 65 & CA$AGE_AT_ONSET > 0 ,]
CA <- CA[!is.na(CA$AGE_AT_ONSET),]
dim(CA)

EOAD.covars <- rbind.data.frame(CO,CA)

write.table(EOAD.covars, "/100/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/01-Bloomfield-preQC/03-PLINK-QC-files/EOAD_samples_CA_65_CO_80_V2.csv", sep ="\t", col.names = T, quote = F, row.names = FALSE)

############################################################
## Demographic based on the updated phenotype from Carlos ##
############################################################
PHENO.C <- read.delim("EOAD_ADSP.csv", header = T, sep = ",", stringsAsFactors = F)
dim(PHENO.C)
PHENO.C$IID <-  paste0("MAP_", PHENO.C$MAP_ID)

sum(PHENO.C$IID %in% covars$CleanedID)
# 434
PHENO.C$STATUS.Carlos <- PHENO.C$new.list
PHENO.C$STATUS.Carlos <- ifelse(grepl("CA",PHENO.C$new.list), 2,1)

# Updating Phenotype for WGS samples that are in Carlos' file 
new.Covars <- merge(covars, PHENO.C, by.x ='CleanedID', by.y ='IID', all.x = T)
new.Covars$STATUS[!is.na(new.Covars$STATUS.Carlos)] <- new.Covars$STATUS.Carlos[!is.na(new.Covars$STATUS.Carlos)]
new.Covars$AGE_AT_ONSET[!is.na(new.Covars$STATUS.Carlos)] <- new.Covars$age_onset[!is.na(new.Covars$STATUS.Carlos)]
new.Covars$AGE_LAST_VISIT[!is.na(new.Covars$STATUS.Carlos)] <- new.Covars$Age_at_last[!is.na(new.Covars$STATUS.Carlos)]

View(t(Get_STATs(covars)))

##########################
## Demographic for GWAS ##
##########################
GWAS.MAP.4843 <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/MAP_GWAS4843_pheno_20220215.csv", sep = ",", header =T)
dim(GWAS.MAP.4843)
head(GWAS.MAP.4843)

GWAS.MAP.4843$ID <- as.character(GWAS.MAP.4843$ID)

# Check with the GWAs ID matrix
GWAS.fam <- read.delim("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/WashU_samples/ID_matrix_hg38_Nov2021.csv", sep = ",")
dim(GWAS.fam)
GWAS.MAP.fam <- GWAS.fam[!is.na(GWAS.fam$MAP),]
GWAS.MAP.fam$MAP_ID <- paste0("MAP_", GWAS.MAP.fam$MAP)

sum(GWAS.MAP.4843$ID %in% GWAS.MAP.fam$MAP_ID)
# 4843

###########################
## Recode STATUS ##
###########################
covars <- GWAS.MAP.4843
colnames(covars) <- c("ID", "MAP_ID", "AGE_LAST_VISIT", "AGE_AT_ONSET", "STATUS", "APOE", "SEX", "cdr_final")
# Recoding CA/CO status
covars$STATUS <- as.character(covars$STATUS)
covars$STATUS[grepl("^CO$|^Neuro_CO$|OT\\(CO\\)", covars$STATUS)] <- 1
#### final_CC_Status categories to include for CA: AD, Neuro_AD, Neuro_AD_DLB. Neuro_AD_FTD, Neuro_AD_PD, Neuro_PreSymptomatic_AD
covars$STATUS[grepl("^AD$|^Neuro_AD_FTD$|^Neuro_AD_PD$|Clin_AD|^CA$|^Neuro_AD$|^Neuro_AD_DLB$|^Neuro_PreSymptomatic_AD$", covars$STATUS)] <- 2
covars$STATUS[!grepl("^1$|^2$", covars$STATUS)] <- -9

##################
## Recode APOE4 ##
##################
## Let's recode APOE4
covars$APOE4ANY <- covars$APOE
sum(grepl("22|23|33", covars$APOE4ANY))
# 2576
covars$APOE4ANY[grepl("22|23|33|32", covars$APOE4ANY)] <- 0
covars$APOE4ANY[grepl("24|34|44|42|43", covars$APOE4ANY)] <- 1

################
## Recode SEX ##
################
covars$SEX <- as.character(covars$SEX)
covars$SEX [covars$SEX == "Female"] <- 2
covars$SEX [covars$SEX == "Male"] <- 1


#######################
## GWAS demographics ##
#######################
demographics <- t(Get_STATs(covars))
View(demographics)


## Find cases and controls of age groups
library(data.table)
covars_controls <- as.data.table(covars[covars$STATUS == 1,])
covars_cases <- as.data.table(covars[covars$STATUS == 2,])

ageGroup.cases <- cut(covars_cases$AGE_AT_ONSET,
                      breaks=c(0, 40, 50, 60, 65, 70, 75, 80, 85, 90, Inf),
                      include.lowest=TRUE,
                      labels=c("<=40", ">40-50", ">50-60", ">60-65", ">65-70",
                               ">70-75", ">75-80", ">80-85", ">85-90", ">90"))

table(ageGroup.cases)

ageGroup.controls <- cut(covars_controls$AGE_LAST_VISIT,
                         breaks=c(0, 40, 50, 60, 65, 70, 75, 80, 85, 90, Inf),
                         include.lowest=TRUE,
                         labels=c("<=40", ">40-50", ">50-60", ">60-65", ">65-70",
                                  ">70-75", ">75-80", ">80-85", ">85-90", ">90"))

table(ageGroup.controls)

AGE.Groups <- setNames(cbind.data.frame(table(ageGroup.cases), table(ageGroup.controls))[c(1,2,4)], c("Age", "CA", "CO"))
AGE.Groups

# write.table(covars, "/100/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/01-Bloomfield-preQC/03-PLINK-QC-files/GWAS_MAP_samples.csv", sep ="\t", col.names = T, quote = F, row.names = FALSE)

## CO
CO <- covars[covars$STATUS == 1,]
CO <- CO[CO$AGE_LAST_VISIT > 80,]
CO <- CO[!is.na(CO$AGE_LAST_VISIT),]
dim(CO)

## CA
CA <- covars[covars$STATUS == 2,]
CA <- CA[CA$AGE_AT_ONSET <= 65 & CA$AGE_AT_ONSET > 0 ,]
CA <- CA[!is.na(CA$AGE_AT_ONSET),]
dim(CA)

GWAS.EOAD.covars <- rbind.data.frame(CO,CA)

write.table(GWAS.EOAD.covars, "/100/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/01-Bloomfield-preQC/03-PLINK-QC-files/GWAS_EOAD_samples_CA_65_CO_80_V2.csv", sep ="\t", col.names = T, quote = F, row.names = FALSE)

####################################
## Extract samples from GWAS data ##
####################################
## CO
CO <- covars[covars$STATUS == 1,]
CO <- CO[CO$AGE_LAST_VISIT > 70,]
CO <- CO[!is.na(CO$AGE_LAST_VISIT),]
dim(CO)

## CA
CA <- covars[covars$STATUS == 2,]
CA <- CA[CA$AGE_AT_ONSET <= 70 & CA$AGE_AT_ONSET > 0 ,]
CA <- CA[!is.na(CA$AGE_AT_ONSET),]
dim(CA)

GWAS.EOAD.covars <- rbind.data.frame(CO,CA)
##########################################################################
## Also adding samples that are mismatch in Vicky's and Fengxian's list ##
##########################################################################
MISMATCHES <- read.delim("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/WashU_samples/EOAD_ADSP-VF=pheno-updated.csv", header = T, sep = ",", stringsAsFactors = F)
dim(MISMATCHES)
MISMATCHES <- MISMATCHES[MISMATCHES$MATCH. == "mismatch",]
sum(MISMATCHES$MAPID %in% covars$ID)
## 113 of the mismatched are in GWAS, so including these as well.

MISMATCHES <- covars[covars$ID %in% MISMATCHES$MAPID,]
GWAS.EOAD.covars <- rbind.data.frame(MISMATCHES, GWAS.EOAD.covars)
write.table(GWAS.EOAD.covars, "/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/WashU_samples/PHENOTYPE_GWAS_EOAD_samples_CA_70_CO_70_V2_plus_mismatches.csv", sep ="\t", col.names = T, quote = F, row.names = FALSE)
######################################################
## Match with the ID matrix and get the IID and FID ##
######################################################
sum(GWAS.MAP.fam$MAP_ID %in% GWAS.EOAD.covars$ID)
## 2953

GWAS.EOAD.covars <- GWAS.MAP.fam[GWAS.MAP.fam$MAP_ID %in% GWAS.EOAD.covars$ID, c(1:2)]
write.table(GWAS.EOAD.covars, "/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/WashU_samples/GWAS_EOAD_samples_CA_70_CO_70_V2_plus_mismatches.csv", sep ="\t", col.names = T, quote = F, row.names = FALSE)


