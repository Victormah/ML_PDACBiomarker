
#title: "Tutorial2 : Gene filtering and analysing gene expression data "
#author: "Tanakamol Mahawan"
#date: "2023-12-12"


  
# Introduction
  
#In this tutorial, we will guide you through the process of filtering and analysing 
#gene expression data from multiple datasets (TCGA, AUS, and CA). 
#The objective is to identify genes with consistent fold change direction and eliminate noise.

#### 1.Loading Required Libraries####

library(stats)
library(dplyr)
library(plyr)
library(Biobase)
#Load all custom functions
source("custom_functions.r")

#### 2.Load data #####

# Read  data and phenotype information
pheno<- read.csv("phenotype_all_5data.csv",header = T, stringsAsFactors = F)
head(pheno)

TCGA <- read.csv('./normLog_TCGA.csv',header=TRUE,stringsAsFactors = F,row.names = 1) 
TCGA$class <- as.factor(TCGA$class)
pheno_TCGA <- pheno[pheno$SAMPLE=="A",c(1,3,5,6,7)]
head(pheno_TCGA)

AUS <- read.csv('./normLog_AUS.csv',header=TRUE,stringsAsFactors = F,row.names = 1) 
AUS$class <- as.factor(AUS$class)
pheno_AUS <- pheno[pheno$SAMPLE=="B",c(1,3,5,6,7)]
head(pheno_AUS)

CA <- read.csv('./normLog_CA.csv',header=TRUE,stringsAsFactors = F,row.names = 1) #read
CA$class <- as.factor(CA$class)
pheno_CA <- pheno[pheno$SAMPLE=="C",c(1,3,5,6,7)]
head(pheno_CA)

# Transpose and convert data frames
TCGA <- as.data.frame(t(TCGA[, -1]))
AUS <- as.data.frame(t(AUS[, -1]))
CA <- as.data.frame(t(CA[, -1]))

#### 3.Find Overlapping Gene IDs ####

#Identifying the common set of genes across datasets is crucial for comparative 
#analysis and downstream processing.

# Find the overlap of gene IDs across the three datasets
ABC_names <- Reduce(intersect, list(row.names(TCGA),
                                    row.names(AUS),
                                    row.names(CA)))
paste0("#overlapping genes :",length(ABC_names)) #How many overlapping genes?

# Filter Data Based on Overlapping Gene IDs

TCGA <- TCGA[ row.names(TCGA) %in% ABC_names  , ]
dim(TCGA)

AUS <- AUS[ row.names(AUS) %in% ABC_names  , ]
dim(AUS)

CA <-  CA[ row.names(CA) %in% ABC_names  , ]
dim(CA)

#### 4. Calculate Mean Expression ####

# Calculate the mean expression for each gene across different classes in each dataset
# TCGA
TCGA_mean <- aggregate(t(TCGA),
                       by = list(pheno_TCGA$CLASS),
                       FUN = mean)

TCGA_mean <-  setNames(data.frame(t(TCGA_mean[,-1])), TCGA_mean[,1])
TCGA_mean$FC <- TCGA_mean$metastasis- TCGA_mean$non_metastasis

#AUS
AUS_mean <- aggregate(t(AUS),
                      by = list(pheno_AUS$CLASS),
                      FUN = mean)
AUS_mean <-  setNames(data.frame(t(AUS_mean[,-1])), AUS_mean[,1])
AUS_mean$FC <- AUS_mean$metastasis- AUS_mean$non_metastasis

#CA
CA_mean <- aggregate(t(CA),
                     by = list(pheno_CA$CLASS),
                     FUN = mean)
CA_mean <-  setNames(data.frame(t(CA_mean[,-1])), CA_mean[,1])
CA_mean$FC <- CA_mean$metastasis- CA_mean$non_metastasis

# Combine fold changes from different datasets into a single dataframe

FC_all <- as.data.frame(cbind(rownames(TCGA_mean),TCGA_mean$FC,AUS_mean$FC,CA_mean$FC ))
names(FC_all)[] <- c("gene","TCGA","AUS","CA")
head(FC_all)

##### 5. Direction of Fold Change #####

#Considering the direction of fold change helps filter genes that consistently change 
#in the same direction across datasets.

# Convert relevant columns to numeric
i <- c(2:4)
FC_all[ , i] <- apply(FC_all[ , i], 2, function(x) as.numeric(as.character(x)))
sapply(FC_all, class)  

# Calculate the direction of fold change for each dataset
FC_all$TCGA2 <- FC_all$TCGA / abs(FC_all$TCGA)
FC_all$AUS2 <- FC_all$AUS / abs(FC_all$AUS)
FC_all$CA2 <- FC_all$CA / abs(FC_all$CA)
FC_all$sum <- FC_all$TCGA2+FC_all$AUS2+FC_all$CA2
FC_all$abs <- abs(FC_all$sum)
head(FC_all)

#### 6. Filter Only Genes with the Same Direction ####

# Filter genes that have the same direction of fold change in all three datasets
fil_genes <- FC_all[FC_all$abs == 3, ]
head(fil_genes)

#### 7. Save the results #####

#fil_genes is gene lists that have save expression patterns ,  
#we will use this in next step which is data integration and variable selection. 

# Save the final filtered fold change data
write.csv(fil_genes, "FC_filtered_ABC.csv")
# Save the R environment
save.image("2_Filter_genes.rdata")


