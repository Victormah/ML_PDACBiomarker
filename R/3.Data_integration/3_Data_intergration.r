

#title: "Tutorial3 : Data integration and batch effect correction "
#author: "Tanakamol Mahawan"
#date: "2023-12-12"

#### Introduction #####
  
#Data integration and batch effect correction are crucial stages in the processing of RNA-seq data. Data integration refers to the procedure of combining distinct datasets derived from diverse sources or experiments into an integrated representation. Batch effect correction involves eliminating non-biological variability in the data, such as technological discrepancies between samples or batches.

#We will use an ARSyNbac function of MultiBaC package to integrate and correct the batch effects of 5 datasets. The ARSyNbac function is part of the MultiBaC package, which is used for multi-omic, multi-batch correction. This function is specifically designed to remove batch effects or unwanted noise from a single omic data matrix.

#The ARSyNbac function offers 4 modes of correction, we will use mode 1 which perform only batch effect correction without noise reduction. 


##### 1.Loading required libraries #####

#Before we begin, make sure to install and load the required packages:
  
#BiocManager::install("MultiBaC")
library(MultiBaC)
library(edgeR)
library(limma)
library(MultiAssayExperiment)
library(dplyr)
library(ggplot2)
library(patchwork)
library(gridExtra)
#Load all custom functions
source("custom_functions.r")

#### 2. Load data ####
#First, we will load all 5 data and phenotype.Then, we need to check dimension if they have all same #genes. 

# Read phenotype data
pheno <- read.csv("phenotype_all_5data.csv", header = TRUE, stringsAsFactors = FALSE)
head(pheno)

# Read and process TCGA data
TCGA <- read.csv('./normLog_TCGA.csv', header = TRUE, stringsAsFactors = FALSE, row.names = 1) 
TCGA$class <- as.factor(TCGA$class)
pheno_TCGA <- pheno[pheno$SAMPLE == "A", c(1, 3, 5, 6, 7)]
head(pheno_TCGA)

# Read and process AUS data
AUS <- read.csv('./normLog_AUS.csv', header = TRUE, stringsAsFactors = FALSE, row.names = 1) 
AUS$class <- as.factor(AUS$class)
pheno_AUS <- pheno[pheno$SAMPLE == "B", c(1, 3, 5, 6, 7)]
head(pheno_AUS)

# Read and process CA data
CA <- read.csv('./normLog_CA.csv', header = TRUE, stringsAsFactors = FALSE, row.names = 1)
CA$class <- as.factor(CA$class)
pheno_CA <- pheno[pheno$SAMPLE == "C", c(1, 3, 5, 6, 7)]
head(pheno_CA)

# Transpose and convert data frames
TCGA <- as.data.frame(t(TCGA[, -1]))
AUS <- as.data.frame(t(AUS[, -1]))
CA <- as.data.frame(t(CA[, -1]))

# Read and process CPTAC data
CPTAC <- read.csv(file = "normLog_CPTAC.csv", row.names = 1, header = TRUE)
CPTAC <- as.data.frame(t(CPTAC[, -1]))

# Read and process GSE data
GSE <- read.csv(file = "normLog_GSE79668.csv", row.names = 1, header = TRUE)
GSE <- as.data.frame(t(GSE[, -1]))

# Fold change data from previous steps
FC_fil_genes <- read.csv("FC_filtered_ABC.csv", row.names=1, header=TRUE)
fil_genes <- FC_fil_genes[FC_fil_genes$abs == 3, ] #filter only same expression pattern
dim(fil_genes) # 2183 


# 3.Filter low counts with Fold Change

ArSyn can handle multiple data with limited sizes e.g. cannot use if data has more than 2,000 genes. We need to reduce amount of genes by filtering very low counts.


abc <- fil_genes[, 1:4]
abc$abs_TCGA <- abs(abc$TCGA)
abc$abs_AUS <- abs(abc$AUS)
abc$abs_CA <- abs(abc$CA)

# Convert to numeric
i <- c(2:7)
abc[, i] <- apply(abc[, i], 2, function(x) as.numeric(as.character(x)))
abc$mean <- apply(abc[, 5:7], 1, mean)

# Filter very low FC genes using 0.1
sel_abc <- abc[abc$mean >= 0.1, ]
dim(sel_abc) # 1038    8
head(sel_abc)

# Selected genes
selected <- sel_abc$gene
length(selected) # total genes = 1038 genes

#### 4.Extract only filtered genes #### 

#We will subset selected genes (n=1038) from the previous step in all data.

A.rna <- TCGA
A.rna <- A.rna[ row.names(A.rna) %in% selected, ]
dim(A.rna) # 1038  142

B.rna <- AUS
B.rna <- B.rna[ row.names(B.rna) %in% selected, ]
dim(B.rna) # 1038   54

C.rna <- CA
C.rna <- C.rna[ row.names(C.rna) %in% selected, ]
dim(C.rna) # 1038 114

D.rna <- CPTAC
D.rna <- D.rna[ row.names(D.rna) %in% selected, ]
dim(D.rna) # 1038 133

E.rna <- GSE
E.rna <- E.rna[ row.names(E.rna) %in% selected, ]
dim(E.rna) #1038   51 

A_pheno <- pheno[pheno$SAMPLE== "A",]
B_pheno <- pheno[pheno$SAMPLE== "B",]
C_pheno <- pheno[pheno$SAMPLE== "C",]
D_pheno <- pheno[pheno$SAMPLE== "D",]
E_pheno <- pheno[pheno$SAMPLE== "E",]

# All rows need to be exactly the same number 1,038

table(pheno$SAMPLE) #Check sample sizes in all data

#### 5.Check data before correction with PCA plots ####

ABC_before <- cbind(A.rna,B.rna,C.rna)
DE_before <- cbind(D.rna,E.rna)

phenoABC <- rbind(A_pheno,B_pheno,C_pheno)
phenoDE <- rbind(D_pheno,E_pheno)


# ABC data
PCA_plot(data = t(ABC_before), group = phenoABC$SAMPLE, 
         title = "ABC before correction")
PCA_plot(data = t(ABC_before), group = phenoABC$CLASS, 
         title = "ABC before correction")

# DE data
PCA_plot(data = t(DE_before), group = phenoDE$SAMPLE, 
         title = "DE before correction")
PCA_plot(data = t(DE_before), group = phenoDE$CLASS, 
         title = "DE before correction")

#### 6. Create list object for ARSyN ####

#This createMbac function creates a list object to be used by MultiBaC function 
#from a set of matrix R objects.

data_RNA<- createMbac(inputOmics = list(A.rna, B.rna, C.rna , D.rna,E.rna), 
                      batchFactor = c("A", "B", "C" , "D","E"),
                      experimentalDesign = list("A" = A_pheno$CLASS,                                                  
                                                "B" = B_pheno$CLASS,
                                                "C" = C_pheno$CLASS,
                                                "D" = D_pheno$CLASS,
                                                "E" = E_pheno$CLASS),
                      omicNames = "RNA")

#Check the data type
summary(data_RNA) # 5 data are "RNA" (all is RNA seq data)

#### 7. Batch effect correction using ARSyN mode 1 ####

# Apply Arsyn
# We will use  arsyn_1 (batch correction only) 
# Importantly, set batchEstimation = TRUE, Interaction = FALSE

arsyn_1 <- ARSyNbac(data_RNA, modelName = "RNA", Variability = 0.95, 
                    batchEstimation = TRUE, Interaction = FALSE)

# Split ABC (train data) and DE(validaion data )
arsyn1_ABC <-cbind(arsyn_1$CorrectedData$A@ExperimentList$RNA,
                   arsyn_1$CorrectedData$B@ExperimentList$RNA,
                   arsyn_1$CorrectedData$C@ExperimentList$RNA)
dim(arsyn1_ABC)
#head(arsyn1_ABC)

arsyn1_DE <-cbind(arsyn_1$CorrectedData$D@ExperimentList$RNA,
                  arsyn_1$CorrectedData$E@ExperimentList$RNA)
dim(arsyn1_DE)
#head(arsyn1_DE)

#Phenotypes
phenoABC <- rbind(A_pheno,B_pheno,C_pheno)
phenoDE <- rbind(D_pheno,E_pheno)

#### 8.Check corrected data with PCA plots ####

# ABC data after correction
PCA_plot(data = t(arsyn1_ABC), group = phenoABC$SAMPLE, 
         title = "arsyn1_ABC after correction")

PCA_plot(data = t(arsyn1_ABC), group = phenoABC$CLASS, 
         title = "arsyn1_ABC after correction")

# DE data after correction
PCA_plot(data = t(arsyn1_DE), group = phenoDE$SAMPLE, 
         title = "arsyn1_DE after correction")

PCA_plot(data = t(arsyn1_DE), group = phenoDE$CLASS, 
         title = "arsyn1_DE after correction")
      
#### 9.Outliers detection and removal ####

#We will find that after correction, that data sources e.g.,
#ABC are not separated as before correction. We have successfully integrated and 
#corrected them together. However,there are some outliers affecting downstream analysis and modelling. 
#We need to remove them. We will use bigutilsr package (robust method) , 
#detecting outliers in genetics is the criterion of being “more than 6 Median absolute deviation(MAD) away from the median.
#(see https://www.r-bloggers.com/2019/08/detecting-outlier-samples-in-pca/)


# ABC data 
result_ABC <- OLEdetect(data = t(arsyn1_ABC), group = phenoABC)
# "C.rna_control_5"  and  "C.rna_metastasis_5" were removed 
arsyn1ABC_rmOLE <- result_ABC[[2]]
pheno_Arsyn1ABC_rmOLE <- result_ABC[[3]]

PCA_plot(data= arsyn1ABC_rmOLE, group = pheno_Arsyn1ABC_rmOLE$CLASS , 
         title = "arsyn1_ABC after correction")

#DE data 
result_DE <- OLEdetect(data = t(arsyn1_DE), group = phenoDE)
# "E.rna_T_07_10_A153a","C3L.02610" ,and "C3L.00819"were removed 
arsyn1DE_rmOLE <- result_DE[[2]]

pheno_Arsyn1DE_rmOLE <- result_DE[[3]]

PCA_plot(data= arsyn1DE_rmOLE, group = pheno_Arsyn1DE_rmOLE$CLASS , 
         title = "arsyn1_DE after correction")

#### 10. Save files and sessions #####

#Save files for the next step , variable selection and modelling. 

write.csv(arsyn1ABC_rmOLE, "arsyn1ABC_rmOLE.csv")
write.csv(arsyn1DE_rmOLE, "arsyn1DE_rmOLE.csv")
write.csv(pheno_Arsyn1ABC_rmOLE,"pheno_Arsyn1ABC.csv")
write.csv(pheno_Arsyn1DE_rmOLE,"pheno_Arsyn1DE.csv")

save.image("3_Data_integration.rdata")



