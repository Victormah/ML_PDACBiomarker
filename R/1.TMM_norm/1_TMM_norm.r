
#title: "Tutorial1 :RNA seq data normalisation by TMM normalisation"
#author: "Tanakamol Mahawan"
#date: "2023-12-09"
#output: html_document

  
# Introduction 
  
#RNA-seq normalisation is crucial for adjusting raw data to account for factors 
#that can affect the accuracy of analysis. These factors include sequencing depth, 
#transcript length, and variability between samples or batches. Normalisation methods help to ensure that gene expression values are comparable across different samples or conditions, improving the reliability of downstream analyses. However, the best normalisation method can vary depending on the data and research question. It’s a vital step in RNA-seq data analysis. 

###### 1.Load Packages#########

#loads necessary R packages for the tutorial requires functions from these packages for data manipulation, statistical analysis, and visualization. You may need to install some of them if they are not installed.


library(edgeR)
library(tmaptools)
library(tidyverse)
library(reshape2)

source("custom_functions.r") # source the built custom functions


###### 2. Input the Data#########

#Reads raw count data and phenotype information from CSV files.
#Make sure that the input data is raw count as EdgeR works on raw data.

#This tutorial will show you only GSE79668 data as an example, you can follow the same steps to run with other data : TCGA-PAAD, PACA-AU and PACA-CA,and CPTAC-PDAC.

#Please ensure that you use the raw RNA seq counts (not processed data e.g., RSEM , log2 normalised counts). 

# Read raw data and phenotype information
data_raw <- read.csv("GSE79668_count.csv", stringsAsFactors = FALSE)
data_raw <- data_raw[,-2]
pheno <- read.csv("pheno_GSE79668.csv", stringsAsFactors = FALSE)
head(pheno)


##### 3.Quality control#######

#We will Conduct basic quality control by inspecting dimensions and structure of the 
#raw data by subseting data to include only selected patients and removes any missing values.
#Quality control ensures that the initial dataset is suitable for analysis and removes any problematic entries.


# Display dimensions and structure of raw data
dim(data_raw)
head(data_raw)

# Subset data for selected patients and remove NAs
data_clean <- data_raw[, colnames(data_raw) %in% pheno$Sample_ID, ]
data_clean <- na.omit(data_clean)
dim(data_clean)


######## 4.Create DGEList Object########

#Next, we create a DGEList object, a specific object used in edgeR for differential gene expression analysis. Then,display the first few counts. 
#We add library sizes that can be used for Normalisation and to account for variability in library sizes.
#The DGEList object is a crucial step for downstream analysis, providing a framework for further statistical testing.


# Create DGEList object
dgList <- DGEList(data_clean)

# View the first few counts
head(dgList[["counts"]])

# Assign sample groups
samp_groups <- as.factor(pheno$CLASS)
dgList[["samples"]]$group <- samp_groups

# Add additional information to samples table
dgList[["samples"]]$lib.size



###### 5.Filtering lowly expressed genes######

#Next,we filter genes based on expression levels and log-transforms counts 
#per million using function filterByExpr. Filtering reduces noise and focuses the analysis on genes with sufficient expression levels.


# Filter genes
keep <- filterByExpr(y = dgList, design = samp_groups)
dgList <- dgList[keep, ]

# Log-transform counts per million (CPM)
log_cpm <- cpm(dgList, log = TRUE)

##### 6. Exploratory Data Analysis####

#Visualisation helps assess the datas distribution and identify potential outliers 
#or patterns.We generate a box plot and PCA plot for visualizing the data before normalisation.


# Box plot of unnormalized data
boxplot(log_cpm, main = "Unnormalized_GSE79668")

# PCA plot
data01 <- t(as.matrix(dgList)) + 1
data01 <- log2(data01)
PCA_plot(data = data01, group = pheno$CLASS, title = "GSE79668 data before normalisation")


#Normalisation using TMM

#Normalisation using the TMM (trimmed mean of M-values) method is done using the calcNormFactors function from the edgeR package. MM Normalisation involves finding a scaling factor based on the trimmed mean of M-values, which helps account for differences in library size between samples.
#The Normalisation step using the TMM method is crucial in preparing the data for downstream analyses, such as differential gene expression analysis. It improves the reliability of results and facilitates meaningful comparisons between different samples.

#More information : Robinson, M.D., Oshlack, A. A scaling normalization method for differential expression analysis of RNA-seq data. Genome Biol 11, R25 (2010). https://doi.org/10.1186/gb-2010-11-3-r25


# Normalize data using TMM
dgList <- calcNormFactors(dgList, method = "TMM")

# Display normalization factors
dgList$samples$norm.factors

# Log-transform normalized counts
normalised_counts_log <- cpm(dgList, log = TRUE, normalised.lib.sizes = TRUE)



########8. Data visualisation ##########

# Box plot of normalized data
boxplot(normalised_counts_log, main = "TMM Normalized_GSE79668")

# PCA plot of normalized data
data02 <- as.matrix(t(normalised_counts_log))
PCA_plot(data = data02, group = pheno$CLASS, title = "GSE79668 data normalized by TMM")


########9.Save Results#########


# Save normalized data
write.csv(normalised_counts_log, "normalised_log_GSE79668.csv")

# Save workspace
save.image("1_TMM_norm.rdata")




