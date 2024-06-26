
#title: "Tutorial4 : Variable selection "
#author: "Tanakamol Mahawan"
#date: "2023-12-12"

# Introduction
  
#Variable selection is the process of selecting the most relevant features for a 
#machine learning model. It’s crucial because it can enhance the model’s prediction 
#performance, provide faster and more cost-effective predictions, and offer a better 
#understanding of the data. 

#We use combination of methods : LASSO selection first (100 models) , 
#then selected genes from LASSO will be further run in Boruta and VarselRF(100 models), 
#The overlapping genes of these 2 methods will be finally selected for modelling. 

#Selection criteria : Genes must be found in at least 80% and 5 splits out of 10,


####### 1.Loading Required Libraries#######

#Before we begin, make sure to install and load the required packages:

library(smotefamily)
library(MLmetrics)
library(mlr3measures)
library(clusterProfiler)
library(org.Hs.eg.db)
#Load all custom functions
source("custom_functions.r")


#### 2.Loading data and phenotype######

# Loads all data using in gene filtering
load("Data_after_integration.rdata")

abc<- ABC_data
abcMeta<-ABC_pheno

abc <- read.csv("arsyn1ABC_rmOLE.csv",row.names = 1, header = T)
abcMeta <- read.csv("pheno_Arsyn1ABC.csv",row.names = 1, header = T)

#Check sample and class distribution 
table(abcMeta$SAMPLE)
table(abcMeta$CLASS)

#prep variables for context
#find the match for the ensembl ids
tableIDs<-bitr(colnames(abc),fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = org.Hs.eg.db)
head(tableIDs)


#### 3.Check data by PCA plot ####

PCA_plot(abc,abcMeta$CLASS, title = "ABC data")


##### 4.Preparing data splits#####

Split data into 10 splits 
```{r , results='hide'}
# Split train 90% and test 10% 
abcSplits<-applyPrepData(data_train = abc,prepData = prepData,factor_var = factor(abcMeta$CLASS),split = 0.9,n = 10)


##### 5.Apply the variable selection algorithms to each of the splits######

#LASSO is used first.

abc_LassoRes<-applyMultiVarSelectionForCrossVal_withModel(inputListwSplits = abcSplits,functionSelection =TwoClassLasso ,thresMod = 80,ntree = 500 )

abc_LassoRes_concRes<-as.data.frame(table(concatenateVarSelec(abc_LassoRes)))
abc_LassoRes_concRes$GeneSymbol<-tableIDs[match(abc_LassoRes_concRes$Var1,tableIDs[,1]),2]

#View(abc_LassoRes_concRes)

write.csv(abc_LassoRes_concRes,"abcLassoResults_FrequencyTable_relax.csv",row.names = F)
head(abc_LassoRes_concRes)

##### 6.Modelling and extracting model metrics of LASSO #####


abc_LassoRes_modelMetExtract<-modelMetrExtract(abc_LassoRes)

View(abc_LassoRes_modelMetExtract)

write.csv(abc_LassoRes_modelMetExtract,"abcLassoResultsModel_10DataSplits_relax.csv")


#### 7.Plot of LASSO selected genes####

#LASSO
selectedGenes <- abc_LassoRes_concRes[abc_LassoRes_concRes$Freq>5,1]
length(selectedGenes) #38 genes

PCA_plot(abc[,colnames(abc)%in%selectedGenes],abcMeta$CLASS, title = "Lasso selected genes (n=38)")

bp_lasso <- boxplotSubsetPDACmetastasis(soma = abc,selVars = selectedGenes[1:9],metadata = abcMeta,tableIDs = tableIDs)
#boxplots showing first 9 genes 


####  Backward selection after LASSO selection  #######

#### 8.Apply the backwards selection algorithms to each of the splits ####

#We use selected genes from LASSO to further split and select by backwards selection
#called varselRF. 

abcSplits_afterLasso<-applyPrepData(data_train = abc[,colnames(abc)%in%selectedGenes],prepData = prepData,factor_var = factor(abcMeta$CLASS),split = 0.9,n = 10)

abc_BackwardsAfterLasso <- applyMultiVarSelectionForCrossVal_withModel(
  inputListwSplits = abcSplits_afterLasso,
  functionSelection = BackSelec_wTimes,
  thresMod = 80,
  ntree = 500
)

abc_Backwards_concRes<-as.data.frame(table(concatenateVarSelec(abc_BackwardsAfterLasso)))
abc_Backwards_concRes$GeneSymbol<-tableIDs[match(abc_Backwards_concRes$Var1,tableIDs[,1]),2]
View(abc_Backwards_concRes)
write.csv(abc_Backwards_concRes,"abc_Backwards_FrequencyTable.csv",row.names = F)

selectedGenes_Backwards<-abc_Backwards_concRes[abc_Backwards_concRes$Freq>=5,1]
length(selectedGenes_Backwards) #15 genes

selectedGenes_df <- data.frame(selectedGenes_Backwards)
selectedGenes_df$GeneSymbol <-tableIDs[match(selectedGenes_Backwards,tableIDs[,1]),2]
print(selectedGenes_df)


##### 9.Modelling of Backwards selected genes #####

abc_BackwardsRes_modelMetExtract<-modelMetrExtract(abc_BackwardsAfterLasso)
View(abc_BackwardsRes_modelMetExtract)
write.csv(abc_BackwardsRes_modelMetExtract,"abc_BackwardsResmodel_10DataSplits.csv")

#### 10.Plot of Backwards selected genes####


PCA_plot(abc[,colnames(abc)%in%selectedGenes_Backwards],abcMeta$CLASS,
         title = "Backwards selected genes(n=15)")

bp_Backwards <- boxplotSubsetPDACmetastasis(soma = abc,selVars = selectedGenes_Backwards[1:9],metadata = abcMeta,tableIDs = tableIDs)
#boxplots showing first 9 genes 
 

##### Boruta after LASSO selection #######

##### 11.apply the Boruta selection algorithms to each of the splits #####

```{r ,eval=FALSE}
abc_Boruta_afterLasso <- applyMultiVarSelectionForCrossVal_withModel(
  inputListwSplits = abcSplits_afterLasso,
  functionSelection = Boruta_wTimes,
  thresMod = 80,
  ntree = 500
)

abc_Boruta_concRes<-as.data.frame(table(concatenateVarSelec(abc_Boruta_afterLasso)))

abc_Boruta_concRes$GeneSymbol<-tableIDs[match(abc_Boruta_concRes$Var1,tableIDs[,1]),2]

View(abc_Boruta_concRes)

write.csv(abc_Boruta_concRes,"abc_Boruta_FrequencyTable.csv",row.names = F)

head(abc_Boruta_concRes)

####### 12.Modelling of Boruta  selected genes #######

#Boruta
abc_BorutaRes_modelMetExtract<-modelMetrExtract(abc_Boruta_afterLasso )
View(abc_BorutaRes_modelMetExtract)
write.csv(abc_BorutaRes_modelMetExtract,"abc_BorutaResmodel_10DataSplits.csv")


##### 13.Some plots of boruta selected genes#####

selectedGenes_Boruta <- abc_Boruta_concRes[abc_Boruta_concRes$Freq>5,1]
#22 genes 

PCA_plot(abc[,colnames(abc)%in%selectedGenes_Boruta],abcMeta$CLASS,title = "Boruta selected genes(n=22)")

bp_Boruta<- boxplotSubsetPDACmetastasis(soma = abc,selVars = selectedGenes_Boruta[1:9],metadata = abcMeta,tableIDs = tableIDs)
#boxplots showing first 9 genes 

##### 14.Overlap of selectedGenes_Boruta and selectedGenes_Backwards####

#The overlapping genes will be used in modelling.


Overlap_genes <- intersect(selectedGenes_Backwards,selectedGenes_Boruta )
length(Overlap_genes) # 15 genes (all backwards genes are found in boruta)
Overlap_genes <- abc_Backwards_concRes[abc_Backwards_concRes$Var1%in%Overlap_genes,]
Overlap_genes <- Overlap_genes[order(Overlap_genes$Freq, decreasing = T),]
print(Overlap_genes)

###### Save results

write.csv(Overlap_genes,"Overlap_genes.csv")
save.image("4_variable_selection.rdata")


