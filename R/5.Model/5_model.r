
library(smotefamily)
library(MLmetrics)
library(mlr3measures)

# Loads all data using in modelling 
load("Data_after_integration.rdata")

##### Model, Train and test model by ABC data ##########
pdac_ABC <- read.csv("arsyn1ABC_rmOLE.csv",header=TRUE,stringsAsFactors = F,row.names = 1) #read csv into a dataframe
dim(pdac_ABC) # 308 1038

pheno_ABC <- read.csv("pheno_Arsyn1ABC",header = T, stringsAsFactors = T)
dim(pheno_ABC) #308   7
head(pheno_ABC)

#Overlapping genes
Overlap_genes <- read.csv("Overlap_genes.csv",header = T)
Overlap_genes <- geneset$Var1
# 15 genes 

pdac_ABC <- data.frame(pheno_ABC$CLASS,pdac_ABC )
names(pdac_ABC)[1] <- "CLASS"
pdac_ABC$CLASS <- as.character(pdac_ABC$CLASS)
pdac_ABC$CLASS <- as.factor(pdac_ABC$CLASS)

res_ABC  <- Model_RF(pdac_ABC = pdac_ABC , geneset = Overlap_genes , nround = 100)
View(res_ABC[,101:103])
write.csv(res_ABC,"res_ABC.csv" )



#######  Model, Train and test model by DE data ##########
pdac_DE <- read.csv("arsyn1DE_rmOLE.csv",header=TRUE,stringsAsFactors = F,row.names = 1) #read csv into a dataframe
dim(pdac_DE) 

pheno_DE <- read.csv("pheno_Arsyn1DE",header = T, stringsAsFactors = T)
dim(pheno_DE) 
head(pheno_DE)

pdac_DE <- data.frame(pheno_DE$CLASS,pdac_DE )
names(pdac_DE)[1] <- "CLASS"
pdac_DE$CLASS <- as.character(pdac_DE$CLASS)
pdac_DE$CLASS <- as.factor(pdac_DE$CLASS)

res_DE  <- Model_RF(pdac_DE = pdac_DE , geneset = Overlap_genes , nround = 100)
View(res_DE[,101:103])
write.csv(res_DE,"res_DE.csv" )

########## Train by ABC data and test model by DE data ############

train_df <- pdac_ABC
test_df <-  pdac_DE
geneset <- Overlap_genes

#convert CLASS to factor 
train_df$CLASS <- as.factor(train_df$CLASS )
test_df$CLASS <- as.factor(test_df$CLASS )

# Build RF model 
train.data  <- train_df[ ,which( colnames(train_df) %in% geneset ) ] 
dim(train.data)

test.data  <- test_df[ ,which( colnames(test_df) %in% geneset) ] 
dim(test.data)

train.data  <- cbind(train_df$CLASS , train.data)
names(train.data)[1] <- "CLASS"

test.data  <- cbind(test_df$CLASS , test.data)
names(test.data)[1] <- "CLASS"

# 5 fold validation 
ctrlspecs <- trainControl(method = "cv" , number = 5,
                          savePredictions = "all",
                          classProbs = TRUE) 

#We use ADASYN to balance minority class (non-metastasis)

genData_ADAS = ADAS(X= train.data[,-1], target = train.data$CLASS ,K=5  )
train_df_adas <- genData_ADAS[["data"]]  
# keep CLASS as a df , last column
CLASS_df <- train_df_adas[ncol(train_df_adas)]
# move CLASS to first column 
train_df_adas <- cbind(train_df_adas[ncol(train_df_adas)], train_df_adas[,-ncol(train_df_adas)] )
names(train_df_adas)[1] <- "CLASS"
train.data <- train_df_adas
dim(train.data)
table(train.data$CLASS) #Check balance of class distribution again

# Build RF model using ranger 
# ranger and caret support parallel processing, but the level of control and 
#ease of use may vary. 
#ranger is explicitly designed for parallel computation, 
# which might make it more straightforward for users who prioritize parallelism.

orig_fit <- caret::train(CLASS ~ .,
                         data = train.data,
                         method = "ranger",
                         verbose = FALSE,
                         trControl = ctrlspecs,
                         metric = 'Accuracy',
                         tuneLength = 10)

# Prediction 
prob_RF <- predict(orig_fit, test.data, type="prob")
pred_RF <- predict(orig_fit, test.data)
res_RF <- Res_mod(prob = prob_RF , pred = pred_RF ,test_data = test.data)
res_validate <- data.frame(rownames(resDE),res_RF)
print(res_validate)

write.csv(res_validate,"res_validate.csv" )

save.image("5_model.rdata")
