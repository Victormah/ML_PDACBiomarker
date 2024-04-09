# ML_PDACBiomarker
[![DOI](https://zenodo.org/badge/748142894.svg)](https://zenodo.org/doi/10.5281/zenodo.10949601)

## Project Title : 
Robust and Consistent Biomarker Identification of Pancreatic Ductal Adenocarcinoma Metastasis by Machine Learning Approach

## Description
We provided the full R codes and data for biomarker discovery using RNA seq data from cancer patients, including R code, Rmarkdown codes, and Tutorials (HTML) in 5 sections. This is a part of our publication titled “Robust and Consistent Biomarker Identification of Pancreatic Ductal Adenocarcinoma Metastasis by Machine Learning Approach” in BMC Medical Informatics and Decision Making Journal. 

## Installation
Please follow the tutorials and R codes; we have a custom function that includes many R packages that need to be installed and packages in each section. 

## Usage
You should follow the R codes by following the sections,
1.TMM_norm : Data cleaning and TMM normalisation
2.Gene_filter : Gene filtering by a criteria
3.Data_integration : Data integration and batch effect correction
4.Variable selection : A combination of variable selection methods
5.Model : Random Forest model (RF) in training and validation 
All data can be accessed in the folder "data_R". Some large files are available here https://drive.google.com/drive/folders/1PapRwJpQDnngSl2IoGPZi7dOGpz_DzPL?usp=sharing , where you can use preloaded data.Rdata files : Data_before_integration.Rdata is used for section 1-3 and Data_after_integration.Rdata is used for section 4-5.
A custom function for all sections is "custom_functions.r"
## Code Explanation
You can find all R scripts in R folder. All parameters used in all R scripts can be appraised in the custom_funciton.r. The details of the parameters are explained in the scripts. Please find a detailed explanation in the “Tutorials” folder.

## Data sources and processing
The raw RNA sequencing data (all in gene levels) used in this project were downloaded and processed individually. The TCGA-PAAD and GSE79668 were already in gene expression table format ( rows = Ensembl Gene IDs and columns = Sample IDs). They could be used in data cleaning and data preprocessing. However, PACA-AU , PACA-CA, and CPTAC-PDAC datasets were not ready to use. We created the gene expression table using dcast (),a function in the reshape2 R package that transforms data from long to wide format. The function takes a data frame and a formula as input and returns a new data frame with the data in the specified format. We also identified the redundant gene IDs in datasets by checking for duplicates in IDs. Then we aggregated the data across all instances of that gene by mean of the expression values across all instances of redundant gene IDs.

## Overview of the workflow
![image](https://github.com/Victormah/ML_PDACBiomarker/assets/54091551/3618319a-a9b1-4424-bddf-c2af2367b07e)

Analysis workflow including data pre-processing, variable selection, and classification model.The model training is shown in the blue line,while model validation is shown in the pink line.

## Output
You can find the output from codes in folder “Outputs” including variable selection results and RF model performance.

## Contact 
More information and questions can be contacted at Dr Eva Caamaño Gutiérrez @EvaCaamano (caamano@liverpool.ac.uk) and Tanakamol Mahawan (m.tanakamol@gmail.com and t.mahawan@liverpoo.ac.uk)



