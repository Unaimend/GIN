library(dplyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(caret)
source("utils.R")
####################
## Alternatively use FVA reactions
obj_metabolicModelsMouse <- readRDS("../data/fva_99_reactions_thr06_MetamouseHQBins20211001.RDS")

##Get list of all reactions contained in the models
lst_allRxns <- c()
for(obj_model in obj_metabolicModelsMouse){
  lst_allRxns <- union(lst_allRxns, obj_model$active)

}

##Build incidence matrix of reactions in each model
mtx_rxnInModels <- matrix(0,length(lst_allRxns),length(obj_metabolicModelsMouse))
rownames(mtx_rxnInModels) <- lst_allRxns
colnames(mtx_rxnInModels) <- names(obj_metabolicModelsMouse)
for(str_modelName in names(obj_metabolicModelsMouse)) {
  obj_model <- obj_metabolicModelsMouse[[str_modelName]]
  mtx_rxnInModels[obj_model$active,str_modelName] <- 1
}

absolute_cecum_counts = read.csv(file = "../data/cecum_mag_counts.csv")
nrow(absolute_cecum_counts) == 161
absolute_colon_counts = read.csv(file = "../data/colon_mag_counts.csv")
nrow(absolute_colon_counts) == 161
absolute_stool_counts = read.csv(file = "../data/stool_mag_counts.csv")
nrow(absolute_stool_counts) == 161
calculateRXNAbundances <- function (organMAGCounts, normalize = F)
{
  # The rows of the counts should be the same as the columns as the  activity matrix
  int = intersect(organMAGCounts$V2, colnames(mtx_rxnInModels))
  length(int) == 124
  mtx_rxnInModels <- as.matrix(as.data.frame(mtx_rxnInModels)[int] %>% mutate_all(as.numeric))
  organFinalOTUCount = organMAGCounts[organMAGCounts$V2 %in% int, ]
  nrow(organFinalOTUCount) == 124
  # 124 Models of the 161 models have reaction information
  print(all.equal(organFinalOTUCount$V2, colnames(mtx_rxnInModels)))
  rownames(organFinalOTUCount) = organFinalOTUCount$V2
  organFinalOTUCount = as.matrix(organFinalOTUCount %>% select(-V2))
  # Small example to check if the math works out
  #mtx_test = mtx_rxnInModels[1:5, 1:5]
  #organCount = organFinalOTUCount[1:5, 1:5]
  #res_test = mtx_test %*% organCount
  organFinalActivateRxNCount = mtx_rxnInModels %*% organFinalOTUCount
  colnames(organFinalActivateRxNCount) = gsub("Sum_", "", colnames(organFinalActivateRxNCount))
  if(normalize) {
    organFinalActivateRxNCount <- apply(organFinalActivateRxNCount, 2, function(col) {col/sum(col)})
  }
  return(organFinalActivateRxNCount)
}

abs_cecum = calculateRXNAbundances(absolute_cecum_counts)
abs_stool = calculateRXNAbundances(absolute_stool_counts)
abs_colon = calculateRXNAbundances(absolute_colon_counts)
# Base on NOT normalized zero variance-filtered MAG counts, also only MAGs only reactionsa are included that have a
# MAG in our list, again cecum only includes 47 Samples. For the reaction calculation only 181 are currently included.
# NOT ZERO VARIANCE FILTERING WAS APPLIED TO MAG COUNTS
write.csv(abs_cecum, "../data/absolute_cecum_rxn_abundance.csv", row.names = T)
write.csv(abs_colon, "../data/absolute_colon_rxn_abundance.csv", row.names = T)
write.csv(abs_stool, "../data/absolute_stool_rxn_abundance.csv", row.names = T)

rel_cecum = calculateRXNAbundances(absolute_cecum_counts, normalize = T)
rel_stool = calculateRXNAbundances(absolute_stool_counts, normalize = T)
rel_colon = calculateRXNAbundances(absolute_colon_counts, normalize = T)
# Base on NORMALIZED zero variance-filtered OTU counts, then mapped to mags, before counts were calculated and normalized
# NOT ZERO VARIANCE FILTERING WAS APPLIED TO MAG COUNTS
write.csv(rel_cecum, "../data/relative_cecum_rxn_abundance.csv", row.names = T)
write.csv(rel_colon, "../data/relative_colon_rxn_abundance.csv", row.names = T)
write.csv(rel_stool, "../data/relative_stool_rxn_abundance.csv", row.names = T)
