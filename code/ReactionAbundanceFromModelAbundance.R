library(dplyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(caret)

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

## Normalize matrix such that the sum for each species is "1" -> species with larger genomes have 
## smaller contribution of individual reactions
#colSums(mtx_rxnInModels)
#mtx_rxnInModels <- prop.table(mtx_rxnInModels,2)
#colSums(mtx_rxnInModels)


metadata <- read.csv("../data/Jena_mouse_clean_RNA.csv")
count_data <- read.csv("../data/otu_count_clean.csv", sep = ",", row.names = 1, check.names = F)

# Returns normalized count data for a specific organ
filter_per_organ2 = function (metadata, countdata, organ) {
  # Get only the organ specific data
  metadata_organ = metadata %>% filter(TissueID == organ)
  # Get only the count data for the organ
  filtered_count_data = countdata[, metadata_organ$PatID]
  # Normalize OTU abundances
  #filtered_count_data <- apply(filtered_count_data, 2, function(col) {col/sum(col)})
  colnames(filtered_count_data) = gsub("_.", "", colnames(filtered_count_data))
  return(filtered_count_data)
}



cecum_counts = filter_per_organ2(metadata, count_data, "Cecum")
nrow(cecum_counts) == 702
cecum_counts = cecum_counts[-nearZeroVar(t(cecum_counts)), ]
nrow(cecum_counts) == 527
colon_counts = filter_per_organ2(metadata, count_data, "Colon")
nrow(colon_counts) == 702
colon_counts = colon_counts[-nearZeroVar(t(colon_counts)), ]
nrow(colon_counts) == 590
stool_counts = filter_per_organ2(metadata, count_data, "stool")
nrow(stool_counts) == 702
stool_counts = stool_counts[-nearZeroVar(t(stool_counts)), ]
nrow(stool_counts) == 575

calculateRXNAbundances <- function (organMAGCounts, normalize = F)
{
  OTUtoMAG <- read.csv("../data/VsearchMap6Out99.tsv", sep = "\t", header = F, check.names = F)
  OTUtoMAG = OTUtoMAG %>% select(V1, V2)
  OTUtoMAG$V2 = gsub("::.*", "", OTUtoMAG$V2)

  # Merge otu and mag information
  organMAGCounts <- merge(organMAGCounts, OTUtoMAG, by.x = 0 , by.y = "V1")
  organMAGCounts$Row.names = organMAGCounts$V2
  organMAGCounts = organMAGCounts %>% select(-V2)
  # Sum up same MAGs (150 MAGs left)
  result <- as.data.frame(organMAGCounts%>%
    group_by(Row.names) %>%
    summarize(across(matches(".*/"), sum, .names = "Sum_{.col}")))
  
  # The rows of the counts should be the same as the columns as the  activity matrix
  int = intersect(result$Row.names, colnames(mtx_rxnInModels))
  mtx_rxnInModels <- as.matrix(as.data.frame(mtx_rxnInModels)[int] %>% mutate_all(as.numeric))
  organFinalOTUCount = result[result$Row.names %in% int, ]
  print(all.equal(organFinalOTUCount$Row.names, colnames(mtx_rxnInModels)))
  rownames(organFinalOTUCount) = organFinalOTUCount$Row.names
  organFinalOTUCount = as.matrix(organFinalOTUCount %>% select(-Row.names))
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

abs_cecum = calculateRXNAbundances(cecum_counts)
abs_stool = calculateRXNAbundances(stool_counts)
abs_colon = calculateRXNAbundances(colon_counts)
# Base on NOT normalized zero variance-filtered MAG counts, also only MAGs only reactionsa are included that have a
# MAG in our list
write.csv(abs_cecum, "../data/absolute_cecum_rxn_abundance.csv", row.names = T)
write.csv(abs_colon, "../data/absolute_colon_rxn_abundance.csv", row.names = T)
write.csv(abs_stool, "../data/absolute_stool_rxn_abundance.csv", row.names = T)

rel_cecum = calculateRXNAbundances(cecum_counts, normalize = T)
rel_stool = calculateRXNAbundances(stool_counts, normalize = T)
rel_colon = calculateRXNAbundances(colon_counts, normalize = T)
# Base on NORMALIZED zero variance-filtered OTU counts, then mapped to mags, before counts were calculated and normalized
write.csv(rel_cecum, "../data/relative_cecum_rxn_abundance.csv", row.names = T)
write.csv(rel_colon, "../data/relative_colon_rxn_abundance.csv", row.names = T)
write.csv(rel_stool, "../data/relative_stool_rxn_abundance.csv", row.names = T)
