

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

library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(lmerTest)
metadata <- read.csv("../data/Jena_mouse_clean_RNA.csv")
count_data <- read.csv("../data/otu_count_clean.csv", sep = ",", row.names = 1, check.names = F)

# Returns normalized count data for a specific organ
filter_per_organ = function (metadata, countdata, organ) {
  metadata_organ = metadata %>% filter(TissueID == organ)
  filtered_count_data = countdata[, metadata_organ$PatID]
  filtered_count_data <- apply(filtered_count_data, 2, function(col) {col/sum(col)})
  colnames(filtered_count_data) = gsub("_.", "", colnames(filtered_count_data))
  return(filtered_count_data)
}



cecum_counts = filter_per_organ(metadata, count_data, "Cecum")
cecum_counts = cecum_counts[-nearZeroVar(t(cecum_counts)), ]
colon_counts = filter_per_organ(metadata, count_data, "Colon")
colon_counts = colon_counts[-nearZeroVar(t(colon_counts)), ]
stool_counts = filter_per_organ(metadata, count_data, "stool")
stool_counts = stool_counts[-nearZeroVar(t(stool_counts)), ]

# Load 16S to MAG mapping
OTUtoMAG <- read.csv("../data/VsearchMap6Out99.tsv", sep = "\t", header = F, check.names = F)
OTUtoMAG = OTUtoMAG %>% select(V1, V2)
OTUtoMAG$V2 = gsub("::.*", "", OTUtoMAG$V2)

cecumMAGCounts <- cecum_counts[rownames(cecum_counts) %in% OTUtoMAG$V1, ]
colonMAGCounts <- colon_counts[rownames(colon_counts) %in% OTUtoMAG$V1, ]
stoolMAGCounts <- stool_counts[rownames(stool_counts) %in% OTUtoMAG$V1, ]
 

calculateRXNAbundances <- function (organMAGCounts)
{
  organOTUtoMAG <- merge(organMAGCounts, OTUtoMAG, by.x = 0 , by.y = "V1") %>% select("Row.names", "V2")
  length(unique(organOTUtoMAG$V2)) == length(organOTUtoMAG$V2)
  organMAGCounts <- merge(organMAGCounts, OTUtoMAG, by.x = 0 , by.y = "V1") 
  organMAGCounts$Row.names = organMAGCounts$V2
  organMAGCounts = organMAGCounts %>% select(-V2)
  result <- as.data.frame(organMAGCounts%>%
    group_by(Row.names) %>%
    summarize(across(matches(".*/"), sum, .names = "Sum_{.col}")))
  
  # Subset incidence matrix
  mtx_rxnInModels <- as.matrix(as.data.frame(mtx_rxnInModels)[, colnames(mtx_rxnInModels) %in% result$Row.names] %>% mutate_all(as.numeric))
  organFinalOTUCount = result[result$Row.names %in% colnames(mtx_rxnInModels), ]
  rownames(organFinalOTUCount) = organFinalOTUCount$Row.names
  organFinalOTUCount = as.matrix(organFinalOTUCount %>% select(-Row.names))
  organFinalActivateOTUCount = mtx_rxnInModels %*% organFinalOTUCount 
  colnames(organFinalActivateOTUCount) = gsub("Sum_", "", colnames(organFinalActivateOTUCount))
  return(organFinalActivateOTUCount)
}

cecum = calculateRXNAbundances(cecum_counts)
stool = calculateRXNAbundances(stool_counts)
colon = calculateRXNAbundances(colon_counts)

write.csv(cecum, "../data/cecum_rxn_abundance.csv")
write.csv(colon, "../data/colon_rxn_abundance.csv")
write.csv(stool, "../data/stool_rxn_abundance.csv")