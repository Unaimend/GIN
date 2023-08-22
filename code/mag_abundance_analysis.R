library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(caret)
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
colon_counts = filter_per_organ(metadata, count_data, "Colon")
stool_counts = filter_per_organ(metadata, count_data, "stool")

calculate_p_values <- function(counts1,counts2, counts3,age = 2)
{
  commons = intersect(rownames(counts1), rownames(counts2))
  commons = intersect(commons, rownames(counts3))
  models =  list()
  ps =  list()
  data = data.frame(OTU = character(), coefficient = numeric(),  p.val = numeric()) 
  for(otu in unique((commons)))
  {
    current_otu = as.data.frame(counts1[otu, ])
    current_otu2 = merge(current_otu,  as.data.frame(counts2[otu, ]) , by.x = 0, by.y = 0 )
    current_otu3 = merge(current_otu2,  as.data.frame(counts3[otu, ]) , by.x = "Row.names", by.y = 0 )
    colnames(current_otu3) <- c("patID", "cecum_count", "colon_count", "stool_count")
    rownames(current_otu3) <- gsub(".*/", "", current_otu3$patID)
    temp_meta = metadata %>% select(c(Age, ID))
    temp_meta$ID =  gsub("_.", "", temp_meta$ID)
    temp_meta = unique(temp_meta)
    current_otu4 = merge(current_otu3, temp_meta,  by.x = 0, by.y = "ID")
    current_otu4 = current_otu4[current_otu4$Age == age, ]
    rownames(current_otu4) = current_otu4$Row.names
    current_otu4 = current_otu4[-c(1)]
    colnames(current_otu4) <- c("ID", "1_","2_","3_","Age")
    long_data <- pivot_longer(
      data = current_otu4,
      cols = ends_with("_"),
      names_to = "Organ",
      values_to = "Count"
    )
    long_data$Organ = gsub("_", "", long_data$Organ)
    long_data$Organ = as.integer(long_data$Organ)
    long_data$ID = as.factor(long_data$ID)
    forbidden_fruits = c("Zotu217")
    if(otu %in% forbidden_fruits) {
      next
    }
    model = lmer(Count ~ Organ +(1|ID) , data = long_data)
    model_sm = summary(model)
    models[[otu]] = model_sm
    ps[[otu]] = model_sm$coefficients[, c("Estimate", "Pr(>|t|)")][2, ]
    data[nrow(data)+1, ] <- c((otu),  ps[[otu]][[1]], ps[[otu]][[2]] )
  }
  adjustNonNaNOTUS = p.adjust(data$p.val, method = "BH")
  data$p.adj = adjustNonNaNOTUS
  return(data)
}

#Two_month <- calculate_p_values(cecum_counts, colon_counts, stool_counts)

# Load 16S to MAG mapping
OTUtoMAG <- read.csv("../data/VsearchMap6Out99.tsv", sep = "\t", header = F, check.names = F)
cecumMAGCounts <- cecum_counts[rownames(cecum_counts) %in% OTUtoMAG$V1, ]
cecumMAGCounts <- apply(cecumMAGCounts, 2, function(col) {col/sum(col)})

colonMAGCounts <- colon_counts[rownames(colon_counts) %in% OTUtoMAG$V1, ]
colonMAGCounts <- apply(colonMAGCounts, 2, function(col) {col/sum(col)})

stoolMAGCounts <- stool_counts[rownames(stool_counts) %in% OTUtoMAG$V1, ]
stoolMAGCounts <- apply(stoolMAGCounts, 2, function(col) {col/sum(col)})

cecumMAGCounts = cecumMAGCounts[-nearZeroVar(t(cecumMAGCounts)), ]
colonMAGCounts = colonMAGCounts[-nearZeroVar(t(colonMAGCounts)), ]
stoolMAGCounts = stoolMAGCounts[-nearZeroVar(t(stoolMAGCounts)), ]

### Check if we have the same MAGs in all three communities
all.equal(rownames(stoolMAGCounts), rownames(colonMAGCounts))
all.equal(rownames(cecumMAGCounts), rownames(colonMAGCounts))

Two_month_MAG <- calculate_p_values(cecumMAGCounts, colonMAGCounts, stool_counts)
Nine_month_MAG <- calculate_p_values(cecumMAGCounts, colonMAGCounts, stool_counts, age = "9")
Fifteen_month_MAG <- calculate_p_values(cecumMAGCounts, colonMAGCounts, stool_counts, age = "15")
TwentyFour_month_MAG <- calculate_p_values(cecumMAGCounts, colonMAGCounts, stool_counts, age = "24")
Thirty_month_MAG <- calculate_p_values(cecumMAGCounts, colonMAGCounts, stool_counts, age = "30")


write.csv(Two_month_MAG, "../data/mag_correlation_2M.csv")
write.csv(Nine_month_MAG, "../data/mag_correlation_9M.csv")
write.csv(Fifteen_month_MAG, "../data/mag_correlation_15M.csv")
write.csv(TwentyFour_month_MAG, "../data/mag_correlation_24M.csv")
write.csv(Thirty_month_MAG, "../data/mag_correlation_30M.csv")

cecumRxNCounts = as.matrix(read.csv("../data/cecum_rxn_abundance.csv", check.names = F, row.names = 1))
cecumRxNCounts <- apply(cecumRxNCounts, 2, function(col) {col/sum(col)})
cecumRxNCounts = cecumRxNCounts[-nearZeroVar(t(cecumRxNCounts)), ]

colonRxNCounts = as.matrix(read.csv("../data/colon_rxn_abundance.csv", check.names = F, row.names = 1))
colonRxNCounts <- apply(colonRxNCounts, 2, function(col) {col/sum(col)})
colonRxNCounts = colonRxNCounts[-nearZeroVar(t(colonRxNCounts)), ]

stoolRxNCounts = as.matrix(read.csv("../data/stool_rxn_abundance.csv", check.names = F, row.names = 1))
stoolRxNCounts <- apply(stoolRxNCounts, 2, function(col) {col/sum(col)})
stoolRxNCounts = stoolRxNCounts[-nearZeroVar(t(stoolRxNCounts)), ]

Two_month_RXN <- calculate_p_values(cecumRxNCounts, colonRxNCounts, stoolRxNCounts)
write.csv(Two_month_RXN, "../data/rxn_correlations_2M.csv")

Nine_month_RXN <- calculate_p_values(cecumRxNCounts, colonRxNCounts, stoolRxNCounts, "9")
write.csv(Nine_month_RXN, "../data/rxn_correlations_9M.csv")

Fifteen_month_RXN <- calculate_p_values(cecumRxNCounts, colonRxNCounts, stoolRxNCounts, "15")
write.csv(Fifteen_month_RXN, "../data/rxn_correlations_15M.csv")

TwentyFour_month_RXN <- calculate_p_values(cecumRxNCounts, colonRxNCounts, stoolRxNCounts, "24")
write.csv(TwentyFour_month_RXN, "../data/rxn_correlations_24M.csv")

Thirty_month_RXN <- calculate_p_values(cecumRxNCounts, colonRxNCounts, stoolRxNCounts, "30")
write.csv(Thirty_month_RXN, "../data/rxn_correlations_30M.csv")