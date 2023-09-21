library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(caret)
library(lmerTest)
library(tidyverse)
source("utils.R")

cecum_counts = filter_per_organ(metadata, count_data, "Cecum")
colon_counts = filter_per_organ(metadata, count_data, "Colon")
stool_counts = filter_per_organ(metadata, count_data, "stool")

absolute_cecum_counts = filter_per_organ(metadata, count_data, "Cecum", normalize = F)
absolute_colon_counts = filter_per_organ(metadata, count_data, "Colon", normalize = F)
absolute_stool_counts = filter_per_organ(metadata, count_data, "stool", normalize = F)
write.csv(absolute_cecum_counts, file = "../data/cecum_mag_counts.csv", row.names = F)
write.csv(absolute_colon_counts, file = "../data/colon_mag_counts.csv", row.names = F)
write.csv(absolute_stool_counts, file = "../data/stool_mag_counts.csv", row.names = F)




calculate_p_values <- function(counts1,counts2, counts3,age = 2)
{
  counts1 = as.data.frame(counts1)
  counts2 = as.data.frame(counts2)
  counts3 = as.data.frame(counts3)
  commons = intersect(rownames(counts1), rownames(counts2))
  commons = intersect(commons, rownames(counts3))
  # Only using reactions that are present in all three sites
  print(age)
  print(length(commons))
  models =  list()
  ps =  list()
  data = data.frame(OTU = character(), coefficient = numeric(),  p.val = numeric()) 
  for(otu in unique((commons)))
  {
    d1 =  as.data.frame(t(counts1[otu, ]))
    d2 =  as.data.frame(t(counts2[otu, ]))
    d3 =  as.data.frame(t(counts3[otu, ]))
    current_otu = d1
    current_otu2 = merge(current_otu, d2   , by.x = 0, by.y = 0 )
    current_otu3 = merge(current_otu2, d3  , by.x = "Row.names", by.y = 0 )
    colnames(current_otu3) <- c("patID", "cecum_count", "colon_count", "stool_count")
    #rownames(current_otu3) <- gsub(".*/", "", current_otu3$patID)
    temp_meta = metadata %>% select(c(Age, ID))
    temp_meta$ID =  gsub("_.", "", temp_meta$ID)
    temp_meta = unique(temp_meta)
    temp_meta$fullID = apply(temp_meta, 1, function (row) { str_trim(paste0(row[[1]], ".", row[[2]]))})
    temp_meta$fullID = as.character(temp_meta$fullID)
    current_otu3$patID = as.character(current_otu3$patID)
    current_otu3$patID
    temp_meta$fullID
    intersect(current_otu3$patID, temp_meta$fullID)
    # MAGS USE DIFFERENT IDS
    current_otu3$patID   = gsub("Sum_", "", current_otu3$patID )
    current_otu3$patID   = gsub("/", ".", current_otu3$patID )
    current_otu4 = merge(current_otu3, temp_meta,  by.x = "patID", by.y = "fullID")
    if(nrow(current_otu4) != nrow(current_otu)) {
      stop()
    }
    current_otu4 = current_otu4[current_otu4$Age == age, ]
    rownames(current_otu4) = current_otu4$Row.names
    current_otu4 = current_otu4[-c(1)]
    colnames(current_otu4) <- c("1_","2_","3_","Age", "ID")
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
    out <- tryCatch(
    {
        model = lmer(Count ~ Organ +(1|ID) , data = long_data)
        model_sm = summary(model)
        models[[otu]] = model_sm
        ps[[otu]] = model_sm$coefficients[, c("Estimate", "Pr(>|t|)")][2, ]
        data[nrow(data)+1, ] <- c((otu),  ps[[otu]][[1]], ps[[otu]][[2]] )
    }, 
    error=function(e)
      {
      print("ERROR")
      print(otu)
      })
  }
  adjustNonNaNOTUS = p.adjust(data$p.val, method = "BH")
  data$p.adj = adjustNonNaNOTUS
  return(data)
}

#Two_month <- calculate_p_values(cecum_counts, colon_counts, stool_counts)

# Load 16S to MAG mapping
# 161
cecumMAGCounts <- cecum_counts
cecumMAGCounts = cecumMAGCounts %>% column_to_rownames("V2")
# 161
colonMAGCounts <- colon_counts
colonMAGCounts = colonMAGCounts %>% column_to_rownames("V2")

# 161
stoolMAGCounts <- stool_counts
stoolMAGCounts = stoolMAGCounts %>% column_to_rownames("V2")

cecumMAGCounts = cecumMAGCounts[-nearZeroVar(t(cecumMAGCounts)), ]
# 150
colonMAGCounts = colonMAGCounts[-nearZeroVar(t(colonMAGCounts)), ]
# 154
stoolMAGCounts = stoolMAGCounts[-nearZeroVar(t(stoolMAGCounts)), ]
# 154

### Check if we have the same MAGs in all three communities
all.equal(rownames(stoolMAGCounts), rownames(colonMAGCounts))
all.equal(rownames(cecumMAGCounts), rownames(colonMAGCounts))

write.csv(cecumMAGCounts, file = "../data/filtered_relative_cecum_mag_counts.csv")
write.csv(colonMAGCounts, file = "../data/filtered_relative_colon_mag_counts.csv")
write.csv(stoolMAGCounts, file = "../data/filtered_relative_stool_mag_counts.csv")

Two_month_MAG <- calculate_p_values(cecumMAGCounts, colonMAGCounts, stoolMAGCounts)
Nine_month_MAG <- calculate_p_values(cecumMAGCounts, colonMAGCounts, stoolMAGCounts, age = "9")
Fifteen_month_MAG <- calculate_p_values(cecumMAGCounts, colonMAGCounts, stoolMAGCounts, age = "15")
TwentyFour_month_MAG <- calculate_p_values(cecumMAGCounts, colonMAGCounts, stoolMAGCounts, age = "24")
Thirty_month_MAG <- calculate_p_values(cecumMAGCounts, colonMAGCounts, stoolMAGCounts, age = "30")


write.csv(Two_month_MAG, "../data/mag_correlation_2M.csv")
write.csv(Nine_month_MAG, "../data/mag_correlation_9M.csv")
write.csv(Fifteen_month_MAG, "../data/mag_correlation_15M.csv")
write.csv(TwentyFour_month_MAG, "../data/mag_correlation_24M.csv")
write.csv(Thirty_month_MAG, "../data/mag_correlation_30M.csv")

cecumRxNCounts = as.matrix(read.csv("../data/relative_cecum_rxn_abundance.csv", check.names = F, row.names = 1))
cecumRxNCounts <- apply(cecumRxNCounts, 2, function(col) {col/sum(col)})
cecumRxNCounts = cecumRxNCounts[-nearZeroVar(t(cecumRxNCounts)), ]

colonRxNCounts = as.matrix(read.csv("../data/relative_colon_rxn_abundance.csv", check.names = F, row.names = 1))
colonRxNCounts <- apply(colonRxNCounts, 2, function(col) {col/sum(col)})
colonRxNCounts = colonRxNCounts[-nearZeroVar(t(colonRxNCounts)), ]

stoolRxNCounts = as.matrix(read.csv("../data/relative_stool_rxn_abundance.csv", check.names = F, row.names = 1))
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
#

