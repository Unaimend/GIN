library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(caret)
library(lmerTest)
library(tidyverse)
setwd("/home/td/Projects/Aging/Baseline/GIN/code")
source("utils.R")

cecum_counts = filter_per_organ(metadata, count_data, "Cecum")
colon_counts = filter_per_organ(metadata, count_data, "Colon")
stool_counts = filter_per_organ(metadata, count_data, "stool")

absolute_cecum_counts = filter_per_organ(metadata, count_data, "Cecum", normalize = F)
absolute_colon_counts = filter_per_organ(metadata, count_data, "Colon", normalize = F)
absolute_stool_counts = filter_per_organ(metadata, count_data, "stool", normalize = F)

# We dont not have all samples for all three organs thus we take the intersection
samples_in_all_three = intersect(colnames(absolute_cecum_counts), colnames(absolute_colon_counts))
samples_in_all_three = intersect(samples_in_all_three, colnames(absolute_stool_counts))
absolute_cecum_counts_intersected = absolute_cecum_counts[samples_in_all_three]
absolute_colon_counts_intersected = absolute_colon_counts[samples_in_all_three]
absolute_stool_counts_intersected = absolute_stool_counts[samples_in_all_three]

# Remove Sub_ so the names are again the same as in the gapseqs
colnames(absolute_cecum_counts_intersected) =  gsub("Sum_", "",   colnames(absolute_cecum_counts_intersected))
colnames(absolute_colon_counts_intersected) =  gsub("Sum_", "",   colnames(absolute_colon_counts_intersected))
colnames(absolute_stool_counts_intersected) =  gsub("Sum_", "",   colnames(absolute_stool_counts_intersected))


colnames(absolute_cecum_counts_intersected) =  gsub("/", "_",   colnames(absolute_cecum_counts_intersected))
colnames(absolute_colon_counts_intersected) =  gsub("/", "_",   colnames(absolute_colon_counts_intersected))
colnames(absolute_stool_counts_intersected) =  gsub("/", "_",   colnames(absolute_stool_counts_intersected))


write.table(absolute_cecum_counts_intersected, file = "../MicrobiomeGS/data/absolute_cecum_counts_intersected.tsv", sep = '\t')
write.table(absolute_colon_counts_intersected, file = "../MicrobiomeGS/data/absolute_colon_counts_intersected.tsv", sep = '\t')
write.table(absolute_stool_counts_intersected, file = "../MicrobiomeGS/data/absolute_stool_counts_intersected.tsv", sep = '\t')

filt_cecum_counts = filter_per_organ(metadata, count_data, "Cecum", zero_var_filt = T)
filt_colon_counts = filter_per_organ(metadata, count_data, "Colon", zero_var_filt = T)
filt_stool_counts = filter_per_organ(metadata, count_data, "stool", zero_var_filt = T)

write.csv(absolute_cecum_counts, file = "../data/absolute_cecum_mag_counts.csv", row.names = T)
write.csv(absolute_colon_counts, file = "../data/absolute_colon_mag_counts.csv", row.names = T)
write.csv(absolute_stool_counts, file = "../data/absolute_stool_mag_counts.csv", row.names = T)

write.csv(filt_cecum_counts, file = "../data/filtered_relative_cecum_mag_counts.csv")
write.csv(filt_colon_counts, file = "../data/filtered_relative_colon_mag_counts.csv")
write.csv(filt_stool_counts, file = "../data/filtered_relative_stool_mag_counts.csv")




#Two_month <- calculate_p_values(cecum_counts, colon_counts, stool_counts)

# Load 16S to MAG mapping
# 124
cecumMAGCounts <- cecum_counts
# 124
colonMAGCounts <- colon_counts

# 124
stoolMAGCounts <- stool_counts

cecumMAGCounts = read.csv(file = "../data/filtered_relative_cecum_mag_counts.csv", row.names = 1)
# 118
colonMAGCounts = read.csv(file = "../data/filtered_relative_colon_mag_counts.csv", row.names = 1)
# 120
stoolMAGCounts = read.csv(file = "../data/filtered_relative_stool_mag_counts.csv", row.names = 1)
# 119

### Check if we have the same MAGs in all three communities
all.equal(rownames(stoolMAGCounts), rownames(colonMAGCounts))
all.equal(rownames(cecumMAGCounts), rownames(colonMAGCounts))





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

