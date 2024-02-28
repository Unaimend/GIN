#setwd("/home/td/Projects/Aging/Baseline/GIN/code")
source("utils.R")
# Unnnormalized, unfiltered (by HQness) ABSOLUTE MAG counts
cecum_counts_inc_mq  <-  filter_per_organ(metadata, count_data, "Cecum", filter_HQ_only = FALSE, normalize = FALSE)
colon_counts_inc_mq  <-  filter_per_organ(metadata, count_data, "Colon", filter_HQ_only = FALSE, normalize = FALSE)
stool_counts_inc_mq  <-  filter_per_organ(metadata, count_data, "stool", filter_HQ_only = FALSE, normalize = FALSE)


# We dont not have all samples for all three organs thus we take the intersection
samples_in_all_three <- intersect(colnames(cecum_counts_inc_mq), colnames(colon_counts_inc_mq))
samples_in_all_three <- intersect(samples_in_all_three, colnames(stool_counts_inc_mq))
cecum_counts_inc_mq_intersected <- cecum_counts_inc_mq[samples_in_all_three]
colon_counts_inc_mq_intersected <- colon_counts_inc_mq[samples_in_all_three]
stool_counts_inc_mq_intersected <- stool_counts_inc_mq[samples_in_all_three]

colnames(cecum_counts_inc_mq_intersected) <- gsub("Sum_", "",   colnames(cecum_counts_inc_mq_intersected))
colnames(colon_counts_inc_mq_intersected) <- gsub("Sum_", "",   colnames(colon_counts_inc_mq_intersected))
colnames(stool_counts_inc_mq_intersected) <- gsub("Sum_", "",   colnames(stool_counts_inc_mq_intersected))

nrow(stool_counts_inc_mq_intersected) == 161
nrow(colon_counts_inc_mq_intersected) == 161
nrow(cecum_counts_inc_mq_intersected) == 161

# Those now also contain the mq mags
# Used in building the communities
write.table(cecum_counts_inc_mq_intersected, sep = "\t", file = "../data/communities/cecum_counts_inc_mq_intersected.tsv")
write.table(colon_counts_inc_mq_intersected, sep = "\t", file = "../data/communities/colon_counts_inc_mq_intersected.tsv")
write.table(stool_counts_inc_mq_intersected, sep = "\t", file = "../data/communities/stool_counts_inc_mq_intersected.tsv")
