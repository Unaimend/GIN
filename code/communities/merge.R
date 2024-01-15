# Merge single MQ gapseqs and HQ gapseqs into a single list to be used in the comm. script
library(sybil)
setwd("/home/td/Projects/Aging/Baseline/GIN/code")
hq <- readRDS("../data/metamouse_HQBins-20230711.RDS")
# Set the path to the folder containing the RDS files
folder_path <- "../data/communities/mq_models"

# List all files in the folder with the extension *.RDS
rds_files <- list.files(path = folder_path, pattern = "\\.RDS$", full.names = TRUE)

# Initialize an empty list to store the loaded RDS objects
r <- list()

# Loop through each file and load the RDS object into the list
for (file in rds_files) {
  r[[file]] <- readRDS(file)
}

mqs  <- list()

i  <-  1
for (model in r) {
  mqs[[model@mod_name]] <- r[[i]]
  i  <-  i + 1
}
# The should be no intersection
intersect(names(mqs), names(hq))

mq_and_hq  <-  c(hq, mqs)
saveRDS(mq_and_hq, file = "../data/metamouse_HQBins-20230711_MQBIns-20240109.RDS")
