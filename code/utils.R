
library(dplyr)
library(tidyverse)
metadata <- read.csv("../data/Jena_mouse_clean_RNA.csv")
count_data <- read.csv("../data/otu_count_clean.csv", sep = ",", row.names = 1, check.names = F)
OTUtoMAG <- read.csv("../data/VsearchMap6Out99.tsv", sep = "\t", header = F, check.names = F) %>% select(V1, V2)
OTUtoMAG$V2 = gsub("::.*", "", OTUtoMAG$V2)

# Returns normalized count data for a specific organ
filter_per_organ = function (metadata, countdata, organ, normalize = T, map_to_mags = TRUE, zero_var_filt = F, filter_HQ_only = T) {
  metadata_organ = metadata %>% filter(TissueID == organ)
  filtered_count_data = countdata[, metadata_organ$PatID]
  # Normalize data
  colnames(filtered_count_data) = gsub("_.", "", colnames(filtered_count_data))
  if(map_to_mags) {
    result = merge(filtered_count_data, OTUtoMAG, by.x = 0, by.y = "V1")
    result2 <- result %>% group_by(V2) %>% summarize(across(matches(".*/"), sum, .names = "Sum_{.col}"))
    result3 = as.data.frame(result2) %>% column_to_rownames("V2")
  } else
  {
    result3 = filtered_count_data
  }
  if(filter_HQ_only) {
    taxonomy  = read.csv("../data/TaxonomyHQDraftMAGs.csv", sep = "\t") %>% select(user_genome)
    filter_hq = rownames(result3) %in% as.vector(taxonomy$user_genome)
    result3 = result3[filter_hq, ]

  }
  if(zero_var_filt) {
    result3 = result3[-nearZeroVar(t(result3)), ]
  }
  if(normalize)
  {
    result3 <- apply(result3, 2, function(col) {col/sum(col)})
  }
  return(as.data.frame(result3))
}

load_and_calc = function(filepath, pattern, grp) {
  if(is.data.frame(filepath)) {
    cecum_mag_month_counts = filepath %>%
      tibble::rownames_to_column(var = "rowname") %>%
      filter(grepl(pattern, rowname)) %>%
      tibble::column_to_rownames(var = "rowname") 
  } else {
    cecum_mag_month_counts = as.data.frame(t(read.table(filepath, sep = ",", header = T, row.names = 1)))  %>%
      tibble::rownames_to_column(var = "rowname") %>%
      filter(grepl(pattern, rowname)) %>%
      tibble::column_to_rownames(var = "rowname") 
  }
  mag_shannon <- as.data.frame(diversity(cecum_mag_month_counts, index = "shannon"))        
  mag_shannon$grp = grp
  colnames(mag_shannon) = c("alpha", "grp")
  return(mag_shannon)
}
