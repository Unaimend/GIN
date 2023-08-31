
library(dplyr)
metadata <- read.csv("../data/Jena_mouse_clean_RNA.csv")
count_data <- read.csv("../data/otu_count_clean.csv", sep = ",", row.names = 1, check.names = F)
OTUtoMAG <- read.csv("../data/VsearchMap6Out99.tsv", sep = "\t", header = F, check.names = F) %>% select(V1, V2)
OTUtoMAG$V2 = gsub("::.*", "", OTUtoMAG$V2)

# Returns normalized count data for a specific organ
filter_per_organ = function (metadata, countdata, organ, normalize = T, map_to_mags = T) {
  metadata_organ = metadata %>% filter(TissueID == organ)
  filtered_count_data = countdata[, metadata_organ$PatID]
  # Normalize data
  if(normalize)
  {
    filtered_count_data <- apply(filtered_count_data, 2, function(col) {col/sum(col)})
  }
  colnames(filtered_count_data) = gsub("_.", "", colnames(filtered_count_data))
  if(map_to_mags) {
    result = merge(filtered_count_data, OTUtoMAG, by.x = 0, by.y = "V1")
    result2 <- as.matrix(result%>%
                           group_by(Row.names) %>%
                           summarize(across(matches(".*/"), sum, .names = "Sum_{.col}")) %>%
                           column_to_rownames("Row.names") 
    )
  } else
  {
    result2 = filtered_count_data
  }
  
  return(result2)
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
