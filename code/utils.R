
library(lmerTest)
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

getTableForOtu <- function(counts1, counts2, counts3, otu, age = 2)
{
  current_otu = as.data.frame(counts1[otu, ])
  current_otu2 = merge(current_otu,  as.data.frame(counts2[otu, ]) , by.x = 0, by.y = 0 )
  current_otu3 = merge(current_otu2,  as.data.frame(counts3[otu, ]) , by.x = "Row.names", by.y = 0 )
  colnames(current_otu3) <- c("patID", "cecum_count", "colon_count", "stool_count")
  rownames(current_otu3) <- gsub(".*/", "", current_otu3$patID)
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
  return(long_data)
}


calculate_p_values <- function(counts1, counts2, counts3, age = 2, fluxes = FALSE) {
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
    current_otu3 <- merge(current_otu2, d3  , by.x = "Row.names", by.y = 0 )
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
    if(fluxes) {
      current_otu3$patID   = gsub(".RDS", "", current_otu3$patID )
      current_otu3$patID   = gsub("_", ".", current_otu3$patID )
    }
    current_otu4 = merge(current_otu3, temp_meta,  by.x = "patID", by.y = "fullID")
    # Ignore this for now But we lost a sample
    # if(nrow(current_otu4) != nrow(current_otu)) {
    #   stop()
    # }
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
    long_data$Count = as.numeric(long_data$Count)
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
      print(e)
      })
  }
  adjustNonNaNOTUS = p.adjust(data$p.val, method = "BH")
  data$p.adj = adjustNonNaNOTUS
  print("TER")
  return(data)
}
