library(dplyr)
library(ggplot2)

metadata <- read.csv("../data/Jena_mouse_clean_RNA.csv")
count_data <- read.csv("../data/otu_count_clean.csv", sep = ",", row.names = 1, check.names = F)

# Returns normalized count data for a specific organ
filter_per_organ = function (metadata, countdata, organ) {
  metadata_organ = metadata %>% filter(TissueID == organ)
  filtered_count_data = countdata[, metadata_organ$PatID]
  filtered_count_data <- apply(filtered_count_data, 2, function(col) {col/sum(col)})
  return(filtered_count_data)
}



cecum_counts = filter_per_organ(metadata, count_data, "Cecum")
cecum_counts = cecum_counts[-nearZeroVar(t(cecum_counts)), ]
colon_counts = filter_per_organ(metadata, count_data, "Colon")
colon_counts = colon_counts[-nearZeroVar(t(colon_counts)), ]
stool_counts = filter_per_organ(metadata, count_data, "stool")
stool_counts = stool_counts[-nearZeroVar(t(stool_counts)), ]

counts1 = colon_counts
calculate_p_values <- function(counts1,counts2, counts3, method = "lm")
{
  models =  list()
  ps =  list()
  if(method == "lm") 
  {
    data = data.frame(OTU = character(), intercept = numeric(), coeff = numeric(), p.val = numeric()) 
  }
  for(otu in unique((rownames(counts))))
  {
    current_otu = as.data.frame(counts[otu, ])
    rownames(current_otu) <- gsub(".*/", "", rownames(current_otu))
    temp_meta = metadata %>% select(c(Age, ID))
    current_otu = merge(current_otu, temp_meta,  by.x = 0, by.y = "ID")
    colnames(current_otu) <- c("ID", "count", "Age")
    model = lm(count ~ Age, data = current_otu)
    model_sm = summary(model)
    models[[otu]] = model_sm
    ps[[otu]] = model_sm$coefficients[, 4][[2]]
    data[nrow(data)+1, ] <- c((otu), coef(model), ps[[otu]] )
  }
  adjustNonNaNOTUS = p.adjust(ps, method = "BH")
  data$p.adj = adjustNonNaNOTUS
  return(data)
}

cecum_p_16S_values <- calculate_p_values(cecum_counts)

# Load 16S to MAG mapping
OTUtoMAG <- read.csv("../data/VsearchMap6Out99.tsv", sep = "\t", header = F, check.names = F)
cecumMAGCounts <- cecum_counts[rownames(cecum_counts) %in% OTUtoMAG$V1, ]
colonMAGCounts <- colon_counts[rownames(colon_counts) %in% OTUtoMAG$V1, ]
stoolMAGCounts <- stool_counts[rownames(stool_counts) %in% OTUtoMAG$V1, ]

### Check if we have the same MAGs in all three communities
all.equal(rownames(stoolMAGCounts), rownames(colonMAGCounts))
all.equal(rownames(cecumMAGCounts), rownames(colonMAGCounts))

cecum_MAG_p_values <- calculate_p_values(cecumMAGCounts)
colon_MAG_p_values <- calculate_p_values(colonMAGCounts)
stool_MAG_p_values <- calculate_p_values(stoolMAGCounts)

current_otu = as.data.frame(cecum_counts["Zotu3", ])
rownames(current_otu) = gsub(".*/", "", rownames(current_otu))
current_otu = merge(current_otu, metadata %>% select(Age, ID), by.x = 0, by.y = "ID")
colnames(current_otu) <- c("ID", "count", "Age")
ggplot(current_otu, aes(x = Age, y = count)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Scatter Plot with Linear Regression Lines by Group", x = "Age", y = "OTU count") 
#  