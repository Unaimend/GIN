cecum_mag_counts <- read.table("../data/absolute_cecum_mag_counts.csv", sep = ",", header = T, row.names = 1)
colon_mag_counts = read.table("../data/absolute_colon_mag_counts.csv", sep = ",", header = T, row.names = 1)
stool_mag_counts = read.table("../data/absolute_stool_mag_counts.csv", sep = ",", header = T, row.names = 1)
taxonomy = read.csv("../data/TaxonomyHQDraftMAGs.csv", sep = "\t")
taxonomy = taxonomy[c(1:8)]
rownames(taxonomy) = taxonomy$user_genome
taxonomy = taxonomy %>% select(-user_genome)
rownames(taxonomy) = paste0("OTU_", rownames(taxonomy))
rownames(cecum_mag_counts) = paste0("OTU_", rownames(cecum_mag_counts))
colnames(cecum_mag_counts) = paste0("Sample_", colnames(cecum_mag_counts))
taxonomy = taxonomy[rownames(taxonomy) %in% rownames(cecum_mag_counts), ]
cecum_mag_counts = cecum_mag_counts[rownames(cecum_mag_counts) %in% rownames(taxonomy), ]
taxonomy = as.data.frame(taxonomy) %>% rownames_to_column(var = "Name") %>%  # Convert row names to a column
  arrange(Name) %>%
  column_to_rownames("Name")
taxonomy = as.matrix(taxonomy)
cecum_mag_counts = as.data.frame(cecum_mag_counts) %>% rownames_to_column(var = "Name") %>%  # Convert row names to a column
  arrange(Name) %>%
  column_to_rownames("Name")

cecum_mag_counts = as.matrix(cecum_mag_counts)

TAX = tax_table(taxonomy)
OTU = otu_table(cecum_mag_counts, taxa_are_rows = T)

all.equal(rownames(taxonomy), rownames(cecum_mag_counts))
all.equal(rownames(TAX), rownames(OTU))


library("phyloseq")
p = phyloseq(OTU, TAX)
p2 = rarefy_even_depth(p, rngseed = T, rngseed = 711)

alpha_diversity <- estimate_richness(p, measures = c("observed", "shannon", "simpson"))
alpha_diversity2 <- estimate_richness(p2, measures = c("observed", "shannon", "simpson"))
alpha_df <- data.frame(sample_names = rownames(alpha_diversity), alpha_diversity)