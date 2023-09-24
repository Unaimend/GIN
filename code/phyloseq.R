library("phyloseq")
cecum_mag_counts <- read.table("../data/absolute_cecum_mag_counts.csv", sep = ",", header = T, row.names = 1)
# REMOVE 30.91 becausre it only has one read
cecum_mag_counts = cecum_mag_counts %>% select(-Sum_30.91)
min(colSums(cecum_mag_counts))
colon_mag_counts = read.table("../data/absolute_colon_mag_counts.csv", sep = ",", header = T, row.names = 1)
min(colSums(colon_mag_counts))
i = intersect(colnames(cecum_mag_counts), colnames(colon_mag_counts))
cecum_mag_counts = cecum_mag_counts[i]
colon_mag_counts = colon_mag_counts[i]
all.equal(colnames(cecum_mag_counts), colnames(colon_mag_counts))
cecum_mag_counts = rbind(cecum_mag_counts, colon_mag_counts)
#stool_mag_counts = read.table("../data/absolute_stool_mag_counts.csv", sep = ",", header = T, row.names = 1)

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

phylo = phyloseq(OTU, TAX)

data_phylo_filt = filter_taxa(phylo, function(x) sum(x > 2) > (0.11 * length(x)), TRUE)
set.seed(1782) # set seed for analysis reproducibility
OTU_filt_rar = rarefy_even_depth(otu_table(data_phylo_filt), rngseed = TRUE, replace = FALSE) # rarefy the raw data using Phyloseq package
data_otu_filt_rar = data.frame(otu_table(OTU_filt_rar))
data_phylo_filt_rar <- phyloseq(OTU_filt_rar, TAX)
dist_bc <- as.matrix(vegdist(data_otu_filt_rar, method = "bray"))
pcoa_bc = ordinate(data_phylo_filt_rar, "PCoA", "bray")

plot_ordination(data_phylo_filt_rar, pcoa_bc, color = "site") +
  geom_point(size = 3)
