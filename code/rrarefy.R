library("phyloseq")
cecum_mag_counts <- read.table("../data/absolute_cecum_mag_counts.csv", sep = ",", header = T, row.names = 1)
# REMOVE 30.91 becausre it only has one read
cecum_mag_counts = cecum_mag_counts %>% select(-Sum_30.91)
min(colSums(cecum_mag_counts))
colon_mag_counts = read.table("../data/absolute_colon_mag_counts.csv", sep = ",", header = T, row.names = 1)
min(colSums(colon_mag_counts))
stool_mag_counts = read.table("../data/absolute_stool_mag_counts.csv", sep = ",", header = T, row.names = 1)
min(colSums(stool_mag_counts))

rarefi = function(mag_counts, cutoff) {
  mag_counts = as.matrix(mag_counts)
  #print(colSums(mag_counts))
  # IS THE TRANSPOSE NECESSARY OR WRONG
  rarefied <- t( rrarefy(t(mag_counts), cutoff))
  #rarefied_diversity <- vegan::diversity(x = rarefied, index="shannon", MARGIN = 2)
  return(rarefi)
}

#
rare_cecum = rarefi(cecum_mag_counts, 789)
rare_colon = rarefi(colon_mag_counts, 1033)
rare_stool = rarefi(stool_mag_counts, 1512)

write.csv(rare_cecum, "../data/rarefied_cecum.csv")
write.csv(rare_colon, "../data/rarefied_colon.csv")
write.csv(rare_stool, "../data/rarefied_stool.csv")
