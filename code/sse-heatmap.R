setwd("/home/td/Projects/Aging/Baseline/GIN/code")
source("utils.R")
# Load counts because we wanna plot them in the heatmap
filtered_rel_cecum = read.csv("../data/filtered_relative_cecum_rxn_abundance.csv", row.names = 1, check.names = F)
filtered_rel_colon = read.csv("../data/filtered_relative_colon_rxn_abundance.csv", row.names = 1, check.names = F)
filtered_rel_stool = read.csv("../data/filtered_relative_stool_rxn_abundance.csv", row.names = 1, check.names = F)


make_matrix = function(path_enrich_up, path_enrich_down) {
  # Load enricher results
  enricher_up = readRDS(path_enrich_up)@result
  enricher_down = readRDS(path_enrich_down)@result
  # We want to include up and downregulated sub systems because we want to see the difference in the plot
  enricher_all = rbind(enricher_up, enricher_down)
  enricher_all = enricher_all[enricher_all$p.adjust < 0.001, ]

  # TODO Does the geneID column only represent reaction that we "put into enricher?"
  # TODO Check if all reactions occour in our significant lists
  ss_to_rxn = enricher_all %>% select(Description, geneID)
  # Split the second column by "/"
  split_values <- strsplit(ss_to_rxn$geneID, "/")
  # Compute
  unique_reactions_across_ss =  unique(as.vector(flatten(split_values)))

  rxn_to_ss = readRDS("../data/df_rxn2subsys20230711MetaMouse.rds")

  present_rxn_to_ss = data.frame(matrix(ncol = 4, nrow = 0))
  colnames(present_rxn_to_ss) = colnames(rxn_to_ss)
  for(r in unique_reactions_across_ss) {
    row = rxn_to_ss[rxn_to_ss$ReactionID == r, ]
    present_rxn_to_ss = rbind(present_rxn_to_ss, row)
  }
  # Now add the count
  rxn_ss_count_accross_all_samples = data.frame(rxnid = as.character(), ss = as.character(), sss = as.character(), count_cecum = as.numeric(), count_colon = as.numeric(), count_stool = as.numeric())
  for(i in 1:nrow(present_rxn_to_ss)) {
    row = present_rxn_to_ss[i, ]
    rxn = row$ReactionID
    ss = row$Subsystem
    sss = row$SubsystemShort
    # TODO CHECK THIS WITH AN EXAMPLE
    counts = getTableForOtu(as.matrix(filtered_rel_cecum), as.matrix(filtered_rel_colon), as.matrix(filtered_rel_stool), rxn) %>%
      select(Organ, Count) %>%
      group_by(Organ) %>%
      summarise(mean(Count))

    aver_cecum_count_of_rxn_across_all_samples = counts %>% filter(Organ == "1_") %>% select('mean(Count)')
    aver_colon_count_of_rxn_across_all_samples = counts %>% filter(Organ == "2_") %>% select('mean(Count)')
    aver_stool_count_of_rxn_across_all_samples = counts %>% filter(Organ == "3_") %>% select('mean(Count)')
    rxn_ss_count_accross_all_samples[nrow(rxn_ss_count_accross_all_samples) + 1, ] = c(rxn, ss, sss,
                                                                                       aver_cecum_count_of_rxn_across_all_samples,
                                                                                       aver_colon_count_of_rxn_across_all_samples,
                                                                                       aver_stool_count_of_rxn_across_all_samples)
  }


  # We only want enrichee SSs
  # TODO FIGURE OUT WHY IS IS NOT |enricher_ALL| before we do th %in%
  rxn_ss_count_accross_all_samples =  rxn_ss_count_accross_all_samples[rxn_ss_count_accross_all_samples $ss %in% enricher_all$Description, ]


  averaged_bs_ss_rxn_ss_count_across_all_samples = rxn_ss_count_accross_all_samples %>%
    group_by(sss) %>%
    summarise(across(matches("count*."), mean, .names = "{.col}"))



  averaged_bs_ss_rxn_ss_count_across_all_samples = averaged_bs_ss_rxn_ss_count_across_all_samples %>% column_to_rownames("sss")
  return(averaged_bs_ss_rxn_ss_count_across_all_samples)
}

 {
 mat2M = make_matrix("../data/enricher_up_2M.RDS",  "../data/enricher_down_2M.RDS")
 pdf("../data/plots/sse_heatmap_2M.pdf")
 map2M =  gplots::heatmap.2(as.matrix(mat2M),
                           Colv = FALSE,
                           scale = "row",
                           trace = "none",
                           dendrogram = "none",
                            main = "2M",
                            ylab =  "Subsystems",
                            xlab =  "GI Site",
                            #offsetRow = -56,
                            key = F,
                            cexRow = 0.5,
                          col = c(colorRampPalette(RColorBrewer::brewer.pal("Reds", n = 9))(50))
 )
 dev.off()
 }

{
#mat9M2 = mat9M# t(scale(t(mat9M), center = TRUE, scale = TRUE))
#dendro <- hclust( dist(mat9M2))
#ddr <- as.dendrogram(dendro)
#ddr <- reorder(ddr, F)
#rowInd <- order.dendrogram(ddr)
#ord =   as.data.frame(dendro[c("labels")], rowInd)
#
#ord$t = map9M$rowInd
#
#mat9M2 = merge(mat9M2, ord, by.x = 0, by.y = "labels")
#
#colnames(mat9M2) = c("ss", "count_cecum", "count_colon", "count_stool", "level")
#
#o = mat9M2 %>% select(ss, level) %>% arrange(level)
#mat9M2$sss = factor(mat9M2$ss, levels = o$ss)
#
#m = as.data.frame(mat9M2) %>% pivot_longer(cols = starts_with("count_"),
#  names_to = "variable",         # Name of the new variable column
#  values_to = "value")
#
#p = ggplot(m, aes(x = variable, y = sss, fill = value)) +
#  geom_tile()
#p
#
#
#ggsave("../data/plots/sse_h2.pdf", p )

}

mat9M = make_matrix("../data/enricher_up_9M.RDS",  "../data/enricher_down_9M.RDS")
pdf("../data/plots/sse_heatmap_9M.pdf")
map9M =  gplots::heatmap.2(as.matrix(mat9M),
                           Colv = FALSE,
                           scale = "row",
                           trace = "none",
                           dendrogram = "none",
                           main = "9M",
                           ylab =  "Subsystems",
                           xlab =  "GI Site",
                           #offsetRow = -56,
                           key = F,
                           cexRow = 0.5,
                           cexCol = 0.5,
                           col = colorRampPalette(RColorBrewer::brewer.pal("Reds", n = 9))(50)
)
dev.off()

mat15M = make_matrix("../data/enricher_up_15M.RDS",  "../data/enricher_down_15M.RDS")
pdf("../data/plots/sse_heatmap_15M.pdf")
map15M =  gplots::heatmap.2(as.matrix(mat15M),
                           Colv = FALSE,
                           scale = "row",
                           trace = "none",
                           dendrogram = "none",
                           main = "15M",
                           ylab =  "Subsystems",
                           xlab =  "GI Site",
                           #offsetRow = -56,
                           key = F,
                           cexRow = 0.5,
                           cexCol = 0.5,
                           col = c(colorRampPalette(RColorBrewer::brewer.pal("Reds", n = 15))(50))
)
dev.off()
