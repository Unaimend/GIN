setwd("/home/td/Projects/Aging/Baseline/GIN/code")
library(vegan)
library(patchwork)
library(ggplot2)
source("utils.R")
library(ggpubr)
library(rstatix)

# MAP TO MAGS AND  AND CHECK IF WE CAN USE RELATIVE COUNTS, AND FILTER DOUBLES
cecum_mag_counts <- t(read.table("../data/cecum_mag_counts.csv", sep = ",", header = T, row.names = 1))
colon_mag_counts = t(read.table("../data/colon_mag_counts.csv", sep = ",", header = T, row.names = 1))
stool_mag_counts = t(read.table("../data/stool_mag_counts.csv", sep = ",", header = T, row.names = 1))


mag= function(cecum_counts, colon_counts, stool_counts,  month,
              show_label = F, kind = "MAG", pattern, lim = c(1,4.5),
              breaks = c(1,2,3,4), step = 0) {
  mag_shannon <- as.data.frame(diversity(cecum_counts, index = "shannon"))                          
  mag_shannon$grp = "Cecum"
  colnames(mag_shannon) = c("alpha", "grp")
  mag_shannon = mag_shannon[grep(pattern, rownames(mag_shannon)), ]
  
  mag_shannon2 <- as.data.frame(diversity(colon_counts, index = "shannon"))                          
  mag_shannon2$grp = "Proximal C."
  colnames(mag_shannon2) = c("alpha", "grp")
  mag_shannon2 = mag_shannon2[grep(pattern, rownames(mag_shannon2)), ]
  
  mag_shannon3 <- as.data.frame(diversity(stool_counts, index = "shannon"))                          
  mag_shannon3$grp = "Distal C."
  colnames(mag_shannon3) = c("alpha", "grp")
  mag_shannon3 = mag_shannon3[grep(pattern, rownames(mag_shannon3)), ]
  
  overall_counts = data.frame("alpha" = as.numeric(), "grp" = as.character())
  overall_counts = rbind(overall_counts, mag_shannon)
  overall_counts = rbind(overall_counts, mag_shannon2)
  overall_counts = rbind(overall_counts, mag_shannon3)
  overall_counts$grp = factor(overall_counts$grp, levels = c("Cecum", "Proximal C.", "Distal C."))
  res.kruskal <- overall_counts %>% kruskal_test(alpha ~ grp)
  pwc <- overall_counts %>%
    dunn_test(alpha ~ grp, p.adjust.method = "BH")
  pwc <- pwc %>% add_xy_position(x = "grp")
  p<-ggplot(overall_counts, aes(y=alpha, x= grp, color = grp)) +
    geom_boxplot() +
    ggtitle(paste0(month, " months")) +
    labs(x = "GI Site") +
    geom_jitter(color="black", size=0.4, alpha=0.9, width = 0.1) +
    ylab(if(show_label) {"Shannon diversity"}) +
    theme(plot.title = element_text(size = 20), axis.text=element_text(size=15)) +
    coord_cartesian(ylim = lim)  +
    scale_y_continuous(breaks = breaks)  +
    guides(color = "none") +
    labs(
      subtitle = get_test_label(res.kruskal, detailed = TRUE)

    ) +
    stat_pvalue_manual(
    pwc, label = "Dunns, p.adj = {p.adj.signif}",
    vjust = 0, bracket.nudge.y = 0.0, step.increase = step
  )
  return(p)
}
{
  p2 = mag(cecum_mag_counts, colon_mag_counts, stool_mag_counts,  "2",   show_label = T,   kind = "MAG","_2\\.", step = 0.03)
  p9 = mag(cecum_mag_counts, colon_mag_counts, stool_mag_counts,  "9",   show_label = F,   kind = "MAG","_9\\.", step = 0.03)
  p15 = mag(cecum_mag_counts, colon_mag_counts, stool_mag_counts,  "15", show_label = F, kind = "MAG","_15\\.", step = 0.02)
  p24 = mag(cecum_mag_counts, colon_mag_counts, stool_mag_counts,  "24", show_label = T, kind = "MAG","_24\\.", step = 0.1)
  p30 = mag(cecum_mag_counts, colon_mag_counts, stool_mag_counts,  "30", show_label = F, kind = "MAG","_30\\.", step = 0.1)

  arranged_plots <- wrap_plots(
    p2 + p9 + p15,
    p24 + p30,
    ncol = 1
  )
  final = arranged_plots
  ggsave(final, filename = "../data/plots/mag_alpha_div_along_sites.pdf", width = 12, height  = 12)
  }
{
cecum_rxn_counts = t(as.matrix(read.csv("../data/cecum_rxn_abundance.csv", check.names = F, row.names = 1)))
colon_rxn_counts = t(as.matrix(read.csv("../data/colon_rxn_abundance.csv", check.names = F, row.names = 1)))
stool_rxn_counts = t(as.matrix(read.csv("../data/stool_rxn_abundance.csv", check.names = F, row.names = 1)))

r2  = mag(cecum_rxn_counts, colon_rxn_counts, stool_rxn_counts,  "2",     show_label = T, kind = "Reaction","2/" , lim = c(6.9, 7.4), breaks = waiver(), step = 0.05)
r9  = mag(cecum_rxn_counts, colon_rxn_counts, stool_rxn_counts,  "9",     show_label = F, kind = "Reaction","9/" , lim = c(6.9, 7.4), breaks = waiver(), step = 0.09)
r15 = mag(cecum_rxn_counts, colon_rxn_counts, stool_rxn_counts,  "15",    show_label = F, kind = "Reaction","15/", lim = c(6.9, 7.4), breaks = waiver(), step = 0.09)
r24 = mag(cecum_rxn_counts, colon_rxn_counts, stool_rxn_counts,  "24",    show_label = T, kind = "Reaction","24/", lim = c(6.8, 7.4), breaks = waiver(), step = 0.07)
r30 = mag(cecum_rxn_counts, colon_rxn_counts, stool_rxn_counts,  "30",    show_label = F, kind = "Reaction","30/", lim = c(6.8, 7.4), breaks = waiver(), step = 0.07)
arranged_plots_rxn <- wrap_plots(
  r2 + r9 + r15,
  r24 + r30,
  ncol = 1
)
final_rxn = arranged_plots_rxn
ggsave(final_rxn, filename = "../data/plots/rxn_alpha_div_along_sites.pdf", width = 12, height  = 12)
}


taxonomy = read.csv("../data/TaxonomyHQDraftMAGs.csv", sep = "\t")

cecum_counts = filter_per_organ(metadata, count_data, "Cecum")
colon_counts = filter_per_organ(metadata, count_data, "Colon")
stool_counts = filter_per_organ(metadata, count_data, "stool")

cecum_count_mag_an = merge(OTUtoMAG, cecum_counts, by.x = "V1", by.y = 0)


cecum_count_mag_an = cecum_count_mag_an %>% select(-V1)
length(cecum_count_mag_an$V2)
length(unique(cecum_count_mag_an$V2))
length(cecum_count_mag_an$V2) == length(unique(cecum_count_mag_an$V2))