setwd("/home/td/Projects/Aging/Baseline/GIN/code")
library(vegan)
library(patchwork)
library(ggplot2)
source("utils.R")
library(ggpubr)
library(rstatix)
library(dplyr)

# MAP TO MAGS AND  AND CHECK IF WE CAN USE RELATIVE COUNTS, AND FILTER DOUBLES
cecum_mag_counts <- t(read.table("../data/rarefied_cecum.csv", sep = ",", header = T, row.names = 1))
colon_mag_counts = t(read.table("../data/rarefied_colon.csv", sep = ",", header = T, row.names = 1))
stool_mag_counts = t(read.table("../data/rarefied_stool.csv", sep = ",", header = T, row.names = 1))
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
#{
#  p2 = mag(cecum_mag_counts, colon_mag_counts, stool_mag_counts,  "2",   show_label = T,   kind = "MAG","_2\\.", step = 0.03)
#  p9 = mag(cecum_mag_counts, colon_mag_counts, stool_mag_counts,  "9",   show_label = F,   kind = "MAG","_9\\.", step = 0.03)
#  p15 = mag(cecum_mag_counts, colon_mag_counts, stool_mag_counts,  "15", show_label = F, kind = "MAG","_15\\.", step = 0.02)
#  p24 = mag(cecum_mag_counts, colon_mag_counts, stool_mag_counts,  "24", show_label = T, kind = "MAG","_24\\.", step = 0.1)
#  p30 = mag(cecum_mag_counts, colon_mag_counts, stool_mag_counts,  "30", show_label = F, kind = "MAG","_30\\.", step = 0.1)
#
#  arranged_plots <- wrap_plots(
#    p2 + p9 + p15,
#    p24 + p30,
#    ncol = 1
#  )
#  final = arranged_plots
#  ggsave(final, filename = "../data/plots/mag_alpha_div_along_sites.pdf", width = 12, height  = 12)
#  }
cecum_rxn_counts = t(as.matrix(read.csv("../data/absolute_cecum_rxn_abundance.csv", check.names = F, row.names = 1)))
colon_rxn_counts = t(as.matrix(read.csv("../data/absolute_colon_rxn_abundance.csv", check.names = F, row.names = 1)))
stool_rxn_counts = t(as.matrix(read.csv("../data/absolute_stool_rxn_abundance.csv", check.names = F, row.names = 1)))

r2  = mag(cecum_rxn_counts, colon_rxn_counts, stool_rxn_counts,  "2",     show_label = T, kind = "Reaction","2." , lim = c(6.9, 7.4), breaks = waiver(), step = 0.05)
r9  = mag(cecum_rxn_counts, colon_rxn_counts, stool_rxn_counts,  "9",     show_label = F, kind = "Reaction","9." , lim = c(6.9, 7.4), breaks = waiver(), step = 0.09)
r15 = mag(cecum_rxn_counts, colon_rxn_counts, stool_rxn_counts,  "15",    show_label = F, kind = "Reaction","15.", lim = c(6.9, 7.4), breaks = waiver(), step = 0.09)
r24 = mag(cecum_rxn_counts, colon_rxn_counts, stool_rxn_counts,  "24",    show_label = T, kind = "Reaction","24.", lim = c(6.8, 7.4), breaks = waiver(), step = 0.07)
r30 = mag(cecum_rxn_counts, colon_rxn_counts, stool_rxn_counts,  "30",    show_label = F, kind = "Reaction","30.", lim = c(6.8, 7.4), breaks = waiver(), step = 0.07)
arranged_plots_rxn <- wrap_plots(
  r2 + r9 + r15,
  r24 + r30,
  ncol = 1
)
final_rxn = arranged_plots_rxn
ggsave(final_rxn, filename = "../data/plots/rxn_alpha_div_along_sites.pdf", width = 12, height  = 12)


cecum_counts = filter_per_organ(metadata, count_data, "Cecum")
colon_counts = filter_per_organ(metadata, count_data, "Colon")
stool_counts = filter_per_organ(metadata, count_data, "stool")





rel_tax = function(cecum, colon, stool, age, pattern, color_desc = F, show_label = F) {
  taxonomy  = read.csv("../data/TaxonomyHQDraftMAGs.csv", sep = "\t")
  phylum = taxonomy %>% select(user_genome, Phylum )
  cecum = cecum %>% select(matches(pattern))
  cecum = merge(cecum, phylum, by.x = 0, by.y = "user_genome")
  cecum <- cecum %>% group_by(Phylum) %>% summarize(across(matches(".*/"), sum, .names = "Sum_{.col}"))
  cecum = cecum %>%  pivot_longer(cols = starts_with("Sum_"),    # Specify the columns to pivot
                                  names_to = "variable",         # Name of the new variable column
                                  values_to = "value")

  bac_order = cecum %>% filter(Phylum == "Bacteroidota") %>% arrange(value)
  cecum$organ = "Cecum"

  colon = colon %>% select(matches(pattern))
  colon = merge(colon, phylum, by.x = 0, by.y = "user_genome")
  colon <- colon %>% group_by(Phylum) %>% summarize(across(matches(".*/"), sum, .names = "Sum_{.col}"))
  colon = colon %>%  pivot_longer(cols = starts_with("Sum_"),    # Specify the columns to pivot
                                  names_to = "variable",         # Name of the new variable column
                                  values_to = "value")
  colon$organ = "Proximal C."


  stool = stool %>% select(matches(pattern))
  stool = merge(stool, phylum, by.x = 0, by.y = "user_genome")
  stool <- stool %>% group_by(Phylum) %>% summarize(across(matches(".*/"), sum, .names = "Sum_{.col}"))
  stool = stool %>%  pivot_longer(cols = starts_with("Sum_"),    # Specify the columns to pivot
                                  names_to = "variable",         # Name of the new variable column
                                  values_to = "value")
  stool$organ = "Distal C."

  overall = data.frame(Phylum = as.character(), variable = as.character(), value = as.character(), organ = as.character())
  overall = overall %>% rbind(cecum) %>% rbind(colon) %>% rbind(stool)
  overall$organ = factor(overall$organ, levels = c("Cecum", "Proximal C.", "Distal C."))
  all_samples = unique(overall$variable)

  missing_in_cecum = setdiff(all_samples, bac_order$variable)
  overall$variable = gsub("Sum_Sum_", "", overall$variable)

  missing_in_cecum =  gsub("Sum_Sum_", "", missing_in_cecum)
  bac_order$variable = gsub("Sum_Sum_", "", bac_order$variable)

  overall$variable = factor(overall$variable, levels = c(bac_order$variable, missing_in_cecum))
  overall = overall %>% filter(variable != "30/91")
  overall = overall %>% filter(variable != "30/92")
  p1 <- overall %>%
    ggplot(aes(x = variable, y = value)) +
    geom_bar(aes(fill = Phylum), stat = "identity", position = "fill") +
    facet_grid(~organ, scales = "free_x") +
    ylab(if(show_label) { "Relative Abundance"}) +
    xlab("Sample") +
    ggtitle(paste0(age, " month")) +
    if(color_desc) {  theme(axis.text.x = element_text(angle = 90))} else { theme(axis.text.x = element_text(angle = 90), legend.position = "none") }
  return(p1)
}

two_month = rel_tax(cecum_counts, colon_counts, stool_counts, "2", "2/", show_label = T)
nine_month = rel_tax(cecum_counts, colon_counts, stool_counts, "9", "9/")
fifteen_month = rel_tax(cecum_counts, colon_counts, stool_counts, "15", "15/", color_desc = T)
twentyfour_month = rel_tax(cecum_counts, colon_counts, stool_counts, "24", "24/", color_desc = F)
thirty_month = rel_tax(cecum_counts, colon_counts, stool_counts, "30", "30/", color_desc = T)
final_tax = wrap_plots( two_month + nine_month + fifteen_month,
                             twentyfour_month + thirty_month,
                        ncol = 1)
ggsave(final_tax, filename = "../data/plots/rel_tax.pdf", width = 12, height = 12)




t = readRDS("../data/df_rxn2subsys20230711MetaMouse.rds")