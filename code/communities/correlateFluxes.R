setwd("/home/td/Projects/Aging/Baseline/GIN/code")
library(dplyr)
library(data.table)
library(janitor)
source("utils.R")
library(caret)
# Load cecal fluxes

calculated_fluxes <- function(cecal, proximal, distal, age = 2) {
  cecal_flux_with_host <- t(read.table(cecal,
                                       sep = ",", head = TRUE)) %>%  row_to_names(1)
  proximal_flux_with_host <- t(read.table(proximal,
                                          sep = ",", head = TRUE)) %>%  row_to_names(1)
  distal_flux_with_host <- t(read.table(distal,
                                        sep = ",", head = TRUE)) %>%  row_to_names(1)

  cecal_flux_with_host[is.na(cecal_flux_with_host)] <- 0
  proximal_flux_with_host[is.na(proximal_flux_with_host)] <- 0
  distal_flux_with_host[is.na(distal_flux_with_host)] <- 0

  # I think we should
  cecal_flux_with_host <- cecal_flux_with_host[-nearZeroVar(t(cecal_flux_with_host)), ]
  proximal_flux_with_host <- proximal_flux_with_host[-nearZeroVar(t(proximal_flux_with_host)), ]
  distal_flux_with_host <- distal_flux_with_host[-nearZeroVar(t(distal_flux_with_host)), ]


  correlated_fluxes_host <- calculate_p_values(cecal_flux_with_host, proximal_flux_with_host, distal_flux_with_host,
                                               age = age, fluxes = TRUE)
  correlated_fluxes_host <- correlated_fluxes_host %>% filter(p.adj < 0.05)
  print("DW")
  return(correlated_fluxes_host)
}


cecal_community_fluxes_host <- "../submodules/MicrobiomeGS/cecum_com/fluxes/flux_with_host.csv"
proximal_community_fluxes_host <- "../submodules/MicrobiomeGS/proximal_com/fluxes/flux_with_host.csv"
distal_community_fluxes_host <- "../submodules/MicrobiomeGS/distal_com//fluxes/flux_with_host.csv" 


cecal_community_fluxes_among_bac <- "../submodules/MicrobiomeGS/cecum_com/fluxes/flux_among_bacteria.csv"
proximal_community_fluxes_among_bac <- "../submodules/MicrobiomeGS/proximal_com/fluxes/flux_among_bacteria.csv"
distal_community_fluxes_among_bac <- "../submodules/MicrobiomeGS/distal_com//fluxes/flux_among_bacteria.csv"


cecal_community_fluxes_internal <- "../submodules/MicrobiomeGS/cecum_com/fluxes/internal_reacts.csv"
proximal_community_fluxes_internal <- "../submodules/MicrobiomeGS/proximal_com/fluxes/internal_reacts.csv"
distal_community_fluxes_internal <- "../submodules/MicrobiomeGS/distal_com/fluxes/internal_reacts.csv"

host_fluxes_2 <- calculated_fluxes(cecal_community_fluxes_host,
                                   proximal_community_fluxes_host,
                                   distal_community_fluxes_host, age = 2)


host_fluxes_9 <- calculated_fluxes(cecal_community_fluxes_host,
                                   proximal_community_fluxes_host,
                                   distal_community_fluxes_host, age = 9)

host_fluxes_15 <- calculated_fluxes(cecal_community_fluxes_host,
                                    proximal_community_fluxes_host,
                                    distal_community_fluxes_host, age = 15)

host_fluxes_24 <- calculated_fluxes(cecal_community_fluxes_host,
                                    proximal_community_fluxes_host,
                                    distal_community_fluxes_host, age = 24)

host_fluxes_30 <- calculated_fluxes(cecal_community_fluxes_host,
                                    proximal_community_fluxes_host,
                                    distal_community_fluxes_host, age = 30)



among_bac_fluxes_2 <- calculated_fluxes(cecal_community_fluxes_among_bac,
                                        proximal_community_fluxes_among_bac,
                                        distal_community_fluxes_among_bac, age = 2)


among_bac_fluxes_9 <- calculated_fluxes(cecal_community_fluxes_among_bac,
                                        proximal_community_fluxes_among_bac,
                                        distal_community_fluxes_among_bac, age = 9)

among_bac_fluxes_15 <- calculated_fluxes(cecal_community_fluxes_among_bac,
                                         proximal_community_fluxes_among_bac,
                                         distal_community_fluxes_among_bac, age = 15)

among_bac_fluxes_24 <- calculated_fluxes(cecal_community_fluxes_among_bac,
                                         proximal_community_fluxes_among_bac,
                                         distal_community_fluxes_among_bac, age = 24)

among_bac_fluxes_30 <- calculated_fluxes(cecal_community_fluxes_among_bac,
                                         proximal_community_fluxes_among_bac,
                                         distal_community_fluxes_among_bac, age = 30)



internal_fluxes_2 <- calculated_fluxes(cecal_community_fluxes_internal,
                                       proximal_community_fluxes_internal,
                                       distal_community_fluxes_internal, age = 2)


internal_fluxes_9 <- calculated_fluxes(cecal_community_fluxes_internal,
                                       proximal_community_fluxes_internal,
                                       distal_community_fluxes_internal, age = 9)

internal_fluxes_15 <- calculated_fluxes(cecal_community_fluxes_internal,
                                        proximal_community_fluxes_internal,
                                        distal_community_fluxes_internal, age = 15)

internal_fluxes_24 <- calculated_fluxes(cecal_community_fluxes_internal,
                                        proximal_community_fluxes_internal,
                                        distal_community_fluxes_internal, age = 24)

internal_fluxes_30 <- calculated_fluxes(cecal_community_fluxes_internal,
                                        proximal_community_fluxes_internal,
                                        distal_community_fluxes_internal, age = 30)
