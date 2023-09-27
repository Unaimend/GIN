library(dplyr)
setwd("/home/td/Projects/Aging/Baseline/GIN/code")
two_up = readRDS("../data/enricher_up_2M.RDS")@result
two_up = two_up %>% filter(p.adjust < 0.05)

nine_up = readRDS("../data/enricher_up_9M.RDS")@result
nine_up = nine_up %>% filter(p.adjust < 0.05)

fifteen_up = readRDS("../data/enricher_up_15M.RDS")@result
fifteen_up = fifteen_up %>% filter(p.adjust < 0.05)

shared_sse_up = Reduce(intersect, list(two_up$Description, nine_up$Description, fifteeen_up$Description))

intersect(two_up$Description, nine_up$Description)
intersect(nine_up$Description, fifteeen_up$Description)
