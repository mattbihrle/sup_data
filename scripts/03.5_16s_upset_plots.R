# Upset plots
library(microeco)
load("output/data/mt_16s.RData")

mt_16s$filter_taxa(
  rel_abund = 0.001
)

mt_16s
venn <- trans_venn$new(mt_16s)
