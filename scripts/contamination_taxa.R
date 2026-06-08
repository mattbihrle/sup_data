# Load libraries

library(microeco)
library(tidyverse)

# Load data 
load("mt_16s_contam.RData")
mag_contam <- read_tsv("mag_contamination.tsv")

#Possible contamination taxa in 16S
view(mt_16s_contam$tax_table)
# Plot of major contaminant OTUs in 16S data grouped by genus

trans_abund$new(mt_16s_contam, taxrank = "otu", ntaxa = 30)$plot_bar()
trans_abund$new(mt_16s_contam, taxrank = "otu", ntaxa = 50, high_level = "g")$plot_bar(ggnested = T)



