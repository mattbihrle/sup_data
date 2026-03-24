# Load Libraries
library(tidyverse)

# Run previous scripts/import data needed

de_qual_rpt <- read_tsv("data/derep_bins/all_bins_quality_report.tsv")
de_coverage <- read_tsv("data/derep_bins/Drep_Bins_coverage.tsv")
  
de_coverage$Genome <- str_replace_all(de_coverage$Genome, pattern = "_MAGScoT_cleanbin_000", replacement = "_")

colnames(de_coverage) <- str_remove_all(colnames(de_coverage), "001_val_.{2}fq.gz ")

