# Load required libraries
library(tidyverse)
library(phyloseq)
library(microeco)


# Load in data
  # First the sample names
  samp_names <- readr::read_tsv("lotus3_out/final_sample_map_4_R.txt")
  samp_names <- setNames(samp_names$SampleID, samp_names$RealNames)
  # Then load in taxonomic data
  tax_df <- readr::read_tsv("lotus3_out/hiera_BLAST.txt") 
  tax_df <- tax_df |> 
    dplyr::mutate_at(colnames(tax_df), ~ str_replace_all(., "\\?", "unknown"))
  # Then the otu data, rename, and join with tax data
  otu_df <- readr::read_tsv("lotus3_out/OTU.txt") |> 
  dplyr::rename("OTU" = "...1", any_of(samp_names)) |> 
    dplyr::full_join(tax_df)
long_otu_df <- pivot_longer(otu_df, cols = str_extract_all(colnames(otu_df), "WM.*"))
  # Load phyloseq object
  p_df <- phyloseq::import_biom("lotus3_out/OTU.biom")
  # Check it out
  ntaxa(p_df)
  rank_names(p_df)
  sample_variables(p_df)
# Initial Graphs
  # plot OTU abunance
p_df |> 
plot_bar(p_df, y = "Abundance", fill = "Rank2")

