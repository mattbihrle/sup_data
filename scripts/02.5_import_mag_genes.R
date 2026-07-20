# Load libraries
library(tidyverse)
library(microeco)

# Load required other data
load("output/data/mt_mag.RData")
# Assuming your existing microtable is named 'mt_mag'
bin_abund <- mt_mag$otu_table
sample_meta <- mt_mag$sample_table

# Create kegg table where genes are rows and bins are columns. Values represent presence/absence.

  # Import kegg_rollup_df

  kegg_df <- read_tsv("output/data/kegg_rollup_all.tsv")
  # Sum the number of times a gene shows up in each bin
  kegg_table <- kegg_df |> 
    group_by(bin, id) |> 
    tally(count)
  
  # Plot a histogram to make sure it mostly is just once
  kegg_table |> 
    ggplot(aes(x = n)) +
    geom_histogram()

  # check to make sure that I have the correct bins
  
  length(unique(kegg_df$bin)) == length(bins)

# Okay it looks like there are not really many issues where need to worry about there being more than 1 gene copy.

kegg_table_wide <- kegg_table |> 
  pivot_wider(names_from = bin, values_from = n) |> 
  column_to_rownames("id") |> 
  # Make it a presence abesnce thing (change things here if I need to)
  mutate(across(everything(), ~ifelse(is.na(.x), 0, 1)))

# -------------------------------------------------------------------
# STEP 3: Calculate Gene Abundances across samples
# -------------------------------------------------------------------
# Matrix multiplication: 
# (Genes x Bins matrix) %*% (Bins x Samples matrix) = (Genes x Samples matrix)
gene_abund_matrix <- as.matrix(kegg_table_wide) %*% as.matrix(bin_abund)

# Convert back to a dataframe for the microtable
gene_abund <- as.data.frame(gene_abund_matrix)


# Create a tax_table, but only keep one of the l1, l2, l3 explanations. Here we are losing some information about what gene pathways genes are being used for. 

tax_table <- kegg_df |> 
  select(bin, id, l1, l2, l3, kegg_path, func) |> 
  distinct() |> 
  group_by(id) |> 
  slice_head() |> 
  mutate(gene = id) |> 
  column_to_rownames("id")


# Create new mt object
mt_gene_mag <- microtable$new(otu_table = gene_abund, 
                              sample_table = sample_meta, 
                              tax_table = tax_table)

# Save microtable object

save(mt_gene_mag, file = "output/data/mt_gene_mag.RData")
obs <- ls() |> 
    str_subset("mt_.*", negate = T)

rm(list = obs)
rm(obs)
