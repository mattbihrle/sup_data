library(microeco)
library(tidyverse)
# Master microtable of all genes
# Initial tables ---------------------------------------------------------------
load("output/data/mt_mag_large.RData")
bin_abund <- mt_mag_large$otu_table
sample_meta <- mt_mag_large$sample_table
# CAZy--------------------------------------------------------------------------
load("output/data/mt_cazy.RData")
mt_cazy$otu_table

# Create a tax_table

cazy_tax <- mt_cazy$tax_table |> 
  mutate(origin = "cazy") |> 
  select(origin, matches("cazy_")) |> 
  rownames_to_column("id")
# Microtrait--------------------------------------------------------------------
load("output/data/mt_microtrait.RData")

# Same thing but for the hmm and hmm

g3_tax <- mt_microtrait$tax_table |> 
  mutate(origin = "microtrait") |> 
  rownames_to_column("id") |> 
  rename_with(.cols = -c(id, origin), ~paste0("g3_", .x))


hmm_table <- mt_mag_large$tax_table |> 
  select(matches("hmm_")) |> 
  mutate(across(everything(), ~ifelse(.x == TRUE, 1, 0))) |> 
  t() |> 
  as.data.frame()

hmm_tax <- hmm_table |> 
  mutate(origin = "microtrait", 
         id = rownames(hmm_table)) |> 
  select(origin, id)
# Multiply

# (Genes x Bins matrix) %*% (Bins x Samples matrix) = (Genes x Samples matrix)
hmm_abund_matrix <- as.matrix(hmm_table) %*% as.matrix(bin_abund)

# Convert back to a dataframe for the microtable
hmm_abund <- as.data.frame(hmm_abund_matrix)

#create a tax_table



# KEGG -------------------------------------------------------------------------
load("output/data/mt_gene_mag.RData")

kegg_tax <- mt_gene_mag$tax_table |> 
  rename_with(.cols = -bin, ~paste0("kegg_", .x)) |> 
  mutate(origin = "kegg") |> 
  rownames_to_column("id")

# Combine it all----------------------------------------------------------------
# Now create a big otu_table

big_otu <- bind_rows(mt_gene_mag$otu_table, mt_microtrait$otu_table, 
                     hmm_abund, mt_cazy$otu_table)

big_tax <- cazy_tax |> 
  full_join(kegg_tax, by = c("id", "origin")) |> 
  full_join(g3_tax, by = c("id","origin")) |> 
  full_join(hmm_tax, by = c("id","origin")) |> 
  column_to_rownames("id")


mt_all_genes <- microtable$new(tax_table = big_tax, sample_table = sample_meta, 
                               otu_table = big_otu)
mt_all_genes
mt_all_genes$tidy_dataset()  
mt_all_genes  

mt_all_genes$cal_abund()
mt_all_genes$cal_alphadiv()
mt_all_genes$cal_betadiv(method = c("bray","aitchison","robust.aitchison", "jaccard"))
# Save microtable object

save(mt_all_genes, file = "output/data/mt_all_genes.RData")
obs <- ls() |> 
  str_subset("mt_all.*", negate = T)

rm(list = obs)
rm(obs)