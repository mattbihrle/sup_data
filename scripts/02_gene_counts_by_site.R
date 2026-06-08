# 12_gene_counts_by_site

# Goals: 
# 1) Create a master table of all the foam rollup genes by site
# 2) import that table into a microtable object
# 3) Run analyses including ancombc2, aldex2, and clique abundance. 


# NOTE - This scripts creates both a kegg and foam list but only uses the foam list of genes. swap out kegg for foam in the bottom part of the script
# create a kegg microtable object

# Setup

## Load Libraries

library(tidyverse)
library(microeco)


# Load in the data

# Create FOAM df -----------------------------------------------------------


files <- list.files("msi_downloads/final_all/foam_rollup", full.names = T)
files[1]


# Test read 1 file

for(i in 1:length(files)){
sample_name <- str_extract(files[i], "WM.{2}_[:alnum:]{2,3}")
  
print(sample_name)
  
  int_df <- read_tsv(files[i]) |> 
      mutate(sample = sample_name, .before = ID) |> 
      rename_with(tolower) |> 
      mutate(func = `function`)
  
  if(i == 1){
    gene_df <- int_df
  } else {
    gene_df <- bind_rows(gene_df, int_df)
  }
}
foam_df <- gene_df

#Create KEGG df ------------------------------------------------------------------

files <- list.files("msi_downloads/final_all/kegg_rollup", full.names = T)
files[1]

for(i in 1:length(files)){

sample_name <- str_extract(files[i], "WM.{2}_[:alnum:]{2,3}")
  
print(sample_name)
  
  int_df <- read_tsv(files[i]) |> 
      mutate(sample = sample_name, .before = ID) |> 
      rename_with(tolower) |> 
      mutate(func = `function`)
  
  if(i == 1){
    gene_df <- int_df
  } else {
    gene_df <- bind_rows(gene_df, int_df)
  }
}
kegg_df <- gene_df

# Microtable------------------------------------------------------------------------------------------------------------------------------


# Now create mt object with genes, each gene is an otu and the levels are the taxa, otu table is counts

# Load required other data
load("output/data/mt_mag.RData")
# Take meta_data from mt_mag
sample_meta <- mt_mag$sample_table

# Create an otu table

otu_table <- foam_df |> 
  select(sample, count, id) |> 
  distinct() |>  
  mutate(sample = str_remove(sample, "_S.{1,2}")) |> 
  pivot_wider(names_from = "sample", values_from = "count") |> 
  column_to_rownames("id")

head(otu_table)
# Create tax table


tax_table <- foam_df |>  
    select(sample, id, l1, l2, l3, l4, func) |> 
    distinct() |> 
    mutate(sample = str_remove(sample, "_S.{1,2}")) |>
    group_by(id) |> 
    slice_head() |> 
    mutate(gene = id) |> 
    column_to_rownames("id")

head(tax_table)


# Create new mt object
mt_gene_sample <- microtable$new(otu_table = otu_table, 
                             sample_table = sample_meta, 
                             tax_table = tax_table)

# Save microtable object

save(mt_gene_sample, file = "output/data/mt_gene_sample.RData")

# # Remove unneeded data----------------------------------------------------
remove <- ls()

rm(list = remove)

rm(remove)


load("output/data/mt_gene_sample.RData")


# Create Trans func and run ancomb
test <- "test"

print(test)
# Run this line only if you need to recreate the gene object, it will take a long time. 
t1_sample_ancomb <- trans_diff$new(dataset = mt_gene_sample, method = "ancombc2", group = "strat_season", fix_formula = "strat_season", alpha = 0.05, taxa_level = 'gene')

print("finished running ancomb, saving")
save(t1_sample_ancomb, file = "output/data/t1_sample_ancomb.RData")
print("finished with everything!!")