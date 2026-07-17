# Load libraries
library(tidyverse)
library(microeco)

# Load required other data
load("output/data/mt_mag.RData")
# Assuming your existing microtable is named 'mt_mag'
bin_abund <- mt_mag$otu_table
sample_meta <- mt_mag$sample_table

# Create kegg table where genes are rows and bins are columns. Values represent presence/absence.

# First make a list of all the bins we want to use
bins <- mt_mag_large$tax_table$bin |>
  str_remove("b__") |>
  str_split_fixed("_", n = 3) |>
  as.data.frame()

bins <- paste0(bins$V1, "_", bins$V2,"_", "MAGScoT_cleanbin_000", bins$V3)
test <- read_tsv(files[1])

# turn that into a list of file names to import 
file_start <- "msi_downloads/rollup/kegg/rollup_prodigal_"
file_end <- "-KOFam_all_KEGG.tsv"


files <- paste0(file_start, bins, file_end)
files[1]

# Create kegg df -----------------------------------------------------------
for(i in 1:length(files)){
bin_name <- str_extract(files[i], "WM.*_.*_.*_.*_[0-9]{6}") |>
  str_remove('MAGScoT_cleanbin_000')
print(bin_name)
  int_df <- read_tsv(files[i]) |>
      mutate(bin = bin_name, .before = everything()) |>
      rename_with(tolower)

  if(i == 1){
    gene_df <- int_df
  } else {
    gene_df <- bind_rows(gene_df, int_df)
  }
}
kegg_df <- gene_df

kegg_table <- kegg_df |> 
  group_by(bin, id) |> 
  tally(count)
kegg_table |> 
  ggplot(aes(x = n)) +
  geom_histogram()

# check to make sure that I have the correct bins

length(unique(kegg_df$bin))

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
  mutate(func = `function`) |> 
  mutate(path = str_extract(l3, "ko[0-9]{5}")) |>
  mutate(l3 = str_remove(l3, "\\s\\[.*\\]")) |> 
  select(bin, id, l1, l2, l3, path, func) |> 
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
