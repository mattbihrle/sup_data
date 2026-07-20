# Import and create microtable of CAZy genes

# load libraries
library(microeco)
library(tidyverse)

# Load in mt_mag_large
load("output/data/mt_mag_large.RData")

# Import all the cazy enzyme data

files <- list.files("msi_downloads/clean_bins_cazy_rollup_all/", full.names = TRUE)

bins <- mt_mag_large$tax_table$bin |>
  str_remove("b__") |>
  str_split_fixed("_", n = 3) |>
as.data.frame()

bins <- paste0(bins$V1, "_", bins$V2,"_", "MAGScoT_cleanbin_000", bins$V3)

file_start <- "msi_downloads/clean_bins_cazy_rollup_all/rollup_prodigal_"
file_end <- "-CAZy.tsv"


files <- paste0(file_start, bins, file_end)

for(i in 1:length(files)){
  # i <- 1
  ii <- read_tsv(files[i]) |>
    mutate(bin = str_extract(files[i],"WM.{2}_.{2,3}.*[0-9]{6}")) |>
    mutate(bin = str_remove(bin, "MAG.*cleanbin_000"))
  
  if(i == 1){
    cazy_table <- ii
  } else {
    cazy_table <- bind_rows(cazy_table, ii)
  }
}


cazy_table <- cazy_table |>
  rename_with(~str_to_lower(.x)) |>
  rename_with(~paste0("cazy_", .x), .cols = -bin) |>
  write_tsv("output/data/cazy_table.tsv")

# Now start making it's own microtable -----------------------------------------
bin_abund <- mt_mag_large$otu_table
sample_meta <- mt_mag_large$sample_table

# Create table where each cazy is a row and bins are columns. Values represent presence/absence.
cazy_presence <- cazy_table |> 
  select(bin, cazy_id, cazy_count) |>
  mutate(cazy_count = 1) |>
  pivot_wider(names_from = bin, values_from = cazy_count, values_fill = 0) |>
  column_to_rownames("cazy_id")
  
rownames(cazy_presence) == rownames(bin_abund)

# -------------------------------------------------------------------
# STEP 3: Calculate Trait Abundances across samples
# -------------------------------------------------------------------


# -------------------------------------------------------------------
# Matrix multiplication: 
# (Genes x Bins matrix) %*% (Bins x Samples matrix) = (Genes x Samples matrix)
gene_abund_matrix <- as.matrix(cazy_presence) %*% as.matrix(bin_abund)

# Convert back to a dataframe for the microtable
gene_abund <- as.data.frame(gene_abund_matrix)


# Create a tax table 

tax_table <- cazy_table |>
  select(-bin, -cazy_count) |>
  unique() |>
  column_to_rownames("cazy_id")



# Create new mt object
mt_cazy <- microtable$new(otu_table = gene_abund, 
                                sample_table = sample_meta, 
                                tax_table = tax_table)


mt_cazy$tidy_dataset()

mt_cazy$cal_abund()
mt_cazy$cal_alphadiv()
mt_cazy$cal_betadiv(method = c("bray","aitchison","robust.aitchison", "jaccard"))
# Save microtable object

save(mt_cazy, file = "output/data/mt_cazy.RData")
obs <- ls() |> 
  str_subset("mt_.*", negate = T)

rm(list = obs)
rm(obs)
