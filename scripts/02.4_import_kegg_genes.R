#02.4_import and clean all genes

# Load Libraries
library(microeco)
library(tidyverse)

# Load in mt_mag for the list of clean_bins
load("output/data/mt_mag.RData")

clean_bins <- mt_mag$tax_table$bin |> 
  str_remove(pattern = "b__") 

clean_bins_df <- clean_bins |> 
  str_split_fixed(pattern = "(?<=S.{1,2}_)", n = 2) |> 
  as.data.frame()

long_bins <- paste0(clean_bins_df$V1, "MAGScoT_cleanbin_000", clean_bins_df$V2)

files <- list.files("msi_downloads/rollup/kegg", full.names = T)

file_start <- str_extract(files[1], ".*(rollup_prodigal_=?)")
file_end <- str_extract(files[1], "(?<=cleanbin_[0-9]{6}).*")

files <- paste0(file_start, long_bins, file_end) |> 
  sort()


kegg_rollup <- files |> 
  map(~ read_tsv(.x) |> 
           mutate(bin = str_extract(.x, "WM.*[0-9]{6}")) |> 
           mutate(bin = str_remove(bin, pattern = "MAGScoT_cleanbin_000")) |> 
           rename_with(str_to_lower, everything()) |> 
           mutate(func = `function`,
                  kegg_path = str_extract(l3, "ko[0-9]{5}"), 
                  l3 = str_remove(l3, "\\[.*\\]")) |> 
           select(-`function`)
  ) |> 
  list_rbind()

# # TEST ------------------------------------------------------------------------
# test <- read_tsv(files[6]) |> 
#   mutate(bin = str_extract(files[6], "WM[0-9]{6}")) |> 
#   mutate(bin = str_remove(bin, pattern = "MAGScoT_cleanbin_000")) |> 
#   rename_with(str_to_lower, everything()) |> 
#   mutate(func = `function`,
#          kegg_path = str_extract(l3, "ko[0-9]{5}"), 
#          l3 = str_remove(l3, "\\[.*\\]")) |> 
#   select(-`function`)
# 
# # TEST -------------------------------------------------------------------------
kegg_rollup <- kegg_rollup |> 
  mutate(across(c(l1:l3, func), str_to_lower))

# Save this file 

kegg_rollup <- kegg_rollup |> 
  write_tsv("output/data/kegg_rollup_all.tsv")

# Okay now create a folder that has a single kegg file for each bin

test <- str_remove(mt_mag$tax_table$bin, "b__" ) |> 
  map(~ kegg_rollup |> 
        filter(bin == .x) |> 
        select(bin, id) |> 
        distinct(.keep_all = TRUE) |> 
        write_tsv(file = paste0("output/data/kegg_output/", .x, "_clean_kegg.tsv")) |> 
        print()
)

paste("Saved a master spreadsheet of all kegg_rollup genes as well as a single",
      "file for each clean bin for import into kegg reconstruct")

obs <- ls()
rm(list = obs)
rm(obs)
  