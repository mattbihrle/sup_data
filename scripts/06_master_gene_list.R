# MAKE MASTER GENE LIST

#Import libraries
library(tidyverse)

#Import core bin list
core_bins <- read_tsv("msi_downloads/Derep_list.tsv", col_names = F) |> 
  mutate(bin = str_remove(X1, "^all_bins/")) |> 
  mutate(bin = str_remove(bin, ".fa$")) |> 
  mutate(X1 = NULL)

core_bins <- as.character(core_bins$bin)

for(i in 1:length(core_bins)){
  bin_name <- core_bins[i]
  #Set file path here
  path <- file.path("msi_downloads", "metacerb", paste0("MetaCerb_",bin_name, ".fa"),
          "step_10-visualizeData",paste0("prodigal_",bin_name), "/")
  path
  # Select the file you want with the regex here
 file <- list.files(path, ".*all.*FOAM.*id.tsv", full.names = T) |> 
  read_tsv() |> 
  mutate(bin = bin_name, .before = Id) |> 
  rename_with(tolower)
  
  if(i == 1){
    gene_df <- file
  } else {
    gene_df <- bind_rows(gene_df, file)
  }
}

gene_df <- gene_df |> 
  mutate(bin = str_remove(bin, "MAGScoT_cleanbin_000")) |> 
  separate_wider_delim(name, ";\t", names = c("name", "notes"))

file_test <- file |> 
  mutate(bin = str_remove(bin, "MAGScoT_cleanbin_000")) |> 
  separate_wider_delim(name, ";", names = c("name", "notes")) |>
  mutate(notes =str_trim(notes, side = "both"))
  # separate_wider_delim(name, ", ", names = c("name"), too_few = "align_start", too_many = "debug")