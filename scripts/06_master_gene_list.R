# MAKE MASTER GENE LIST

#Import libraries
library(tidyverse)

files <- list.files("msi_downloads/metacerb/all_files/", full.names = T)
files[1]

foam_files <- files |>
  str_extract(".*FOAM.*")
  na.omit()

kegg_files <- files |> 
  str_extract(".*KEGG.*") |> 
  na.omit()

length(foam_files) == length(kegg_files)

# Create FOAM df -----------------------------------------------------------
for(i in 1:length(foam_files)){
bin_name <- str_extract(foam_files[i], "WM.*_.*_.*_.*_[0-9]{6}") |> 
  str_remove('MAGScoT_cleanbin_000')
  
print(bin_name)
  
  int_df <- read_tsv(foam_files[i]) |> 
      mutate(bin = bin_name, .before = Id) |> 
      rename_with(tolower)
  # bin_name <- core_bins[i]
  # #Set file path here
  # path <- file.path("msi_downloads", "metacerb", paste0("MetaCerb_",bin_name, ".fa"),
  #         "step_10-visualizeData",paste0("prodigal_",bin_name), "/")
  # path
  # Select the file you want with the regex here
  
  if(i == 1){
    gene_df <- int_df
  } else {
    gene_df <- bind_rows(gene_df, int_df)
  }
}
foam_df <- gene_df

#Create KEGG df ------------------------------------------------------------------
for(i in 1:length(kegg_files)){
bin_name <- str_extract(kegg_files[i], "WM.*_.*_.*_.*_[0-9]{6}") |> 
  str_remove('MAGScoT_cleanbin_000')
  
print(bin_name)
  
  int_df <- read_tsv(kegg_files[i]) |> 
      mutate(bin = bin_name, .before = Id) |> 
      rename_with(tolower)
  # bin_name <- core_bins[i]
  # #Set file path here
  # path <- file.path("msi_downloads", "metacerb", paste0("MetaCerb_",bin_name, ".fa"),
  #         "step_10-visualizeData",paste0("prodigal_",bin_name), "/")
  # path
  # Select the file you want with the regex here
  
  if(i == 1){
    gene_df <- int_df
  } else {
    gene_df <- bind_rows(gene_df, int_df)
  }
}
kegg_df <- gene_df

# CAZy bins -----------------------------------------------------------

# Import core bin list
core_bins <- read_tsv("msi_downloads/Derep_list.tsv", col_names = F) |> 
  mutate(bin = str_remove(X1, "^all_bins/")) |> 
  mutate(bin = str_remove(bin, ".fa$")) |> 
  mutate(X1 = NULL)

core_bins <- as.character(core_bins$bin)

for(i in 1:length(core_bins)) {

  bin_name <- core_bins[i]

  bin_name
  #Set file path here
  path <- file.path("msi_downloads", "metacerb", "core_bins_cazy", paste0("cazy_",bin_name, ".fa"),
          "step_10-visualizeData",paste0("prodigal_",bin_name), "/")
  path
  # Select the file you want with the regex here
int_df <- list.files(path, pattern = "*id.tsv", full.names = T) |> 
  read_tsv() |> 
  mutate(bin = bin_name, .before = everything()) |> 
      rename_with(tolower)
  
    if(i == 1){
    gene_df <- int_df
  } else {
    gene_df <- bind_rows(gene_df, int_df)
  }
}
cazy_df <- gene_df

cazy_df <- cazy_df |> 
  str_remove(bin, "MAGScoT_cleanbin_000")
# Multiple levels for KEGG and FOAM --------------------------------------

# First create a function to make things a bit nicer looking
bind_files <- function(path) {
#Given a file path create a list of all files, exctract the sample name 
  # and bind all files into one R object

  files <- list.files(path, full.names = T)

for(i in 1:length(files)){
bin_name <- str_extract(files[i], "WM.*_.*_.*_.*_[0-9]{6}") |> 
  str_remove('MAGScoT_cleanbin_000')
  
print(bin_name)
  
  int_df <- read_tsv(files[i]) |> 
      mutate(bin = bin_name, .before = everything()) |> 
      rename_with(tolower)
  # bin_name <- core_bins[i]
  # #Set file path here
  # path <- file.path("msi_downloads", "metacerb", paste0("MetaCerb_",bin_name, ".fa"),
  #         "step_10-visualizeData",paste0("prodigal_",bin_name), "/")
  # path
  # Select the file you want with the regex here
  
  if(i == 1){
    gene_df <- int_df
  } else {
    gene_df <- bind_rows(gene_df, int_df)
  }
  }
  return(gene_df)
}
 

# Import all rollup files ----------------------------------------------------------
bin_names <- list.files("msi_downloads/rollup/foam/") |> 
  str_extract("WM.*[0-9]{6}")
  # str_remove("MA.*_000")

length(bin_names)

# # Import core bin list
# core_bins <- read_tsv("msi_downloads/Derep_list.tsv", col_names = F) |> 
#   mutate(bin = str_remove(X1, "^all_bins/")) |> 
#   mutate(bin = str_remove(bin, ".fa$")) |> 
#   mutate(X1 = NULL)

# core_bins <- as.character(core_bins$bin)

# core_bins <- core_bins |> str_remove("MAGScoT_cleanbin_000")

for(i in 1:length(bin_names)) {

  bin_name <- bin_names[i]

  bin_name
  #Set file path here
  path <- file.path("msi_downloads", "rollup", "foam")
  path
  # Select the file you want with the regex here
int_df <- list.files(path, pattern = paste0(".*", bin_name, ".*", ".tsv"), full.names = T) |> 
  read_tsv() |> 
  mutate(bin = bin_name, .before = everything()) |> 
      rename_with(tolower)
  
    if(i == 1){
    gene_df <- int_df
  } else {
    gene_df <- bind_rows(gene_df, int_df)
  }
}
foam_rollup_df <- gene_df

foam_rollup_df$bin <- foam_rollup_df$bin |> 
  str_remove(pattern = "MAGScoT_cleanbin_000")

# Then for KEGG -----------------------------------------------------------------

for(i in 1:length(bin_names)) {

  bin_name <- bin_names[i]

  bin_name
  #Set file path here
  path <- file.path("msi_downloads", "rollup", "kegg")
  path
  # Select the file you want with the regex here
int_df <- list.files(path, pattern = paste0(".*", bin_name, ".*", ".tsv"), full.names = T) |> 
  read_tsv() |> 
  mutate(bin = bin_name, .before = everything()) |> 
      rename_with(tolower)
  
    if(i == 1){
    gene_df <- int_df
  } else {
    gene_df <- bind_rows(gene_df, int_df)
  }
}
kegg_rollup_df <- gene_df

kegg_rollup_df$bin <- kegg_rollup_df$bin |> 
  str_remove(pattern = "MAGScoT_cleanbin_000")
# Multiple levels for KEGG and FOAM --------------------------------------

# First create a function to make things a bit nicer looking
bind_files <- function(path) {
#Given a file path create a list of all files, exctract the sample name 
  # and bind all files into one R object

  files <- list.files(path, full.names = T)

for(i in 1:length(files)){
bin_name <- str_extract(files[i], "WM.*_.*_.*_.*_[0-9]{6}") |> 
  str_remove('MAGScoT_cleanbin_000')
  
print(bin_name)
  
  int_df <- read_tsv(files[i]) |> 
      mutate(bin = bin_name, .before = everything()) |> 
      rename_with(tolower)
  # bin_name <- core_bins[i]
  # #Set file path here
  # path <- file.path("msi_downloads", "metacerb", paste0("MetaCerb_",bin_name, ".fa"),
  #         "step_10-visualizeData",paste0("prodigal_",bin_name), "/")
  # path
  # Select the file you want with the regex here
  
  if(i == 1){
    gene_df <- int_df
  } else {
    gene_df <- bind_rows(gene_df, int_df)
  }
  }
  return(gene_df)
}


# -------------------------------------------------------------------------
#Run the function on each folder and save
foam_1_df <- bind_files("msi_downloads/foam/lvl_1") |> 
  write_tsv("output/data/metabolism/foam_1.tsv")

foam_2_df <- bind_files("msi_downloads/foam/lvl_2") |> 
  write_tsv("output/data/metabolism/foam_2.tsv")

foam_3_df <- bind_files("msi_downloads/foam/lvl_3") |> 
  write_tsv("output/data/metabolism/foam_3.tsv")

kegg_1_df <- bind_files("msi_downloads/kegg/lvl_1") |> 
  write_tsv("output/data/metabolism/kegg_1.tsv")

kegg_2_df <- bind_files("msi_downloads/kegg/lvl_2") |> 
  write_tsv("output/data/metabolism/kegg_2.tsv")

kegg_3_df <- bind_files("msi_downloads/kegg/lvl_3") |> 
  write_tsv("output/data/metabolism/kegg_3.tsv")

#Save outputs to tsvs  
write_tsv(foam_df, "output/data/metabolism/foam_genes.tsv")
write_tsv(kegg_df, "output/data/metabolism/kegg_genes.tsv")
write_tsv(cazy_df, "output/data/metabolism/cazy_genes.tsv")
write_tsv(foam_rollup_df, "output/data/metabolism/foam_rollup.tsv")
write_tsv(kegg_rollup_df, "output/data/metabolism/kegg_rollup.tsv")


# Reimport Cazy data to update and rewrite

cazy_df <- read_tsv("output/data/cazy_genes.tsv") |> 
  mutate(bin = str_remove(bin, "MAGScoT_cleanbin_000")) |> 
  write_tsv("output/data/metabolism/cazy_genes.tsv")


# gene_df <- gene_df |> 
#   mutate(bin = str_remove(bin, "MAGScoT_cleanbin_000")) |> 
#   separate_wider_delim(name, ";\t", names = c("name", "notes"))

# file_test <- file |> 
#   mutate(bin = str_remove(bin, "MAGScoT_cleanbin_000")) |> 
#   separate_wider_delim(name, ";", names = c("name", "notes")) |>
#   mutate(notes =str_trim(notes, side = "both"))
#   # separate_wider_delim(name, ", ", names = c("name"), too_few = "align_start", too_many = "debug")

  