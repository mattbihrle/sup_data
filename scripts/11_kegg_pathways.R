# Install libraries that I may need

# BiocManager::install("pathview")

# BiocManager::install("KEGGREST")
BiocManager::install("gage")
#Load Libraries
library(tidyverse)
library(pathview)
library(KEGGREST)
library(gage)


# Read in the foam rollup data

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
  read_tsv(show_col_types = F) |> 
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

# Load in the mag_data

load("output/data/mt_mag.RData")

core_bins <- mt_mag$tax_table$bin |> 
  str_remove("b__")

head(core_bins)
head(foam_rollup_df$bin)
# Create a new column to filter out bins that I have abundance data for
foam_rollup_df <- foam_rollup_df |> 
  mutate(core_bin = foam_rollup_df$bin %in% core_bins)

# Now start following tutorial from pathview

foam_int <- foam_rollup_df |> 
  filter(core_bin == "TRUE") |> 
  group_by(bin, id) |> 
  tally(name = "count")

keggFind(database = "pathway", query = "Nitrogen Cycle")

n_cycle_id <- 01310

foam_path <- foam_int |> 
  select(bin, id, count) |> 
  pivot_wider(names_from = "bin", 
  values_from = "count", 
values_fill = NA) |> 
  column_to_rownames("id") |> 
  as.matrix()

pv_bin_0 <- pathview(
  gene.data = foam_path[, "WM01_S1_001"],
  pathway.id = n_cycle_id,
  species = "ko",
  out.suffix = "pv_bin_0"
)

mt_chitin <- clone(mt_mag)

mt_chitin$tax_table <- mt_chitin$tax_table |> 
  slice(str_which(mt_chitin$tax_table$f, ".*Chitin"))

mt_chitin$tidy_dataset()$cal_abund()
plot <- trans_abund$new(mt_chitin, taxrank = "bin", ntaxa = "30", high_level = "g")$
  plot_bar(ggnested = TRUE, xtext_angle = 30, facet = "strat_season")
plot
plotly::plotly_build(plot)
