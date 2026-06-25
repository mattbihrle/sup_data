# Load libraries
library(tidyverse)
library(microeco)

# Load required other data
load("output/data/mt_mag.RData")
# Assuming your existing microtable is named 'mt_mag'
bin_abund <- mt_mag$otu_table
sample_meta <- mt_mag$sample_table

# Create foam table where genes are rows and bins are columns. Values represent precense/absence.
foam_rollup_df <- read_tsv("output/data/metabolism/foam_rollup_all.tsv")

length(unique(foam_rollup_df$bin))
# files <- list.files("msi_downloads/metacerb/all_files/", full.names = T)
# files[1]

# foam_files <- files |>
#   str_extract(".*FOAM.*") |> 
#   na.omit()

# # Create FOAM df -----------------------------------------------------------
# for(i in 1:length(foam_files)){
# bin_name <- str_extract(foam_files[i], "WM.*_.*_.*_.*_[0-9]{6}") |> 
#   str_remove('MAGScoT_cleanbin_000')
  
# print(bin_name)
  
#   int_df <- read_tsv(foam_files[i]) |> 
#       mutate(bin = bin_name, .before = Id) |> 
#       rename_with(tolower)
#   # bin_name <- core_bins[i]
#   # #Set file path here
#   # path <- file.path("msi_downloads", "metacerb", paste0("MetaCerb_",bin_name, ".fa"),
#   #         "step_10-visualizeData",paste0("prodigal_",bin_name), "/")
#   # path
#   # Select the file you want with the regex here
  
#   if(i == 1){
#     gene_df <- int_df
#   } else {
#     gene_df <- bind_rows(gene_df, int_df)
#   }
# }
# foam_df <- gene_df
# Cut out to retain the bins that are in the mag tax table
foam_df_cut <- foam_rollup_df |> 
    slice(which(foam_rollup_df$bin %in% rownames(mt_mag$tax_table)))

  
  foam_table <- foam_df_cut |> 
    group_by(bin, id) |> 
      tally(count)
    foam_table |> 
      ggplot(aes(x = n)) +
        geom_histogram()
      
# check to make sure that I have the correct bins
      
length(unique(foam_df_cut$bin))
      
# Okay it looks like there are not really many issues where need to worry about there being more than 1 gene copy.

foam_table_wide <- foam_table |> 
    pivot_wider(names_from = bin, values_from = n) |> 
    column_to_rownames("id") |> 
    # Make it a presence abesnce thing (change things here if I need to)
    mutate(across(everything(), ~ifelse(is.na(.x), 0, 1)))

# -------------------------------------------------------------------
# STEP 3: Calculate Gene Abundances across samples
# -------------------------------------------------------------------
# Matrix multiplication: 
# (Genes x Bins matrix) %*% (Bins x Samples matrix) = (Genes x Samples matrix)
gene_abund_matrix <- as.matrix(foam_table_wide) %*% as.matrix(bin_abund)

# Convert back to a dataframe for the microtable
gene_abund <- as.data.frame(gene_abund_matrix)


# Create a tax_table, but only keep one of the l1, l2, l3 explanations. Here we are losing some information about what gene pathways genes are being used for. 

tax_table <- foam_df_cut |> 
  mutate(func = `function`) |> 
    select(bin, id, l1, l2, l3, l4, func) |> 
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
obs <- ls()

rm(list = obs)
rm(obs)
