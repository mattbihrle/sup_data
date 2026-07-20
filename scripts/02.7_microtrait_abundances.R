# Create a microtable object of abundances of each microtrait object

# Load libraries
library(tidyverse)
library(microeco)

# Load required other data
load("output/data/mt_mag_large.RData")
load("output/data/genomeset_results.rds")
bin_abund <- mt_mag_large$otu_table
sample_meta <- mt_mag_large$sample_table

# Create table where each g3 trait is a row and bins are columns. Values represent presence/absence.
g3_traits <- genomeset_results |>
  pluck("trait_matrixatgranularity3") |>
  rename_with(~str_to_lower(.x)) |>
  column_to_rownames("id") |>
  t() |>
  as.data.frame() |> 
  mutate(across(everything(), ~ifelse(.x == 0, 0, 1)))

# -------------------------------------------------------------------
# STEP 3: Calculate Trait Abundances across samples
# -------------------------------------------------------------------


# -------------------------------------------------------------------
# Matrix multiplication: 
# (Genes x Bins matrix) %*% (Bins x Samples matrix) = (Genes x Samples matrix)
gene_abund_matrix <- as.matrix(g3_traits) %*% as.matrix(bin_abund)

# Convert back to a dataframe for the microtable
gene_abund <- as.data.frame(gene_abund_matrix)


# Create a tax_table of the heirarchy of the g3 traits from microtrait using
# the hierarchy pdf from the microtrait github

# First create a function for the split list
# Function to split and repeat the last object
split_and_fill <- function(x, n_cols) {
  parts <- str_split(x, ":")[[1]]
  
  # If we have too few parts, pad with the last part
  if (length(parts) < n_cols) {
    last_part <- tail(parts, 1)
    parts <- c(parts, rep(last_part, n_cols - length(parts))) |>
      as_vector()
  }
  
  return(parts)
}

  tax_table <- genomeset_results |>
    pluck("trait_matrixatgranularity3") |>
    select(-id) |>
    rename_with(~str_to_lower(.x)) |>
    t() |>
    as.data.frame() |>
    rownames_to_column("g3_traits") |>
    select("g3_traits") |>
    mutate(granularity = "g3") |>
    slice(1:189) |> 
    mutate(split_list = map_chr(g3_traits, ~paste(split_and_fill(.x, n_cols = 7), 
                              collapse = ":"))) |>
    separate_wider_delim(cols = split_list, 
                         delim = ":", 
                         names = c("l1", "l2", "l3", "l4", "l5", "l6", "l7"), 
                         too_few = "error", 
                         too_many = "error") |>
    column_to_rownames("g3_traits") |>
    as.data.frame()



# Create new mt object
mt_microtrait <- microtable$new(otu_table = gene_abund, 
                              sample_table = sample_meta, 
                              tax_table = tax_table)

mt_microtrait$tidy_dataset()

mt_microtrait$cal_abund()
mt_microtrait$cal_alphadiv()
mt_microtrait$cal_betadiv(method = c("bray","aitchison","robust.aitchison", "jaccard"))
# Save microtable object

save(mt_microtrait, file = "output/data/mt_microtrait.RData")
obs <- ls() |> 
  str_subset("mt_.*", negate = T)

rm(list = obs)
rm(obs)



