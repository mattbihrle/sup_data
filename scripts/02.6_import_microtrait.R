#13_microtrait 

# Microtrait and grodon run previously on msi servers and microtrait_results.rds downloaded with filezilla


# load("msi_downloads/microtrait_results.rds")
library(microtrait)
library(tidyverse)
# # test <- microtrait_results
# for(i in 1:length(microtrait_results)){
#   microtrait_results[[i]]$rds_file <- microtrait_results[[i]]$rds_file |> 
#     str_replace_all("../clean_bins_6_19", "msi_downloads/rds_files") |> 
#     str_replace_all( "/", "\\\\") 
#   microtrait_results[[i]]$microtrait_result$id <- 
#     str_remove(microtrait_results[[i]]$microtrait_result$id, 
#       "MAG.*_000")
# }
# # microtrait_results[[1]]$microtrait_result$id
# Keep only the bins that made it through contamination check in mt_mag
load("output/data/mt_mag.RData")

# bins <- rownames(mt_mag$tax_table)
# head(bins)

# test <- rep(NA, length(microtrait_results))
# for(i in 1:length(microtrait_results)){
#   if(microtrait_results[[i]]$microtrait_result$id %in% bins) {
#     test[i] <- TRUE
#   } else {
#     test[i] <- FALSE
#   }
# }
# table(test)
# # Extract individual file names
# rds_files = unlist(parallel::mclapply(microtrait_results, "[[", "rds_file"))

# rds_files <- rds_files[which(test)]


# # create a single output 
# # # Skip here if already run ----------------------------------------------------
# # genomeset_results = make.genomeset.results(rds_files = rds_files,
# #                                            ids = bins,
# #                                            ncores = 1)
# # # --------------------------------------------------------------------------

# # # Next add metadata

# metadata <-bind_cols(mt_mag$tax_table, mt_mag$otu_table) 
# metadata <- metadata |> 
#   mutate(id = rownames(metadata))
# head(metadata)


# genomeset_results <- microtrait::add.metadata(genomeset_results, metadata, genome_metadata_idcol = "id")
# save(genomeset_results, file = "output/data/genomeset_results.rds")
load("output/data/genomeset_results.rds")

#Normalize genome results based on the genome length
genomeset_results_norm <- genomeset_results |> 
  trait.normalize(normby = "genome_length")

# Remove genomeset_results that are all zero
test <- genomeset_results_norm |> 
  # map(~column_to_rownames("id")) 
  map( 
    ~select(.x, where(~ !all(.x == 0)))) |> 
  # remove the anything that is not part of the original table 
  map(~select(.x, !d:WM16)
  )
# Remove the columns that are metadata

colnames(test[[1]])[which(grepl(".*:.*", colnames(test[[1]])))] <- 
  paste0("g1_", colnames(test[[1]])[which(grepl(".*:.*", colnames(test[[1]])))])

colnames(test[[2]])[which(grepl(".*:.*", colnames(test[[2]])))] <- 
  paste0("g2_", colnames(test[[2]])[which(grepl(".*:.*", colnames(test[[2]])))])

colnames(test[[3]])[which(grepl(".*:.*", colnames(test[[3]])))] <- 
  paste0("g3_", colnames(test[[3]])[which(grepl(".*:.*", colnames(test[[3]])))])

test[[4]] <- test[[4]] |> 
  rename_with(~paste0("hmm_", .x), .cols = !matches("^(id|genome_length|optimumT|mingentime)")
  )

colnames(test[[4]])



test[[5]] <- test[[5]] |> 
  rename_with(~paste0("rule_", .x), .cols = !matches("^(id|genome_length|optimumT|mingentime)")
  )
colnames(test[[5]])


# view(test$trait_matrixatgranularity1)
# view(mt_mag$tax_table)


mt_mag$tax_table <- mt_mag$tax_table |> 
  rownames_to_column("id") |> 
  left_join(test$trait_matrixatgranularity1, by = "id") |> 
  left_join(test$trait_matrixatgranularity2, by = c("id", "genome_length", "optimumT", "mingentime")) |> 
  left_join(test$trait_matrixatgranularity3, by = c("id", "genome_length", "optimumT", "mingentime")) |> 
  left_join(test$hmm_matrix, by = c("id", "genome_length", "optimumT", "mingentime")) |> 
  left_join(test$rule_matrix, by = c("id", "genome_length", "optimumT", "mingentime")) |> 
  column_to_rownames("id")

# view(test$hmm_matrix)
# view(genomeset_results_norm$hmm_matrix)
# view(mt_mag$tax_table)

mt_mag_large <- mt_mag

mt_mag_large$tax_table <- mt_mag_large$tax_table |> 
  rename_with(~str_replace_all(., " ", "_")) |> 
  rename_with(~str_replace_all(., ":", "__")) |>
  rename_with(~str_to_lower(.)) |> 
  rename_with(~str_replace_all(., "[^[a-z][0-9]_.]", "_"))


mt_mag_large$tax_table |> 
  select(contains("denitrification")) |> 
  colnames()
save(mt_mag_large, file = "output/data/mt_mag_large.RData")

# Make a quick df of the column names for later lookup

microtrait_columns <- mt_mag_large$tax_table |> 
  select(matches("^(g[0-9]{1}|hmm|rule)")) |> 
  colnames() |>
  as.data.frame() |>
  rename("columns" = matches(".*")) |> 
  write_csv(file = "output/data/microtrait_columns.csv")

# Okay from here let me start to look for autotrophs and such
library(microeco)
mt_mag_auto <- clone(mt_mag_large)


mt_mag_auto$tax_table <- mt_mag_auto$tax_table |> 
  select(d, p, c, o , f, g, s, bin, matches("g1_.*resource.*")) |> 
  mutate(poss_chemotroph = if_any(matches("g1_.*chemolith.*"), ~.x != 0), 
         poss_phototroph = if_any(matches("g1_.*phototro.*"), ~.x != 0), 
         poss_heterotroph = if_any(matches("g1_.*organohetero.*"), ~.x != 0))


# view(mt_mag_auto$tax_table)

# create a little df of those new columns

metabolism_df <- mt_mag_auto$tax_table |> 
  select(bin, matches("poss.*"))

#Add that back into the main microtable 

mt_mag_large$tax_table <- mt_mag_large$tax_table |> 
  left_join(metabolism_df, by = "bin")
# Turn these all into presence/absence

mt_mag_large$tax_table <- mt_mag_large$tax_table |> 
  mutate(across(matches(c("g1", "g2", "g3", "hmm", "rule")), ~ifelse(. == 0, FALSE, TRUE)))

#Lastly make sure that the row names are included in the tax table
mt_mag_large$tax_table <- mt_mag_large$tax_table |> 
  mutate(rowname = str_remove(bin, "b__")) |> 
  column_to_rownames("rowname")

mt_mag_large$tidy_dataset()

mt_mag_large$cal_abund()

# Calc Alpha Diversity
mt_mag_large$cal_alphadiv()

# Calc Bray Curtis Dissimilarity
mt_mag_large$cal_betadiv(method = c("bray","aitchison","robust.aitchison", "jaccard"), unifrac = T)
head(mt_mag_large$beta_diversity)

save(mt_mag_large, file = "output/data/mt_mag_large.RData")

remove <- ls() |> 
  str_subset("mt_.*", negate = T)

rm(list = remove)
rm(remove)
