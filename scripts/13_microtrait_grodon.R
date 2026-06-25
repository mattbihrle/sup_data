#13_microtrait 

# Microtrait and grodon run previously on msi servers and microtrait_results.rds downloaded with filezilla


load("msi_downloads/microtrait_results.rds")
library(microtrait)
library(tidyverse)
# test <- microtrait_results
for(i in 1:length(microtrait_results)){
  microtrait_results[[i]]$rds_file <- microtrait_results[[i]]$rds_file |> 
    str_replace_all("../clean_bins_6_19", "msi_downloads/rds_files") |> 
    str_replace_all( "/", "\\\\") 
  microtrait_results[[i]]$microtrait_result$id <- 
    str_remove(microtrait_results[[i]]$microtrait_result$id, 
      "MAG.*_000")
}
microtrait_results[[1]]$microtrait_result$id
# Keep only the bins that made it through contamination check in mt_mag
load("output/data/mt_mag.RData")

bins <- rownames(mt_mag$tax_table)
head(bins)

test <- rep(NA, length(microtrait_results))
for(i in 1:length(microtrait_results)){
  if(microtrait_results[[i]]$microtrait_result$id %in% bins) {
    test[i] <- TRUE
  } else {
    test[i] <- FALSE
  }
}
table(test)
# Extract individual file names
rds_files = unlist(parallel::mclapply(microtrait_results, "[[", "rds_file"))

rds_files <- rds_files[which(test)]


# create a single output 
# # Skip here if already run ----------------------------------------------------
# genomeset_results = make.genomeset.results(rds_files = rds_files,
#                                            ids = bins,
#                                            ncores = 1)
# # --------------------------------------------------------------------------

# # Next add metadata

# metadata <-bind_cols(mt_mag$tax_table, mt_mag$otu_table) |> 
#   mutate(id = rownames(metadata))

# genomeset_results <- microtrait::add.metadata(genomeset_results, metadata, genome_metadata_idcol = "id")
# save(genomeset_results, file = "output/data/genomeset_results.rds")
load("output/data/genomeset_results.rds")

# Normalize count traits

genomeset_results_norm <- genomeset_results |> 
  trait.normalize(normby = "Genome_Size")

traits = traits_listbygranularity[[3]] %>%
  dplyr::select(`microtrait_trait-name`) %>%
  dplyr::filter(`microtrait_trait-name` != "Resource Use:Chemotrophy:chemolithoautotrophy:anaerobic ammonia oxidation") %>%
  dplyr::pull(`microtrait_trait-name`)
trait_matrixatgranularity3 = genomeset_results_wmetadata_norm[["trait_matrixatgranularity3"]] %>%
  dplyr::filter(`Ecosystem_Category` %in% c("Terrestrial", "Plants")) %>%
  dplyr::filter(`Ecosystem_Type` %in% c("Soil", "Rhizoplane", "Rhizosphere", "Roots")) %>%
  dplyr::select(c("id", traits, "mingentime", "ogt")) %>%
  #dplyr::slice(1:100) %>%
  dplyr::filter(`mingentime` < 100) %>%   # max out mingentime at 25 days to avoid errors
  dplyr::filter(!is.na(`mingentime`) & !is.na(`ogt`))
# Test
base_dir = system.file("extdata/precomputed",package = "microtrait")
dataset = "environmentalgenomes"
test_genomeset_results = readRDS(file.path(base_dir, paste0(dataset, ".microtraitresults.rds")))
lapply(test_genomeset_results, dim)

library(dplyr)
genome_metadata = readRDS(file.path(base_dir, paste0(dataset, ".metadata.rds")))
genomeset_results_wmetadata = microtrait::add.metadata(test_genomeset_results, genome_metadata, genome_metadata_idcol = "IMG Taxon ID") |> convert_traitdatatype(binarytype = "logical")
# saveRDS(genomeset_results_wmetadata, file.path(base_dir, paste0(dataset, ".microtraitresults.rds")))

microtrait_results[[1]]
microtrait_results[[1]]$rds_file
