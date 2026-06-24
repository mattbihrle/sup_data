#13_microtrait 

# Microtrait and grodon run previously on msi servers and microtrait_results.rds downloaded with filezilla


load("msi_downloads/microtrait_results.rds")
library(microtrait)
library(tidyverse)
test <- microtrait_results
for(i in 1:length(microtrait_results)){
  test[[i]]$rds_file <- microtrait_results[[i]]$rds_file |> 
    str_replace_all("../clean_bins_6_19", "msi_downloads/rds_files") |> 
    str_replace_all( "/", "\\\\") 
}
# Extract individual file names
rds_files = unlist(parallel::mclapply(test, "[[", "rds_file"))
# create a single output 
genomeset_results = make.genomeset.results(rds_files = rds_files,
                                           ids = sub(".microtrait.rds", "", basename(rds_files)),
                                           ncores = 1)

save(genomeset_results, file = "genomeset_results.rds")

# Start here
# Test

base_dir = system.file("extdata/precomputed",package = "microtrait")
dataset = "environmentalgenomes"
genomeset_results = readRDS(file.path(base_dir, paste0(dataset, ".microtraitresults.rds")))
lapply(genomeset_results, dim)

library(dplyr)
genome_metadata = readRDS(file.path(base_dir, paste0(dataset, ".metadata.rds")))
genomeset_results_wmetadata = microtrait::add.metadata(genomeset_results, genome_metadata, genome_metadata_idcol = "IMG Taxon ID") |> convert_traitdatatype(binarytype = "logical")
saveRDS(genomeset_results_wmetadata, file.path(base_dir, paste0(dataset, ".microtraitresults.rds")))

microtrait_results[[1]]
microtrait_results[[1]]$rds_file
