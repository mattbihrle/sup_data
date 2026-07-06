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

colnames(test[[1]])[which(grepl(".*:.*", colnames(test[[1]])))] <- 
  paste0("g1_", colnames(test[[1]])[which(grepl(".*:.*", colnames(test[[1]])))])

colnames(test[[2]])[which(grepl(".*:.*", colnames(test[[2]])))] <- 
  paste0("g2_", colnames(test[[2]])[which(grepl(".*:.*", colnames(test[[2]])))])

colnames(test[[3]])[which(grepl(".*:.*", colnames(test[[3]])))] <- 
  paste0("g3_", colnames(test[[3]])[which(grepl(".*:.*", colnames(test[[3]])))])

colnames(test[[4]])[which(grepl(".*:.*", colnames(test[[4]])))] <- 
  paste0("hmm_", colnames(test[[4]])[which(grepl(".*:.*", colnames(test[[4]])))])

colnames(test[[5]])[which(grepl(".*:.*", colnames(test[[5]])))] <- 
  paste0("rule_", colnames(test[[5]])[which(grepl(".*:.*", colnames(test[[5]])))])

view(test$trait_matrixatgranularity1)
view(mt_mag$tax_table)


# Remove the columns that are metadata

test_nometa <- test |> 
  map(~select(.x ,"id":"optimumT", "mingentime"))
mt_mag$tax_table <- mt_mag$tax_table |> 
  rownames_to_column("id") |> 
  left_join(test_nometa$trait_matrixatgranularity1, by = "id") |> 
  left_join(test_nometa$trait_matrixatgranularity2, by = c("id", "genome_length", "optimumT", "mingentime")) |> 
  left_join(test_nometa$trait_matrixatgranularity3, by = c("id", "genome_length", "optimumT", "mingentime")) |> 
  left_join(test_nometa$hmm_matrix, by = c("id", "genome_length", "optimumT", "mingentime")) |> 
  left_join(test_nometa$rule_matrix, by = c("id", "genome_length", "optimumT", "mingentime")) |> 
  column_to_rownames("id")

view(test$hmm_matrix)
view(genomeset_results_norm$hmm_matrix)
view(mt_mag$tax_table)

mt_mag_large <- mt_mag

save(mt_mag_large, file = "output/data/mt_mag_large.RData")


# Okay from here let me start to look for autotrophs and such
library(microeco)
mt_mag_auto <- clone(mt_mag_large)


mt_mag_auto$tax_table <- mt_mag_auto$tax_table |> 
  select(d, p, c, o , f, g, s, bin, matches("g1_.*Resource.*")) |> 
  mutate(poss_chemotroph = if_any(matches("g1_.*Chemolith.*"), ~.x != 0), 
         poss_phototroph = if_any(matches("g1_.*Phototro.*"), ~.x != 0), 
         poss_heterotroph = if_any(matches("g1_.*organohetero.*"), ~.x != 0))


view(mt_mag_auto$tax_table)

# create a little df of those new columns

metabolism_df <- mt_mag_auto$tax_table |> 
  select(bin, matches("poss.*"))

#Add that back into the main microtable 

mt_mag_large$tax_table <- mt_mag_large$tax_table |> 
  left_join(metabolism_df, by = "bin")

# Turn these all into presence/absence

mt_mag_large$tax_table <- mt_mag_large$tax_table |> 
  mutate(across(matches(c("g1", "g2", "g3", "hmm", "rule")), ~ifelse(. == 0, FALSE, TRUE)))

mt_mag_large$cal_abund()



mt_mag_auto <- clone(mt_mag_large)

mt_mag_auto$tax_table |> 
  select(matches("g2")) |> 
  colnames()







# Normalize count traits


traits <- traits_listbygranularity[[3]] %>%
  dplyr::select(`microtrait_trait-name`) %>%
  # dplyr::filter(`microtrait_trait-name` != "Resource Use:Chemotrophy:chemolithoautotrophy:nitrite oxidation" ) %>%
  dplyr::pull(`microtrait_trait-name`)
summer_trait_matrixatgranularity3 <- genomeset_results_norm[["trait_matrixatgranularity3"]] %>%
  filter(strat_season_2 == "summer") |> 
  dplyr::select(c("id", traits, "mingentime", "optimumT")) %>%
  #dplyr::slice(1:100) %>%
  dplyr::filter(`mingentime` < 100) %>%   # max out mingentime at 25 days to avoid errors
  dplyr::filter(!is.na(`mingentime`) & !is.na(`optimumT`))

summer_matrixatgranularity3_bin <- summer_trait_matrixatgranularity3 |> 
  microtrait::trait.continuous2binary()

trait2trait_corr(summer_matrixatgranularity3_bin, verbose = TRUE, idcol = "id", outdir = "output/microtrait", dataset = "test")

summer_matrixatgranularity3_bin_log <- summer_matrixatgranularity3_bin |> 
  mutate(across(2:189, ~ifelse(.x == 1, TRUE, FALSE)))
genomset_distances <- summer_matrixatgranularity3_bin_log |> 
  microtrait::calc_mixeddist(idcol = "id", col2ignore = c("mingentime", "optimumT"), method = "wishart", binarytype = "logical", byrow = 1, verbose = TRUE)

# Fix the broken microtrait function here
compute.prevalence.mb <- function (feature_matrix, type, idcol = "id") 
{
    if (type == "hmm_matrix") {
        prevalence = feature_matrix %>% tibble::as_tibble() %>% 
            dplyr::select(c(idcol, intersect(colnames(feature_matrix), 
                hmms_fromrules %>% pull(`microtrait_hmm-name`)))) %>% 
            tidyr::pivot_longer(!id, names_to = "hmm", values_to = "presence") %>% 
            dplyr::select(-id) %>% dplyr::count(hmm, presence, 
            .drop = FALSE) %>% dplyr::group_by(hmm) %>% dplyr::mutate(percent = n/sum(n) * 
            100) %>% dplyr::filter(presence == 0) %>% dplyr::mutate(percent = 100 - 
            percent) %>% dplyr::select(c("hmm", "percent")) %>% 
            dplyr::inner_join(hmms_fromrules, by = c(hmm = "microtrait_hmm-name"))
    }
    if (type == "rule_matrix") {
        prevalence = feature_matrix %>% tibble::as_tibble() %>% 
            dplyr::select(c(idcol, rules %>% pull(`microtrait_rule-name`))) %>% 
            tidyr::pivot_longer(!id, names_to = "rule", values_to = "presence") %>% 
            dplyr::select(-id) %>% dplyr::count(rule, presence, 
            .drop = FALSE) %>% dplyr::group_by(rule) %>% dplyr::mutate(percent = n/sum(n) * 
            100) %>% dplyr::filter(presence == 0) %>% dplyr::mutate(percent = 100 - 
            percent) %>% dplyr::select(c("rule", "percent")) %>% 
            dplyr::inner_join(rules, by = c(rule = "microtrait_rule-name"))
    }
    if (type == "trait_matrixatgranularity3") {
        prevalence = feature_matrix %>% tibble::as_tibble() %>% 
            dplyr::select(c(idcol, intersect(colnames(feature_matrix), 
                traits_listbygranularity[[3]] %>% dplyr::pull(`microtrait_trait-name`)))) %>% 
            dplyr::mutate_at(vars(!starts_with("id")), funs(case_when(. >= 
                1 ~ factor(1, levels = c(0, 1)), TRUE ~ factor(0, 
                levels = c(0, 1))))) %>% tidyr::pivot_longer(!idcol, 
            names_to = "microtrait_trait-name", values_to = "presence") %>% 
            dplyr::select(-idcol) %>% dplyr::count(`microtrait_trait-name`, 
            presence, .drop = FALSE) %>% dplyr::group_by(`microtrait_trait-name`) %>% 
            dplyr::mutate(percent = n/sum(n) * 100) %>% dplyr::filter(presence == 
            1) %>% dplyr::select(c("microtrait_trait-name", "n", 
            "percent"))
    }
    prevalence
}

# Now run the prevalence
summer_prevalence <- compute.prevalence.mb(summer_trait_matrixatgranularity3, type="trait_matrixatgranularity3")


# Do the same thing but for the fall/winter samples
traits <- traits_listbygranularity[[3]] %>%
  dplyr::select(`microtrait_trait-name`) %>%
  # dplyr::filter(`microtrait_trait-name` != "Resource Use:Chemotrophy:chemolithoautotrophy:nitrite oxidation" ) %>%
  dplyr::pull(`microtrait_trait-name`)
fall_trait_matrixatgranularity3 <- genomeset_results_norm[["trait_matrixatgranularity3"]] %>%
  filter(strat_season_2 == "winter") |> 
  dplyr::select(c("id", traits, "mingentime", "optimumT")) %>%
  #dplyr::slice(1:100) %>%
  dplyr::filter(`mingentime` < 100) %>%   # max out mingentime at 25 days to avoid errors
  dplyr::filter(!is.na(`mingentime`) & !is.na(`optimumT`))

fall_matrixatgranularity3_bin <- fall_trait_matrixatgranularity3 |> 
  microtrait::trait.continuous2binary()

trait2trait_corr(fall_matrixatgranularity3_bin, verbose = TRUE, idcol = "id", outdir = "output/microtrait", dataset = "test")

fall_matrixatgranularity3_bin_log <- fall_matrixatgranularity3_bin |> 
  mutate(across(2:189, ~ifelse(.x == 1, TRUE, FALSE)))
genomset_distances <- fall_matrixatgranularity3_bin_log |> 
  microtrait::calc_mixeddist(idcol = "id", col2ignore = c("mingentime", "optimumT"), method = "wishart", binarytype = "logical", byrow = 1, verbose = TRUE)

# Fix the broken microtrait function here
compute.prevalence.mb <- function (feature_matrix, type, idcol = "id") 
{
    if (type == "hmm_matrix") {
        prevalence = feature_matrix %>% tibble::as_tibble() %>% 
            dplyr::select(c(idcol, intersect(colnames(feature_matrix), 
                hmms_fromrules %>% pull(`microtrait_hmm-name`)))) %>% 
            tidyr::pivot_longer(!id, names_to = "hmm", values_to = "presence") %>% 
            dplyr::select(-id) %>% dplyr::count(hmm, presence, 
            .drop = FALSE) %>% dplyr::group_by(hmm) %>% dplyr::mutate(percent = n/sum(n) * 
            100) %>% dplyr::filter(presence == 0) %>% dplyr::mutate(percent = 100 - 
            percent) %>% dplyr::select(c("hmm", "percent")) %>% 
            dplyr::inner_join(hmms_fromrules, by = c(hmm = "microtrait_hmm-name"))
    }
    if (type == "rule_matrix") {
        prevalence = feature_matrix %>% tibble::as_tibble() %>% 
            dplyr::select(c(idcol, rules %>% pull(`microtrait_rule-name`))) %>% 
            tidyr::pivot_longer(!id, names_to = "rule", values_to = "presence") %>% 
            dplyr::select(-id) %>% dplyr::count(rule, presence, 
            .drop = FALSE) %>% dplyr::group_by(rule) %>% dplyr::mutate(percent = n/sum(n) * 
            100) %>% dplyr::filter(presence == 0) %>% dplyr::mutate(percent = 100 - 
            percent) %>% dplyr::select(c("rule", "percent")) %>% 
            dplyr::inner_join(rules, by = c(rule = "microtrait_rule-name"))
    }
    if (type == "trait_matrixatgranularity3") {
        prevalence = feature_matrix %>% tibble::as_tibble() %>% 
            dplyr::select(c(idcol, intersect(colnames(feature_matrix), 
                traits_listbygranularity[[3]] %>% dplyr::pull(`microtrait_trait-name`)))) %>% 
            dplyr::mutate_at(vars(!starts_with("id")), funs(case_when(. >= 
                1 ~ factor(1, levels = c(0, 1)), TRUE ~ factor(0, 
                levels = c(0, 1))))) %>% tidyr::pivot_longer(!idcol, 
            names_to = "microtrait_trait-name", values_to = "presence") %>% 
            dplyr::select(-idcol) %>% dplyr::count(`microtrait_trait-name`, 
            presence, .drop = FALSE) %>% dplyr::group_by(`microtrait_trait-name`) %>% 
            dplyr::mutate(percent = n/sum(n) * 100) %>% dplyr::filter(presence == 
            1) %>% dplyr::select(c("microtrait_trait-name", "n", 
            "percent"))
    }
    prevalence
}

# Now run the prevalence
fall_prevalence <- compute.prevalence.mb(fall_trait_matrixatgranularity3, type="trait_matrixatgranularity3")


# OTHER STUFFFFFF
library(grid)
library(ComplexHeatmap)
library(tibble)
A4_ratio = 11.75/8.25
width = grid::unit(8.25*4, "inches")
height = grid::unit(8.25*4*A4_ratio, "inches")
cluster_traitmatrix(trait_matrix = matrixatgranularity3_bin,
                    idcol = "id", annot_cols = c("mingentime", "optimumT"), granularity = 3,
                    clustering_distance_rows = genomeset_distances, clustering_distance_cols = "binary",
                    width = width, height = height,
                    heatmap_width = width*0.8, heatmap_height = height*0.70,
                    row_dend_width = width*0.05, column_dend_height = height*0.02,
                    rightannotation_width = width*0.1, topannotation_height = height*0.02,
                    bottomannotation_height = height*0.005,
                    outdir = "output/microtrait", dataset = "test", pdf = TRUE)
# Test
base_dir = system.file("extdata/precomputed",package = "microtrait")
dataset = "environmentalgenomes"
test_genomeset_results = readRDS(file.path(base_dir, paste0(dataset, ".microtraitresults.rds")))
lapply(test_genomeset_results, dim)

library(dplyr)
genome_metadata = readRDS(file.path(base_dir, paste0(dataset, ".metadata.rds")))
genomeset_results_wmetadata = microtrait::add.metadata(test_genomeset_results, genome_metadata, genome_metadata_idcol = "IMG Taxon ID") |> convert_traitdatatype(binarytype = "logical")
# saveRDS(genomeset_results_wmetadata, file.path(base_dir, paste0(dataset, ".microtraitresults.rds")))


