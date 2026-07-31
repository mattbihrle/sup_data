
# List packages to load
packages <- c("microeco", "file2meco", "micRoclean", "tidyverse", "googlesheets4")
# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
# Packages loading
invisible(lapply(packages, library, character.only = TRUE))
# Remove the list of packages and installed packages from the environment
rm(installed_packages, packages)

# Import metadata ----------
load("output/data/sample_table.RData")
# Import 16S data----------------------------------------------------------------------------
  # First the sample names
  names_df <- readr::read_table("msi_downloads/lotus3_out/final_sample_map_4_R.txt")
  names_vec <- setNames(names_df$SampleID, names_df$RealNames)


  # Load microeco object
    mt_16s <- phyloseq::import_biom("msi_downloads/lotus3_out/OTU.biom") |> 
        phyloseq2meco()
      # Rename columns
      mt_16s$otu_table <- mt_16s$otu_table |> 
        rename(any_of(names_vec)) |> 
        # Remove the S__ from the end of the otu columns
        rename_with(~ str_extract(.x,".*(?=_S[0-9]{1,2})" ))
      # Reorder columns
      mt_16s$otu_table <- mt_16s$otu_table |> 
       dplyr::select(order(colnames(mt_16s$otu_table)))
      colnames(mt_16s$otu_table)
      # rename rows
    mt_16s$sample_table <- mt_16s$sample_table |> 
      rownames_to_column("sample_seq") |> 
      dplyr::left_join(names_df, by = join_by(sample_seq == SampleID)) |> 
        mutate(full_id = str_extract(RealNames, ".*(?=_S[0-9]{1,2})")) |> 
        mutate(rows = full_id) |> 
      select(-RealNames)
    
    # join the 16s sample table and the one created earlier

mt_16s$sample_table <- mt_16s$sample_table |> 
  left_join(sample_table, by = "full_id") |> 
  column_to_rownames("rows")
colnames(mt_16s$sample_table)
# Reorder rows alphabetically
      mt_16s$sample_table <- mt_16s$sample_table |> 
        slice(order(rownames(mt_16s$sample_table)))
    print(rownames(mt_16s$sample_table))
    # Rename Taxa
    taxa_names <- setNames( c("k", "p", "c", "o", "f", "g", "s"), colnames(mt_16s$tax_table))
    colnames(mt_16s$tax_table) <- taxa_names
    # Add an "otu" column for later use
    mt_16s$tax_table <- mt_16s$tax_table |> 
      mutate(otu = rownames(mt_16s$tax_table))
    # check to make sure it worked
    rownames(mt_16s$tax_table) == mt_16s$tax_table$otu
    colnames(mt_16s$tax_table) 
  # Replace question marks with 'unknown'
    mt_16s |> tidy_taxonomy(pattern = "\\?", replacement = "unknown")
  # Replace blanks with "unknown"
    mt_16s$tax_table <- mt_16s$tax_table |> 
      mutate(across(everything(), ~str_replace(.x, "__$", "__unknown")))
    mt_16s$tax_table[1:5,]

  # Remove mitochondria and chloroplasts
  print("removing mitochontria and chloroplast")
  mt_16s
  mt_16s$filter_pollution(taxa = c("mitochondria", "chloroplast", "metagenome"))
  mt_16s$tidy_dataset()
  mt_16s     

  # Remove anything not Archaea or Bacteria
  mt_16s$tax_table <- mt_16s$tax_table |> 
  dplyr::slice(
    stringr::str_which(mt_16s$tax_table$k, ".*Bacteria|.*Archaea")
  )
# Check to be sure we removed them
  mt_16s$tax_table$k |> 
    unique()
mt_16s
mt_16s$tidy_dataset()
mt_16s

# Remove 'Syn Mock' sample 
mt_16s$sample_table <- mt_16s$sample_table |> 
  slice(
    str_which(mt_16s$sample_table$fastqFile, "Syn.*", negate = T)
  )

mt_16s$tidy_dataset()
# This should also remove OTUs only found in the synthetic community
mt_16s

mt_16s$tax_table <- mt_16s$tax_table |> 
  mutate(otu_row = rownames(mt_16s$tax_table)) |> 
  pivot_longer(cols = c(k, p, c, o, f, g, s), names_to = "rank", values_to = "val") |> 
  mutate(known_val = if_else(str_detect(val, "unknown"), NA_character_, val)) |> 
  group_by(otu_row) |> 
  fill(known_val, .direction = "down") |> 
  mutate(val = if_else(str_detect(val, "unknown") & !is.na(known_val), 
          paste0(val, "_", known_val), val)) |> 
  select(-known_val) |> 
  pivot_wider(names_from = rank, values_from = val) |> 
  ungroup() |> 
  column_to_rownames("otu_row")

# Label all blanks as unknowns
# # set the strat_season to be a factor
# mt_16s$sample_table$strat_season <- factor(mt_16s$sample_table$strat_season, 
#                                            levels = c("summer", "fall", "winter", "spring"), 
#                                            labels = c("Summer", "Fall", "Winter", "Spring"),
#                                            ordered = T)


# Okay here start working with the micRoclean package for contamination------------------------------------------------------------

# Now add the control and sample type
control_samps <- mt_16s$sample_table$full_id[which(mt_16s$sample_table$pumped == FALSE)]

mt_16s$sample_table <- mt_16s$sample_table |> 
  mutate(is_control = ifelse(full_id %in% control_samps, TRUE, FALSE), 
         sample_type = ifelse(is_control == TRUE, "blank", "DNA")
  )

# mt_16s$sample_table <- mt_16s$sample_table |> 
#   mutate( 
#          is_control = ifelse(sample == "SynM", TRUE, FALSE),
#         sample_type = ifelse(is_control == TRUE, "blank", "DNA"), 
#       batch = "a")
# Now extract that as the required metadata file

meta <- mt_16s$sample_table |> 
  dplyr::select(sample_type, is_control, sample_well, batch) |>
  rownames_to_column("full_id") |> 
  arrange(batch, desc(is_control), full_id) |> 
  column_to_rownames("full_id")



head(meta)
# Okay, now pull out the count data

count <- mt_16s$otu_table |> 
  t() |> 
  as.data.frame() |>
  rownames_to_column("full_id") |>
  slice(match(rownames(meta), full_id)) |>
  column_to_rownames("full_id") 
all(rownames(count) == rownames(meta))

head(count)
# Run it!

mclean_results <-  micRoclean(counts = count, 
                              meta = meta, 
                              research_goal = 'orig.composition', 
                              control_name = control_samps)

mclean_results$blank <- "all"
paste("filtering loss:",mclean_results$filtering_loss) |> 
  print()

# remove the taxa that are all 0s 
mclean_results$decontaminated_count <- 
  mclean_results$decontaminated_count |> 
  as.data.frame() |> 
  setNames(colnames(count)) |> 
  dplyr::select(where(~ !all(.x == 0)))

# Filter to keep only the taxa in the "decontaminated count" table
  # First save the mt_16s before doing any filtering
full_mt_16s <- clone(mt_16s)
mt_16s$otu_table <- mt_16s$otu_table |> 
  filter(rownames(mt_16s$otu_table) %in% colnames(mclean_results$decontaminated_count))

# Then remove the blank columns (WM17-WM24)
mt_16s$sample_table <- mt_16s$sample_table |>
  filter_out(full_id %in% control_samps)
# Verify that we have 11 samples and less OTUs in otu_table
mt_16s
#Tidy dataset
mt_16s$tidy_dataset()
mt_16s
save(mt_16s, file = "output/data/mt_16s.RData")

# -------------------------------------------------------------------------------
# Import list of contaminants from Sheik et al 2018 Supplemental Material and from Eisenhofer 2019
load("output/data/mt_16s.RData")
mt_16s
contam_genus <- read_tsv("sheik_2018_contamination_table.txt") |>
  rename_with(tolower) |>  
  filter_out(extraction == FALSE & water == FALSE) |> 
  dplyr::select(g) |> 
  as_vector()

e_contam <- read_tsv("eisenhofer_contam.txt") |> 
as_vector() |> 
  sort()

contam_genus <- c(contam_genus, e_contam) |> 
  unique()

# Label which taxa are found in that list

mt_16s$tax_table <- mt_16s$tax_table |> 
  mutate(sheik_contam = ifelse(g %in% paste0("g__", contam_genus), TRUE, FALSE))

# view(mt_16s$tax_table)



# hmm okay I want to try and plot how the abundance of the sheik contaminents changes
# n = 100
# rank = "g"
# abund <- trans_abund$new(mt_16s, ntaxa = n, taxrank = rank, high_level = "sheik_contam")$
#   plot_bar(others_color = "grey70", ggnested = T, xtext_keep = TRUE, legend_text_italic = FALSE, 
#            facet = "strat_season") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=9)) +
#   ggtitle("16S Data", subtitle = paste(n, "most abundant", rank, "brown, possible contaminent from 16s data,", 
#                                        "green, not possible contaminent")) +
#   theme(legend.position = "none",
#         # axis.text.x = element_blank())
# abund
# plotly::plotly_build(abund)

# Okay now I want to look at just the taxa that are sheik_contams

sheik_16s <- clone(mt_16s)
sheik_16s$tax_table <- sheik_16s$tax_table |> 
  filter(sheik_contam == TRUE)

sheik_16s$tidy_dataset()
sheik_16s$cal_abund()

# n = 66
# rank = "g"
# abund <- trans_abund$new(sheik_16s, ntaxa = n, taxrank = rank)$
#   plot_bar(others_color = "grey70", xtext_keep = TRUE, legend_text_italic = FALSE, 
#            facet = "strat_season") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=9))
# # ggtitle("16S Data", subtitle = paste(n, "most abundant", rank, "brown, possible contaminent from 16s data,", 
# # "green, not possible contaminent")) +
# # theme(legend.position = "none")
# abund
# plotly::plotly_build(abund)
genus <- sheik_16s$tax_table |> 
  mutate(across(everything(), ~str_remove(.x, "[a-z]{1}__"))) |> 
  arrange(g)
  # slice(unique(str_which(g, "g__Caulobacter")):nrow(sheik_16s$tax_table))
# view(genus)
# sheik_16s$tax_table |>
#   filter(g == "g__Bacillus") |> 
#   arrange(s) |> 
#   view()


test <- clone(mt_16s)
test$tax_table <- test$tax_table |>
  filter_out(g == "g__unknown")
test$cal_abund()
trans_abund$new(test, ntaxa = 50, taxrank = "s")$
  plot_bar()
test$cal_abund()

geni <- trans_abund$new(test, ntaxa = 630 , taxrank = "g")$
  plot_bar()
geni$data$Taxonomy
top_tax <- geni$data |> 
  select(Taxonomy) |> 
  distinct()

top_tax <- top_tax |> 
  mutate(tax_rank = 630 - as.numeric(rownames(top_tax))) |> 
  mutate(tax = as.character(Taxonomy)) |> 
  arrange(tax)

genus <- genus |> 
  left_join(top_tax, by = join_by("g" == "tax"))
  
# trans_abund$new(mt_16s, ntaxa = 200, taxrank = "s")$
#   plot_bar() |> 
#   plotly::plotly_build()


# trans_abund$new(test, ntaxa = 200, taxrank = "g")$
#   plot_bar() |> 
#   plotly::plotly_build()
# Hand removing taxa based on the rules:
# Found in either the Sheik contamination or Eisenhofer contamination list AND at least one of:
# 1) If it is a "small" abundance (not included in top 200 genuses by mean abundance or determined as "small" by visual inspection)
# 2) If the genus or species are not recorded to be found in water samples
# 3) Based on other conclusions from quick lit search and relative abundance (will be explain in text)
# When in doubt, will air on the side of removing more

bottom_tax <- top_tax |> 
  filter(tax_rank >= 200)

mt_16s$tax_table <- mt_16s$tax_table |> 
  # filter out contaminants that are not in the top 200 genuses
  filter_out(sheik_contam == TRUE & g %in% paste0("g__", bottom_tax$tax)) |> 
  # Acidiovorax seems to be primarily a plant genus so removing that
  # Remove Acidiovorax from the tax_table
  filter(g != "g__Acidiovorax") |> 
  # Achromobacter
  filter( g != "g__Achromobacter") |> 
  filter(g != "g__Acinetobacter") |> 
  # Keep aeromicrobium in because it is found in marine systems
  # Bacillus seems to be primarily human-related so should remove those as well
  filter(g != "g__Bacillus") |> 
  # Get ride of afipia
  filter(g != "g__Afipia") |> 
  # Some of the aquabacterium seem to be a mix of human derived and other. Because they already have a small abundance,
  # I will remove them
  filter(g != "g__Aquabacterium") |>
  # Bacillus mar... found in human stool removing that.
  filter(s != "s__Bacillus marasmi") |> 
  # As far as the other bacillus they seem to some found in water and some human related.
  # I will take out each species that is not water related. 
  # I cannot find any information that cohnii is water related
  filter(s != "s__Bacillus cohnii") |> 
  # I am going to keep Bacillus cerus in because it is found all over but I am 
  # keeping an eye on it because it could also be a contaminant
  # Bosea is found in freshwater so I'll keep it in.
  # Brevibacillus is also in water so I'll keep it in
  # Looks like bradyrhyzobium is primarily soil and nitrogen fixing, I think this would be one to get rid of as well
  filter(g != "g__Bradyrhizobium") |> 
  # Brevundimonas seems to be pretty ubiquitous so I will let it stay, it's pretty small anyways
  # Brevebacterium is found in food and human skin so I'll take it out
  filter(g != "g__Brevibacterium") |>
  # Caulobacter is found in FW, keeping it in
  # comamonas is found in freshwater, keeping it in
  # g__Chryseobacterium is found in freshwater and isn't huge so I'll keep it in
  # Corynebacterium is found on humans but the two species identified are found in marine environments so I think I should keep them
  # Cupriavidas is found all over so might be a contaminant, removing it
  filter(g != "g__Cupriavidus") |> 
  # Curvibacter are from well water so keep them in
  # s__Deinococcus yunweiensis was originally found as a contaminant, remove it and the rest. Their defining characteristic is being resistant to UV radiation which makes me think they just resisted the UV, cleaning and freezing
  filter(g != "g__Deinococcus") |> 
  # g__Devosia is mostly soil bacteria, remove this as well
  filter(g != "g__Devosia") |> 
  # Dietzia looks to be mostly human related, removing that as well
  filter(g != "g__Dietzia") |> 
  # Dyadobacteria is such a small abundance I am going to remove it
  filter(g != "g__Dyadobacter") |> 
  # flavobacteriums seem to be all over in freshwater, keeping in the ones that are "unknown" for species
  # Looks like sp VMW seems to be originally from freshwater so I'll keep it
  # s__Flavobacterium branchiophilum creates a gill disease in fish so likely is a 'real' microbe
  #s__Flavobacterium swingsii originally isolated from a river so keep it in
  # s__antartic bacterium was isolated from soil and isn't very abundant. going to remove it
  filter(s != "s__Antarctic bacterium") |> 
  # Hydrotalea has been found in aquatic environments so I will keep it in
  # Janthinobacteria are found in freshwater so I will keep them in
  # Leptothrix found in freshwater so will keep
  # limnobacter seems to be found in lake sediments so I'll keep it in
  # Massilia are found in freshwater so I should keep them in, but one to keep an eye on
  # Mesorhizobium are found in freshwater so I'll keep them in
  # Microbacterium is not found once I swap to the new kits so I am going to call it a contaminant and remove it
  filter(g != "g__Microbacterium") |> 
  filter(g != "g__Microlunatus") |> 
  filter(g != "g__Nevskia") |> 
  # Novosphingo seem to be pretty much just contaminents removing anything that is unknown
  filter_out(g == "g__Novosphingobium" & s == "s__unknown") |> 
  # g__Paenibacillus has so little abundance and seems to be more soil related
  filter(g != "g__Paenibacillus") |> 
  # Patulibacter are mostly just soil and wastewater samples
  filter(g != "g__Patulibacter") |> 
  # pedobacter is just a contaminent
  filter(g != "g__Pedobacter") |> 
  # I am going to keep Polaromonas in as it is found all over the place and is found in most of my samples but in varying proportions
# pseudomonas seems to be all over the place in water and in soil. and it is most of the samples. I will keep it in. 
  # Ralstonia is found in just the WM2122 year. That's odd, I'll keep it in. It is primarily a bacterial pathogen which is interesting
  # Rhodococcus is found in water and doesn't have any suspicious patterns. Keep it in
  # Roseomonas is primarily found in a group of samples that was extracted at the same time. Of those taxa, 
  # the most common is one that is found in drinking water suggesting it isn't really a freshwater taxa, plus it is not very abundant overall. I am going to remove it. 
  filter(g != "g__Roseomonas") 
  # Sediminibacterium. This one has a lot of OTUs and is quite abundant. The candidatus Jacksonbacteria
    # is found in some water samples. I will keep it in. The goheungnse species was originally 
    # cultured from a freshwater reservior so keep that in as well. For the other "unknown" species I might as well keep them in too.
  # Sphingobium are found everywhere. I will keep them in
  # Sphingomonas is found in all the samples and species of it are found in freshwater. I am a little worried abou the pattern where they are most abundnant in the samples that are oldest suggesting that they are gradually taking over.
    # I might run it by Cody and see what he thinks. Looking at the breakdown of OTUs it seems like most of the abundance comes from OTU4 and OTU82. I should look more into these to see what they are similar to. 
    # For now keeping sphinogomonas in
  # Sphynogopyxis have been found in freshwater and are pretty abundant. Keeping them in
  # Stenotrophomonas is found all over and is most abundant in the Lake Erie sample so it seems to be fine. 
  # Sulfuritaleta was found in freshwater lakes, keep it in
  # unidbacterium is around in alpine lakes so maybe I will keep them in as well
  # Variovorax are found in freshwater so keep them in too

# Okay now that I've gone through all of that, tidy mt_16s once more

mt_16s$tidy_dataset()
mt_16s

# Create a tax table of just the taxa that were removed from the full dataset 

contam_16s <- clone(full_mt_16s)
contam_16s
contam_16s$tax_table <- contam_16s$tax_table |> 
  filter_out(otu %in% mt_16s$tax_table$otu)
contam_16s
contam_16s$tidy_dataset()
contam_16s
save(contam_16s, file = "output/data/contam_16s.RData")
save(full_mt_16s, file = "output/data/full_16s.RData")
## Calculate relative abunance and diversity metrics-----------------------------------------

mt_16s$cal_abund()
mt_16s$taxa_abund$p[1:5,1:5]

# Calc Alpha Diversity
mt_16s$cal_alphadiv()

# Calc Bray Curtis Dissimilarity
mt_16s$cal_betadiv(method = c("bray","aitchison","robust.aitchison", "jaccard"))
mt_16s$beta_diversity$bray


save(mt_16s, file = "output/data/mt_16s.RData")

obs <- ls()
rm(list = obs)
rm(obs)

