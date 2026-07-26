# Load libraries
library(microeco)
library(file2meco) # for importing 16S phyloseq objects
library(micRoclean) # for cleaning out contamination
library(tidyverse)
library(googlesheets4)

# Import meta data from 2024-2025 year -----------------------------------------
sw_meta <- read_csv("output/data/metadata_supwinter.csv") |> 
  dplyr::select(!matches(c("12h_med", "round_date")))|> 
  distinct() |> 
  mutate(full_id = paste0("s_2425_", sample))

# Import master spreadsheet from google
# gs4_auth()
ext_sheet <- read_sheet("https://docs.google.com/spreadsheets/d/1Iq5AIxiH-EBa2lT83cKCUsLwYCvvCl0iHlCbkJw9rQo/edit?gid=1015683691#gid=1015683691")
# Create a sample table for import into the microtable
sample_table <- sw_meta |> 
  full_join(ext_sheet, by = "full_id") |> 
  mutate(rows = full_id) |> 
    column_to_rownames("rows") |> 
    as.data.frame()

sample_table <- sample_table |> 
  mutate(strat_season = fct_reorder(strat_season, date, .na_rm = F), 
         strat_season_2 = fct_reorder(strat_season_2, date, .na_rm = F),
         solar_season = fct_reorder(solar_season, date, .na_rm = F), 
         mixing = fct_relevel(mixing, "stratified_std", "mixed", "stratified_inverse"))


# Import 16S data----------------------------------------------------------------------------
  # First the sample names
  names_df <- readr::read_table("msi_downloads/lotus3_out/final_sample_map_4_R.txt")
  names_vec <- setNames(names_df$SampleID, names_df$RealNames)


  # Load microeco object
    mt_16s <- phyloseq::import_biom("msi_downloads/lotus3_out/OTU.biom") |> 
        phyloseq2meco()
      # Rename columns
      mt_16s$otu_table <- mt_16s$otu_table |> 
        rename(any_of(names_vec))
      # Reorder columns
      mt_16s$otu_table <- mt_16s$otu_table |> 
       dplyr::select(order(colnames(mt_16s$otu_table)))
      colnames(mt_16s$otu_table)
      # rename rows
    mt_16s$sample_table <- mt_16s$sample_table |> 
      rownames_to_column("sample") |> 
      dplyr::left_join(names_df, by = join_by(sample == SampleID)) |> 
      column_to_rownames("RealNames") 
    
    # Create a column of sample names for later metadata
    mt_16s$sample_table <- mt_16s$sample_table|> 
      mutate(full_id = str_extract(rownames(mt_16s$sample_table), ".*(?=_S[0-9]{1,2})"))
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
  ) |> 
  rownames_to_column(var = "rownames") |> 
  # Add in the other metadata
  left_join(y = sample_table, by = "full_id") |> 
  # move rownames back
  column_to_rownames(var = "rownames")

mt_16s$tidy_dataset()
# This should also remove OTUs only found in the synthetic community
mt_16s
# set the strat_season to be a factor
mt_16s$sample_table$strat_season <- factor(mt_16s$sample_table$strat_season, 
                                           levels = c("summer", "fall", "winter", "spring"), 
                                           labels = c("Summer", "Fall", "Winter", "Spring"),
                                           ordered = T)


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
# save(mt_16s, file = "output/data/mt_16s.RData")
# stop()
# -------------------------------------------------------------------------------
# Import list of contaminants from Sheik et al 2018 Supplemental Material
load("output/data/mt_16s.RData")
contam_genus <- read_tsv("sheik_2018_contamination_table.txt") |> 
  rename_with(tolower) |>  
  filter_out(extraction == FALSE & water == FALSE) |> 
  dplyr::select(g) |> 
  as_vector()

# Label which taxa are found in that list

mt_16s$tax_table <- mt_16s$tax_table |> 
  mutate(sheik_contam = ifelse(g %in% paste0("g__", contam_genus), TRUE, FALSE))

view(mt_16s$tax_table)



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
  arrange(g)

sheik_16s$tax_table |>
  filter(g == "g__Bacillus") |> 
  arrange(s) |> 
  view()


test <- clone(mt_16s)
test$tax_table <- test$tax_table |> 
  filter_out(s == "s__unknown")

test$cal_abund()

trans_abund$new(test, ntaxa = 100, taxrank = "s")$
  plot_bar() |> 
  plotly::plotly_build()

trans_abund$new(mt_16s, ntaxa = 100, taxrank = "s")$
  plot_bar() |> 
  plotly::plotly_build()
# Acidiovorax seems to be primarily a plant genus so removing that
# Remove Acidiovorax from the tax_table
mt_16s$tax_table <- mt_16s$tax_table %>%
  filter(g != "g__Acidiovorax") |> 
  # Bacillus seems to be primarily human-related so should remove those as well
  filter(g != "g__Bacillus") |> 
  # Get ride of afipia
  filter(g != "g__Afipia") |> 
  # Some of the aquabacterium seem to be a mix of human derived and other. Because they already have a small abundance,
  # I will remove them
  filter(g != "g__Aquabacterium") |>
  # Of the artrhobacter genus one species is particularly found in blood, removing that
  filter(s != "s__Arthrobacter woluwensis") |> 
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
  # START HERE
  # Looks like bradyrhyzobium is primarily soil and nitrogen fixing, I think this would be one to get rid of as well
  filter(g != "g__Bradyrhizobium") |> 
  # Brevundimonas seems to be pretty ubiquidus so I will let it stay, it's pretty small anyways
  # g__Chryseobacterium is found in freshwater and isn't huge so I'll keep it in
  # s__Deinococcus yunweiensis was originally found as a contaminant, remove it!
  filter(g != "g__Deinococcus") |> 
  # g__Devosia is mostly soil bacteria, remove this as well
  filter(g != "g__Devosia") |> 
  # flavobacteriums seem to be all over in freshwater, keeping in the ones that are "unknown" for species
  # Looks like sp VMW seems to be originally from freshwater so I'll keep it
  # s__Flavobacterium branchiophilum creates a gill disease in fish so likely is a 'real' microbe
  #s__Flavobacterium swingsii originally isolated from a river so keep it in
  # s__antartic bacterium was isolated from soil and isn't very abundant. going to remove it
  filter(s != "s__Antarctic bacterium") |> 
  # limnobacter seems to be found in lake sediments so I'll keep it in
  # Novosphingo seem to be pretty much just contaminents removing anything that is unknown
  filter(g != "g__Novosphingobacter") |> 
  # g__Paenibacillus has so little abundance and seems to be more soil related
  filter(g != "g__Paenibacillus") |> 
  # pedobacter is just a contaminent
  filter(g != "g__Pedobacter") |> 
  # Polaromonas is around in just that one WM03 sample. Makes me think it is an error
  filter(g != "g__Polaromonas") 
# pseudomonas seems to be all over the place in water and in soil. and it is most of the samples. I will keep it in. 
# unidbacterium is around in alpine lakes so maybe I will keep them in as well


# Okay now that I've gone through all of that, tidy mt_16s once more

mt_16s$tidy_dataset()
mt_16s
stop()
# Create a tax table of just the taxa that were removed from the full dataset 

contam_16s <- clone(full_mt_16s)
contam_16s
contam_16s$tax_table <- contam_16s$tax_table |> 
  filter_out(otu %in% mt_16s$tax_table$otu)
contam_16s
contam_16s$tidy_dataset()
save(contam_16s, file = "output/data/contam_16s.RData")
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
