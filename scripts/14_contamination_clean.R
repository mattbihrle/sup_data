# Contamination 

library(microeco)
library(tidyverse)

# MAG contamination -------------------------------------------------------

# Load in the MAG data

bac_df <- readr::read_tsv("data/gtdb_bin_tax/gtdbtk.bac120.summary.tsv", na = "N/A")
arc_df <- read_tsv("data/gtdb_bin_tax/gtdbtk.ar53.summary.tsv", na = "N/A")
all_bins_qual <- read_tsv("data/derep_bins/all_bins_quality_report.tsv") |> 
  rename_with(tolower)

mag_df <- bind_rows(bac_df, arc_df) |> 
  left_join(all_bins_qual, by = join_by(user_genome == name)) |> 
  mutate(sample = str_extract(user_genome, "WM.{2}"), .after = user_genome) |> 
    mutate(genome = str_remove(user_genome, "MAGScoT_cleanbin_000")) |> 
  separate_wider_delim(
  cols = classification, delim = ";", names = clade_cols, too_few = "align_start") |> 
  select(-user_genome)

# Create a list of contaminated samples
con_samp <- unique(mag_df$sample) |> 
  str_subset("WM(2.{1}|1[7-9]{1})")
con_samp
mag_contam <- mag_df |> 
  filter(sample %in% con_samp) |> 
  select(sample, genome, completeness, contamination, everything())

mag_contam <- mag_contam |> 
  write_tsv("output/data/mag_contamination.tsv")

# 16S contamination ------------------------------------------------------------------------------
# Import 16S data----------------------------------------------------------------------------
  # First the sample names
  names_df <- readr::read_tsv("data/16s/lotus3_out/final_sample_map_4_R.txt")
  names_vec <- setNames(names_df$SampleID, names_df$RealNames)


  # Load microeco object
    mt_16s <- phyloseq::import_biom("data/16s/lotus3_out/OTU.biom") |> 
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
      mutate(sample = str_extract(rownames(mt_16s$sample_table), ".{4}"))
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
    
# Look for OTUs only in the synmock 
    syn <- clone(mt_16s)

    syn_taxa <- syn$otu_table |> 
      select(matches("syn", ignore.case = TRUE)) |> 
      mutate(across(everything(), ~ifelse(.x == 0, NA, .x))) |> 
      drop_na() |> 
      rownames()

syn_taxa <- paste0("o__", syn_taxa)
head(syn_taxa)
syn$tax_table <- syn$tax_table |> 
  filter(otu %in% syn_taxa)

syn$tidy_dataset()
syn$cal_abund()
view(syn$taxa_abund$otu)

# syn$taxa_abund$otu <- syn$taxa_abund$otu |> 
#   rowwise() |> 
#   mutate(mean_abund = sum(c_across(starts_with("WM"))/19)) |> 
#   ungroup() |> 
#   select(matches("Syn"), mean_abund, everything())
# view(syn$taxa_abund$otu)

# Calculate test statistics and thngs
syn_table <- syn$otu_table |> 
  rownames_to_column("otu") |> 
  rowwise() |> 
  mutate(mean_abund = sum(c_across(starts_with("WM")))/19, 
        mean_abund_samp = sum(c_across(WM01_S2:WM16_S12))/11,
        mean_abund_contam = sum(c_across(WM17_S13:WM24_S20))/8)  |>
  mutate(
    t_test_diff = t.test(c_across(WM01_S2:WM16_S12), 
                        c_across(WM17_S13:WM24_S20), alternative = "g")$p.value
                      ) |> 
  ungroup() |> 

  select(otu, t_test_diff, matches("Syn"), matches("mean_abund"), everything())
view(syn_table)

# Now remove taxa where the mean abundance of reads across samples is larger than the reads
# in the synthetic community

syn_table <- syn_table |> 
  filter_out(mean_abund > `SynMock-Sheik_S1`)
syn_table |> 
  view()

# Now filter out taxa that had a significant t_test statistic (there is a significant difference
# between the real samples and the contamination samples)

syn_table <- syn_table |> 
  filter(t_test_diff > 0.05)
nrow(syn_table)

view(syn_table)

# Okay save this list as possible contaminants

syn_contam_table <- syn$tax_table |> 
  filter(otu %in% paste0("o__", syn_table$otu))

syn_contam <- syn_contam_table |> 
  select(otu) |> 
  as_vector() |> 
  str_remove("o__")

head(syn_contam)
# ----------------------------------------------------


# Now back to the original 
# Clean unneeded data
  # Remove 'Syn Mock' sample 
  mt_16s$sample_table <- mt_16s$sample_table |> 
    slice(
    str_which(mt_16s$sample_table$fastqFile, "Syn.*", negate = T)
      ) |> 
    rownames_to_column(var = "rownames") |> 
      # Add in the other metadata
    left_join(y = sw_meta) |> 
      # move rownames back
    column_to_rownames(var = "rownames")

  nrow(mt_16s$sample_table)
# Remove anything not Archaea or Bacteria
   mt_16s$tax_table <- mt_16s$tax_table |> 
    dplyr::slice(
    stringr::str_which(mt_16s$tax_table$k, ".*Bacteria|.*Archaea")
   )
# Check to be sure we removed them
    mt_16s$tax_table$k |> 
      unique()
  # Remove mitochondria and chloroplasts
  print("removing mitochontria and chloroplast")
  mt_16s
  mt_16s$filter_pollution(taxa = c("mitochondria", "chloroplast", "metagenome"))
  mt_16s      
mt_16s$tidy_dataset()
mt_16s
mt_16s_contam <- clone(mt_16s)
## Remove contamination ---------------------------------------------------------------
# Create vector of sites that were not pumped
  con_site_cols <- site_cols[12:19]

# Create df of all taxa found in the samples that weren't pumped to act as potential contamination
con_df <- mt_16s_contam$otu_table |> 
  select(!matches("WM01"): matches("WM16")) |> 
  mutate(across(everything(), ~ifelse(.x == 0.00, NA, .x))) |> 
  drop_na() |> 
  rownames_to_column("otu") |> 
  mutate(id = otu) |> 
  column_to_rownames("id")

# Double check that there are no dupliated species
nrow(con_df) == length(unique(con_df$otu))

# Great, now put all those taxa in a vector
con_vec <- unique(con_df$otu)

# Check to see what is similar between con_vec and syn_con


otu_contam <- c(con_vec, syn_contam)

length(otu_contam)

length(unique(otu_contam))

contam_16s <- unique(otu_contam)
# Remove taxa that are NOT contamination
mt_16s_contam$tax_table <- mt_16s_contam$tax_table |> 
  filter(otu %in% paste0("o__", contam_16s))

# Now remove the samples that DID filter
# mt_16s_contam$sample_table <- mt_16s_contam$sample_table |> 
#   filter(sample %in% con_site_cols)
# -------------------------------------------------------------------

  # Tidy dataset
  print("Tidying dataset")
  mt_16s_contam
  mt_16s_contam$tidy_dataset()
  mt_16s_contam

mt_16s_contam$sample_table$strat_season <- factor(mt_16s_contam$sample_table$strat_season, 
  levels = c("summer", "fall", "winter", "spring"), 
  labels = c("Summer", "Fall", "Winter", "Spring"),
   ordered = T)
mt_16s_contam$cal_abund()
save(mt_16s_contam, file = "output/data/mt_16s_contam.RData")

trans_abund$new(mt_16s_contam, taxrank = "otu", ntaxa = 30)$plot_bar()
  trans_abund$new(mt_16s_contam, taxrank = "otu", ntaxa = 50, high_level = "g")$plot_bar(ggnested = T)


mt_16s_contam <- mt_16s_contam$tax_table |> 
  write_tsv("output/data/16s_contam_taxa.tsv")
