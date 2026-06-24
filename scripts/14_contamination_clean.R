# Contamination 
# devtools::install_github("rachelgriffard/micRoclean")
# devtools::load_all("packages/micRoclean-main/micRoclean-main/micRoclean.Rproj")
# devtools::load_all("packages/SCRuB-Latest/SCRuB-Latest/SCRuB.Rproj")
# devtools::inst
library(microeco)
library(file2meco)
library(tidyverse)
library(micRoclean)
# MAG contamination -------------------------------------------------------

# # Load in the MAG data

# bac_df <- readr::read_tsv("data/gtdb_bin_tax/gtdbtk.bac120.summary.tsv", na = "N/A")
# arc_df <- read_tsv("data/gtdb_bin_tax/gtdbtk.ar53.summary.tsv", na = "N/A")
# all_bins_qual <- read_tsv("data/derep_bins/all_bins_quality_report.tsv") |> 
#   rename_with(tolower)

# mag_df <- bind_rows(bac_df, arc_df) |> 
#   left_join(all_bins_qual, by = join_by(user_genome == name)) |> 
#   mutate(sample = str_extract(user_genome, "WM.{2}"), .after = user_genome) |> 
#     mutate(genome = str_remove(user_genome, "MAGScoT_cleanbin_000")) |> 
#   separate_wider_delim(
#   cols = classification, delim = ";", names = clade_cols, too_few = "align_start") |> 
#   select(-user_genome)

# # Create a list of contaminated samples
# con_samp <- unique(mag_df$sample) |> 
#   str_subset("WM(2.{1}|1[7-9]{1})")
# con_samp
# mag_contam <- mag_df |> 
#   filter(sample %in% con_samp) |> 
#   select(sample, genome, completeness, contamination, everything())

# mag_contam <- mag_contam |> 
#   write_tsv("output/data/mag_contamination.tsv")

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

  # Remove mitochondria and chloroplasts
  print("removing mitochontria and chloroplast")
  mt_16s
  mt_16s$filter_pollution(taxa = c("mitochondria", "chloroplast", "metagenome"))
  mt_16s      

# Okay here start working with the micRoclean package ------------------------------------------------------------
      # First thing is to create a sequence of well numbers
wells <- c(LETTERS[1:8], LETTERS[1:8], LETTERS[1:4])
wells_num <- rep(1:3, each = 8)[1:length(wells)] |> 
  str_pad(width = 2, side = "left", pad = 0)
wells <- paste0(wells, wells_num)

wells <- c(
  "D03",  # Position 1: SynMock control
  paste0(rep(LETTERS[1:8], 2), rep(sprintf("%02d", 1:2), each = 8)),  # Positions 2-17
  paste0(LETTERS[1:3], "03")  # Positions 18-20
)
# Now add that to the mt_16s$sample_table with the D4 being the synmock
mt_16s$sample_table <- mt_16s$sample_table |> 
  mutate(sample_well = wells)

# Now add the control and sample type
 control_samps <- mt_16s$sample_table$sample[13:20]

mt_16s$sample_table <- mt_16s$sample_table |> 
  mutate(is_control = ifelse(sample %in% control_samps, TRUE, FALSE), 
         is_control = ifelse(sample == "SynM", TRUE, is_control),
        sample_type = ifelse(is_control == TRUE, "blank", "DNA"), 
        batch = "a")


# mt_16s$sample_table <- mt_16s$sample_table |> 
#   mutate( 
#          is_control = ifelse(sample == "SynM", TRUE, FALSE),
#         sample_type = ifelse(is_control == TRUE, "blank", "DNA"), 
#       batch = "a")
# Now extract that as the required metadata file

meta <- mt_16s$sample_table |> 
  select(sample_type, is_control, sample_well, batch) |>
  rownames_to_column("sample") |> 
  arrange(desc(is_control), sample) |> 
  column_to_rownames("sample")



head(meta)
# Okay, now pull out the count data

count <- mt_16s$otu_table |> 
  t() |> 
  as.data.frame() |>
  rownames_to_column("sample") |>
  slice(match(rownames(meta), sample)) |>
  column_to_rownames("sample") 
all(rownames(count) == rownames(meta))

head(count)
# Run it!
control_names <- c(rownames(count)[1], rownames(count)[13:20])

mclean_results <-  micRoclean(counts = count, 
      meta = meta, 
      research_goal = 'orig.composition', 
      control_name = rownames(count)[1:9])
mclean_results$blank <- "all"

control_names <- c(rownames(count)[13:20])
new_count <- count |> 
  slice(-1)
new_meta <- meta |> 
  slice(-1)
mclean_results_no_syn <-  micRoclean(counts = new_count, 
      meta = new_meta, 
      research_goal = 'orig.composition', 
      control_name = rownames(count)[1:9])
mclean_results_no_syn$blank <- "all_no_syn"

# remove the taxa that are all 0s 
  mclean_results_no_syn$decontaminated_count <- 
    mclean_results_no_syn$decontaminated_count |> 
    as_data_frame() |> 
    setNames(colnames(count)) |> 
    select(where(~ !all(.x == 0)))
# Here I want to set it up so it will run through the control samps one by one
# so there is not so much contamination taken out. 

results_list <- list()

for(i in 1:8){
  counts_new <- count |> 
    t() |> 
    as_data_frame() |> 
    select(control_names[1], matches(control_samps[i]), !matches(control_samps)) |> 
    t()

  meta_new <- meta |> 
    rownames_to_column("sample") |> 
      filter(sample %in% rownames(counts_new)) |> 
        column_to_rownames("sample")
      results_list[[i]] <-
        micRoclean(counts = counts_new, 
      meta = meta_new, 
      research_goal = 'orig.composition', 
      control_name = rownames(count_new)[1:2])

# List what the blank was
  results_list[[i]]$blank <- c(control_samps[i], "SynMock")
# Filter out the taxa that are all zeros across the board
  results_list[[i]]$decontaminated_count <- 
    results_list[[i]]$decontaminated_count |> 
    as_data_frame() |> 
    setNames(colnames(count)) |> 
    select(where(~ !all(.x == 0)))

}

for(i in 1:8){
  counts_new <- count |> 
    t() |> 
    as_data_frame() |> 
    select(matches(control_samps[i]), !matches(control_samps), -control_names[1]) |> 
    t()

  meta_new <- meta |> 
    rownames_to_column("sample") |> 
      filter(sample %in% rownames(counts_new)) |> 
        column_to_rownames("sample")
      results_list[[i+8]] <-
        micRoclean(counts = counts_new, 
      meta = meta_new, 
      research_goal = 'orig.composition', 
      control_name = rownames(count_new)[1:2])

# List what the blank was
  results_list[[i+8]]$blank <- c(control_samps[i])
# Filter out the taxa that are all zeros across the board
  results_list[[i+8]]$decontaminated_count <- 
    results_list[[i+8]]$decontaminated_count |> 
    as_data_frame() |> 
    setNames(colnames(count)) |> 
    select(where(~ !all(.x == 0)))

}
# Create new rows where TRUE is that the taxa was not removed when the sample in the column title
# was set as the blank
mt_16s$tax_table <- mt_16s$tax_table |> 
  mutate(
    wm17_blank = ifelse(otu %in% paste0("o__", colnames(results_list[[1]]$decontaminated_count)), TRUE, FALSE),
    wm18_blank = ifelse(otu %in% paste0("o__", colnames(results_list[[2]]$decontaminated_count)), TRUE, FALSE),
    wm19_blank = ifelse(otu %in% paste0("o__", colnames(results_list[[3]]$decontaminated_count)), TRUE, FALSE),
    wm20_blank = ifelse(otu %in% paste0("o__", colnames(results_list[[4]]$decontaminated_count)), TRUE, FALSE),
    wm21_blank = ifelse(otu %in% paste0("o__", colnames(results_list[[5]]$decontaminated_count)), TRUE, FALSE),
    wm22_blank = ifelse(otu %in% paste0("o__", colnames(results_list[[6]]$decontaminated_count)), TRUE, FALSE),
    wm23_blank = ifelse(otu %in% paste0("o__", colnames(results_list[[7]]$decontaminated_count)), TRUE, FALSE),
    wm24_blank = ifelse(otu %in% paste0("o__", colnames(results_list[[8]]$decontaminated_count)), TRUE, FALSE),
    wm17_blank_no_syn = ifelse(otu %in% paste0("o__", colnames(results_list[[9]]$decontaminated_count)), TRUE, FALSE),
    wm18_blank_no_syn = ifelse(otu %in% paste0("o__", colnames(results_list[[10]]$decontaminated_count)), TRUE, FALSE),
    wm19_blank_no_syn = ifelse(otu %in% paste0("o__", colnames(results_list[[11]]$decontaminated_count)), TRUE, FALSE),
    wm20_blank_no_syn = ifelse(otu %in% paste0("o__", colnames(results_list[[12]]$decontaminated_count)), TRUE, FALSE),
    wm21_blank_no_syn = ifelse(otu %in% paste0("o__", colnames(results_list[[13]]$decontaminated_count)), TRUE, FALSE),
    wm22_blank_no_syn = ifelse(otu %in% paste0("o__", colnames(results_list[[14]]$decontaminated_count)), TRUE, FALSE),
    wm23_blank_no_syn = ifelse(otu %in% paste0("o__", colnames(results_list[[15]]$decontaminated_count)), TRUE, FALSE),
    wm24_blank_no_syn = ifelse(otu %in% paste0("o__", colnames(results_list[[16]]$decontaminated_count)), TRUE, FALSE),
    all_blank = ifelse(otu %in% paste0("o__", colnames(mclean_results$decontaminated_count)), TRUE, FALSE)
  )

# Okay now filter out all the taxa that are FALSE for all of the blank samples
clean_mt_16s <- clone(mt_16s)

clean_mt_16s$tax_table <- clean_mt_16s$tax_table |> 
  filter(if_any(matches("blank")) == TRUE)

# Import list of contaminants from Sheik et al 2018 Supplemental Material

contam_genus <- read_tsv("sheik_2018_contamination_table.txt") |> 
  rename_with(tolower) |>  
  filter_out(extraction == FALSE & water == FALSE) |> 
  select(g) |> 
  as_vector()

clean_mt_16s
clean_mt_16s$tax_table <- clean_mt_16s$tax_table |> 
  filter_out(g %in% paste0("g__", contam_genus))
clean_mt_16s


# Okay now finish cleaning out the rest of the table

# ----------------------------------------------------------------------------------------------
# Now back to the original 
# Clean unneeded data
  # Remove 'Syn Mock' sample and add other metadata 
  clean_mt_16s$sample_table <- clean_mt_16s$sample_table |> 
    slice(
    str_which(mt_16s$sample_table$fastqFile, "Syn.*", negate = T)
      ) |> 
    rownames_to_column(var = "rownames") |> 
      # Add in the other metadata
    left_join(y = sw_meta) |> 
      # move rownames back
    column_to_rownames(var = "rownames")

  nrow(clean_mt_16s$sample_table)
# Remove anything not Archaea or Bacteria
   clean_mt_16s$tax_table <- clean_mt_16s$tax_table |> 
    dplyr::slice(
    stringr::str_which(clean_mt_16s$tax_table$k, ".*Bacteria|.*Archaea")
   )
# Check to be sure we removed them
    clean_mt_16s$tax_table$k |> 
      unique()
   # Remove the samples that are blanks
   clean_mt_16s$sample_table <- clean_mt_16s$sample_table |> 
     filter_out(sample %in% control_samps)
clean_mt_16s
clean_mt_16s$tidy_dataset()
clean_mt_16s

# --------------------------------------------------------
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

# Calculate test statistics and things
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
# view(syn_table)

# Now remove taxa where the mean abundance of reads across samples is larger than the reads
# in the synthetic community

syn_table <- syn_table |> 
  filter_out(mean_abund > `SynMock-Sheik_S1`)
# syn_table |> 
  # view()

# Now filter out taxa that had a significant t_test statistic (there is a significant difference
# between the real samples and the contamination samples)
syn_table <- syn_table |> 
  filter(t_test_diff > 0.05)
nrow(syn_table)

# view(syn_table)

# Okay save this list as possible contaminants
syn_contam <- syn_table |> 
  mutate(otu = paste0("o__", otu)) |> 
  select(otu) |> 
  as_vector()

head(syn_contam)
# ----------------------------------------------------

# Try micRoclean --------------------------------------------------------------------------------
library(micRoclean)



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
