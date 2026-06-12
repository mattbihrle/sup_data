# install extra packages

# remotes::install_github("rachelgriffard/micRoclean", build_vignettes = TRUE, force = TRUE)
# pak::pak("rachelgriffard/micRoclean")

# Load Packages ------------------------------
  # For analysis
library(vegan) 
library(microeco)
# For importing .biom file
library(phyloseq)
library(file2meco)
library(ape) # for reading in tree data
library(GUniFrac) # for unifrac diversity
# Extra packages for plotting
# library(ggalluvial)
# library(ggnested) # Needs to be installed from github not CRAN
# library(ggarrow)
# library(Polychrome) 
# library(forcats)
# library(paletteer)
# library(plotly)
library(micRoclean)
library(tidyverse) #tidyverse last


# Import meta data ----------------------------------------------------------
sw_meta <- read_csv("output/data/metadata_supwinter.csv") |> 
  dplyr::select(!matches(c("12h_med", "round_date")))|> 
  distinct()

# Import list of contaminants from Sheik et al 2018 Supplemental Material

contam_genus <- read_tsv("sheik_2018_contamination_table.txt") |> 
  rename_with(tolower) |>  
  filter_out(extraction == FALSE & water == FALSE) |> 
  dplyr::select(g) |> 
  as_vector()
# Load in and clean sylph data-----------------------------------------------------------------------------------

b_df <- readr::read_tsv("data/superior_sylph/Sylph_TaxAbund_out.tsv")

    # Clean

# Create vectors of column names
clade_cols <- stringr::str_extract_all(b_df$clade_name[8], "[a-z]{1}(?=_)")

clade_cols <- as.vector(clade_cols[[1]]) %>% print()

site_cols <- 
  stringr::str_extract(colnames(b_df), "(?<=/).{4}") %>%   
  na.omit() %>% 
  print() 

colnames(b_df)[2:ncol(b_df)] <- site_cols

#Spread columns out, clean names

b_df <- b_df %>% tidyr::separate_wider_delim(cols = clade_name, delim = "|", names = clade_cols, too_few = "align_start") %>% 
  # Drop all rows that don't get down to a taxa 
  tidyr::drop_na("t") %>%
  #remove additional characters and make data columns numeric
  mutate(across(clade_cols, ~ str_remove(., "[a-z]__")))

# Quick check to see if abundances are all still 100 
colSums(b_df[9:ncol(b_df)])

# Create df of all taxa found in the samples that weren't pumped to act as potential contamination
con_df <- b_df |> 
  dplyr::select(!WM01:WM16) |> 
  mutate(across(any_of(site_cols), ~ifelse(.x == 0.00, NA, .x))) |> 
  drop_na()

# Double check that there are no dupliated species
nrow(con_df) == length(unique(con_df$t))

# Great, now put all those taxa in a vector
con_vec <- unique(con_df$t)

#Create long df 
long_b_df <- b_df %>% pivot_longer(cols = all_of(site_cols), 
  names_to = "sample", values_to = "abundance")

long_b_df <- long_b_df %>% 
  mutate(date = sw_meta$date[match(long_b_df$sample, sw_meta$sample)]
  )
long_b_df$date <- ymd(long_b_df$date)

# Turn the Sylph data into a microeco object ('mt_sylph')----------------------------------------------------------------------------------

otu_table <- b_df |> 
    rownames_to_column("otu") |> 
    dplyr::select(all_of(site_cols), "otu") |> 
    column_to_rownames("otu") |> 
    as.data.frame()

tax_table <- b_df |> 
    rownames_to_column("otu") |> 
    dplyr::select(all_of(clade_cols), "otu") |> 
    column_to_rownames("otu") |> 
    as.data.frame()

sample_table <- sw_meta |> 
    column_to_rownames("sample") |> 
    as.data.frame()
# Add sample back as a column
sample_table <- sample_table |> 
  mutate(sample = rownames(sample_table))
# Create mt object
mt_sylph <-microtable$new(otu_table = otu_table, sample_table = sample_table, tax_table = tax_table)
# Reorder factors of stratification
mt_sylph$sample_table$strat_season <- factor(mt_sylph$sample_table$strat_season, 
  levels = c("summer", "fall", "winter", "spring"), 
  labels = c("Summer", "Fall", "Winter", "Spring"),
   ordered = T)

# # Remove contamination ---------------------------------------------------------------
# # Create vector of sites that were not pumped
#   con_site_cols <- site_cols[12:19]

# # Create df of all taxa found in the samples that weren't pumped to act as potential contamination
# con_df <- b_df |> 
#   dplyr::select(!WM01:WM16) |> 
#   mutate(across(any_of(site_cols), ~ifelse(.x == 0.00, NA, .x))) |> 
#   drop_na()

# # Double check that there are no dupliated species
# nrow(con_df) == length(unique(con_df$t))

# # Great, now put all those taxa in a vector
# con_vec <- unique(con_df$t)

# # Remove taxa that are contamination
# mt_sylph$filter_pollution(con_vec)

# # Now remove the samples that did not filter
# mt_sylph$sample_table <- mt_sylph$sample_table |> 
#   filter_out(sample %in% con_site_cols)
# # -------------------------------------------------------------------
# mt_sylph$tidy_dataset()
#  # Calculations for later
# # Calculate relative abunance
# mt_sylph$cal_abund()
# mt_sylph$taxa_abund$p[1:5,1:5]

# # Calc Alpha Diversity
# mt_sylph$cal_alphadiv()

# # Calc Bray Curtis Dissimilarity
# mt_sylph$cal_betadiv(method = c("bray","aitchison","robust.aitchison", "jaccard"))
# mt_sylph$beta_diversity$bray

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
# Okay here start working with the micRoclean package for contamination------------------------------------------------------------
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
  dplyr::select(sample_type, is_control, sample_well, batch) |>
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

# remove the taxa that are all 0s 
  mclean_results$decontaminated_count <- 
    mclean_results$decontaminated_count |> 
    as_data_frame() |> 
    setNames(colnames(count)) |> 
    dplyr::select(where(~ !all(.x == 0)))
# Here I want to set it up so it will run through the control samps one by one
# so there is not so much contamination taken out. 

results_list <- list()

for(i in 1:8){
  counts_new <- count |> 
    t() |> 
    as_data_frame() |> 
    dplyr::select(control_names[1], matches(control_samps[i]), !matches(control_samps)) |> 
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
    dplyr::select(where(~ !all(.x == 0)))

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
    all_blank = ifelse(otu %in% paste0("o__", colnames(mclean_results$decontaminated_count)), TRUE, FALSE)
  )
mt_16s$tax_table <- mt_16s$tax_table |> 
  rownames_to_column("rows") |> 
  rowwise() |> 
  dplyr::mutate(true_sum = sum(c_across(matches("blank")))) |> 
  ungroup() |> 
  column_to_rownames("rows")

# Okay now filter out all the taxa that are FALSE for all of the blank samples
clean_mt_16s <- clone(mt_16s)

clean_mt_16s$tax_table <- clean_mt_16s$tax_table |> 
  filter(if_any(matches("blank")) == TRUE)

# Import list of contaminants from Sheik et al 2018 Supplemental Material

contam_genus <- read_tsv("sheik_2018_contamination_table.txt") |> 
  rename_with(tolower) |>  
  filter_out(extraction == FALSE & water == FALSE) |> 
  dplyr::select(g) |> 
  as_vector()

clean_mt_16s
clean_mt_16s$tax_table <- clean_mt_16s$tax_table |> 
  filter_out(g %in% paste0("g__", contam_genus))
clean_mt_16s

# Okay, now back to the original cleaning and things

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


clean_mt_16s$sample_table$strat_season <- factor(clean_mt_16s$sample_table$strat_season, 
  levels = c("summer", "fall", "winter", "spring"), 
  labels = c("Summer", "Fall", "Winter", "Spring"),
   ordered = T)

   # Remove the samples that are blanks
   clean_mt_16s$sample_table <- clean_mt_16s$sample_table |> 
     filter_out(sample %in% control_samps)

# Remove the samples that are
clean_mt_16s
clean_mt_16s$tidy_dataset()
clean_mt_16s
## Calculate relative abunance-------------------------------------------------------------------------------------

clean_mt_16s$cal_abund()
clean_mt_16s$taxa_abund$p[1:5,1:5]

# Calc Alpha Diversity
clean_mt_16s$cal_alphadiv()

# Calc Bray Curtis Dissimilarity
clean_mt_16s$cal_betadiv(method = c("bray","aitchison","robust.aitchison", "jaccard"))
clean_mt_16s$beta_diversity$bray


# Import MAG data -----------------------------------------------------
bac_df <- readr::read_tsv("data/gtdb_bin_tax/gtdbtk.bac120.summary.tsv", na = "N/A")
arc_df <- read_tsv("data/gtdb_bin_tax/gtdbtk.ar53.summary.tsv", na = "N/A")

mag_df <- bind_rows(bac_df, arc_df) |> 
  mutate(sample = str_extract(user_genome, "WM.{2}"), .after = user_genome) |> 
    mutate(genome = str_remove(user_genome, "MAGScoT_cleanbin_000")) |>
    dplyr::select(genome, classification)
  
  # Create vectors of column names
  clade_cols <- stringr::str_extract_all(mag_df$classification[3], "[a-z]{1}(?=_)")
  
  clade_cols <- as.vector(clade_cols[[1]]) %>% print()

#Spread columns out, clean names

mag_df_clean <- mag_df %>% tidyr::separate_wider_delim(
  cols = classification, delim = ";", names = clade_cols, too_few = "align_start")

# Find possible contamination from the 16s data
contam_tax_16s <- mt_16s$tax_table |> 
  filter_out(mt_16s$tax_table$otu %in% clean_mt_16s$tax_table$otu)

# Add a possible_contaminant column to the df

mag_df_clean <- mag_df_clean |> 
  mutate(poss_contaminent = ifelse(mag_df_clean$f %in% contam_tax_16s$f & mag_df_clean$g %in% contam_tax_16s$g, TRUE, FALSE))
# Import coverage/abundance data
de_coverage <- read_tsv("msi_downloads/Drep_Bins_coverage_04_23.tsv")
  
de_coverage$genome <- str_replace_all(de_coverage$Genome, pattern = "_MAGScoT_cleanbin_000", replacement = "_")

colnames(de_coverage) <- str_remove_all(colnames(de_coverage), "001_val_.{2}fq.gz ")

# Recalculate relative abundances
   # First extract only forward reads and relative abundance to keep
col_keep <- colnames(de_coverage[str_which(colnames(de_coverage), 
  ".*1_Relative.*")])
col_keep
# Remove columns we dont want and the first 'unmapped' row
de_coverage_clean <- de_coverage |> 
  dplyr::select(genome, all_of(col_keep)) |> 
  slice(-1)

# Create a vector of sites
site_cols <- 
  stringr::str_extract(col_keep, "^.{4}") %>%   
  na.omit() %>% 
  unique() |> 
  print() 

# Rename long names to match'site_cols'
de_coverage_clean <- de_coverage_clean |> 
  rename_with(
    ~site_cols, 
    .cols = matches("WM.*")
  )
# # Calculate relative abundances  
# de_coverage_clean <- de_coverage_clean |> 
#   mutate(across(str_which(colnames(de_coverage_clean), "WM.*"),
#    ~ .x/sum(.x, na.rm = T)*100))

# Quick check to see if abundances are all still 100 
colSums(de_coverage_clean[2:ncol(de_coverage_clean)])

# import tree data
tree_files <- list.files("data/gtdb_bin_tax/classify/", pattern = ".tree$", full.names = T) 
tree <- ape::read.tree(tree_files[9])

# Update tree tip labels
tree$tip.label <- tree$tip.label |> 
  str_remove("MAGScoT_cleanbin_000")
# plot(tree[[1]])
# Turn the mag data into a microeco object ('mt_mag')----------------------------------------------------------------------------------

  # Use the de_rep coverage data for otus and abundance
otu_table <- de_coverage_clean |> 
    column_to_rownames("genome") |> 
  as.data.frame()
head(otu_table)

tax_table <- mag_df_clean |> 
  mutate(bin = genome) |> 
    column_to_rownames("genome") |> 
    as.data.frame() |> 
  tidy_taxonomy()
head(tax_table)

# Use sample table from above here

head(sample_table)
mt_mag <-microtable$new(otu_table = otu_table, sample_table = sample_table, 
                        tax_table = tax_table, phylo_tree = tree)
# mt_mag <- trans_norm$new(mt_mag)
# mt_mag <- mt_mag$norm(method = "rclr")
mt_mag_clone <- clone(mt_mag)
mt_mag_clone$tidy_dataset()
## Remove contamination ---------------------------------------------------------------
# Create vector of sites that were not pumped
  con_site_cols <- site_cols[12:19]

# Create df of all taxa found in the samples that weren't pumped to act as potential contamination
con_df <- mt_mag$otu_table |> 
  dplyr::select(!matches("WM01"): matches("WM16")) |> 
  mutate(across(everything(), ~ifelse(.x == 0.00, NA, .x))) |> 
  drop_na() |> 
  rownames_to_column("bin") |> 
  mutate(id = bin) |> 
  column_to_rownames("id")

# Double check that there are no dupliated species/otu/bin
nrow(con_df) == length(unique(con_df$bin))

# Great, now put all those taxa in a vector
con_vec <- unique(con_df$bin)

# Remove taxa that are contamination
mt_mag$filter_pollution(con_vec)
mt_mag$filter_pollution(con_site_cols)

# also remove any bins that are WM17-WM24
mt_mag$tax_table <- mt_mag$tax_table |> 
  filter_out(bin %in% con_site_cols)

# Now remove the bins from samples that didn't filter (eg WMxx_S17_021)
  # First create a vector of those names
contam_bins <- rownames(mt_mag$otu_table) |> 
  str_subset("WM(2[0-9]{1}|1[7-9]{1}).*") 
  # then filter_pollution to remove

mt_mag$filter_pollution(contam_bins)

#-----------------------------------------------------

mt_mag$sample_table$strat_season <- factor(mt_mag$sample_table$strat_season, 
  levels = c("summer", "fall", "winter", "spring"), 
  labels = c("Summer", "Fall", "Winter", "Spring"),
   ordered = T)
mt_mag
mt_mag$tidy_dataset()
mt_mag
 # Calculations for later
# Calculate relative abunance
mt_mag$cal_abund()
mt_mag$taxa_abund

# Calc Alpha Diversity
mt_mag$cal_alphadiv()

# Calc Bray Curtis Dissimilarity
mt_mag$cal_betadiv(method = c("bray","aitchison","robust.aitchison", "jaccard"), unifrac = T)
head(mt_mag$beta_diversity)

# Save the tables we care about
mt_16s <- clean_mt_16s
save(mt_mag, file = "output/data/mt_mag.RData")
save(mt_sylph, file = "output/data/mt_sylph.RData")
save(mt_16s, file = "output/data/mt_16s.RData")
# # Remove unneeded data----------------------------------------------------
remove <- ls() |> 
  str_subset("mt_.*", negate = T)

rm(list = remove)

rm(remove)
