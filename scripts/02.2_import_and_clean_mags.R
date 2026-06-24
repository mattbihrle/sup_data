# Load Libraries
# install.packages("stringdist")
library(microeco)
library(stringdist)
library(tidyverse)


# Import sample table -------------------------------------------------------------
sw_meta <- read_csv("output/data/metadata_supwinter.csv") |> 
  dplyr::select(!matches(c("12h_med", "round_date")))|> 
  distinct()
# Create a sample table for import into the microtable
sample_table <- sw_meta |> 
    column_to_rownames("sample") |> 
    as.data.frame()
# Add sample back as a column
sample_table <- sample_table |> 
  mutate(sample = rownames(sample_table))

# Create all the factors that I need to 

sample_table <- sample_table |> 
  mutate(strat_season = fct_reorder(strat_season, date), 
         strat_season_2 = fct_reorder(strat_season_2, date),
         solar_season = fct_reorder(solar_season, date), 
         mixing = fct_relevel(mixing, "stratified_std", "mixed", "stratified_inverse"))

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

# Pull out all the bins that are WM17 and larger
bad_bins <- mag_df_clean |> 
  filter(grepl("WM(2.{1}|1[7-9]{1}).*", genome))

# Remove those bins from the master list
mag_df_clean <- mag_df_clean |> 
  filter_out(genome %in% bad_bins$genome)
# Add a possible_contaminant column to the df
  # first load in the 16s contaminants
load("output/data/contam_16s.RData")
contam_tax_16s <- contam_16s$tax_table
mag_df_clean <- mag_df_clean |> 
  # Match create a column based on exact matches and on fuzzy matches to look through by hand
  mutate(poss_16s_contaminant = ifelse(mag_df_clean$f %in% contam_tax_16s$f & mag_df_clean$g %in% contam_tax_16s$g, TRUE, FALSE),
        poss_fuzzy_contaminant = ifelse(ain(mag_df_clean$f,  contam_tax_16s$f, method = "jaccard", maxDist = 0.2) & ain(mag_df_clean$g, contam_tax_16s$g,method = "jaccard", maxDist = 0.2), TRUE, FALSE))
table(mag_df_clean$poss_16s_contaminant)
table(mag_df_clean$poss_fuzzy_contaminant)
# Filter out where the exact match and the fuzzy match differ
manual_check <- mag_df_clean |> 
  filter(poss_16s_contaminant != poss_fuzzy_contaminant)

# view(manual_check)

# Listing out genomes that do match genus and family in contam_tax_16s to label as poss contaminants
extra_contam_genomes <- 
  c("WM03_S3_018", 
    "WM01_S1_022", 
    "WM01_S1_030", 
    "WM02_S2_019", 
    "WM04_S4_030", 
    "WM05_S5_022",
    "WM12_S7_045",
    "WM15_S10_032", 
    "WM04_S4_059",
    "WM14_S9_009",
    "WM15_S10_006",
    "WM16_S11_015", 
    "WM03_S3_045",
    "WM01_S1_004",
    "WM05_S5_024",
    "WM16_S11_016",
    "WM01_S1_015",
    "WM03_S3_043",
    "WM16_S11_035")

# Set those genomes as TRUE in the possible contaminant column 
mag_df_clean <- mag_df_clean |> 
  mutate(poss_16s_contaminant = ifelse(genome %in% extra_contam_genomes, TRUE, poss_16s_contaminant))

mag_df_clean_small <- mag_df_clean |> 
  filter(poss_16s_contaminant == FALSE)


# Okay from here go through the sheik table once more to see what phyla and genuses could still be contamination
sheik_contam <- read_tsv("sheik_2018_contamination_table.txt") |> 
  mutate(p = paste0("p__", p), 
         c = paste0("c__", c), 
         g = paste0("g__", g))

# From that look, I think we have everything from the Sheik list removed. 

# Extract the final MAG list and save it
mag_df_clean_small <- mag_df_clean_small |> 
  dplyr::select("genome") |> 
  as_vector() |> 
  writeLines("output/data/clean_bins.txt")


# Import the dereplicated bins list
de_rep_bins <- read_tsv("msi_downloads/Derep_list_clean_6_19.tsv", col_names = "bin") |> 
  mutate(bin = str_extract(bin, "WM.*_[0-9]{6}")) |> 
  mutate(bin = str_remove(bin, "MAG.*cleanbin_000"))

# Import coverage/abundance data from coverM
de_coverage <- read_tsv("msi_downloads/Drep_Bins_coverage_clean_6_19.tsv")
  
de_coverage$genome <- str_replace_all(de_coverage$Genome, pattern = "_MAGScoT_cleanbin_000", replacement = "_")

colnames(de_coverage) <- str_remove_all(colnames(de_coverage), "001_val_.{2}fq.gz ")

# Recalculate relative abundances
   # First extract mean and relative abundance to keep
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

# Quick check to see how column abundances change from WM17-WM24
colSums(de_coverage_clean[2:ncol(de_coverage_clean)])
# Check to see which bins are shared between the coverage and the dereplicated list

test <- mag_df_clean |> 
  mutate(derep_bin = ifelse(genome %in% de_rep_bins$bin, TRUE, FALSE), 
        coverage_bin = ifelse(genome %in% de_coverage_clean$genome, TRUE, FALSE)
) |> 
  select(genome, derep_bin, coverage_bin, everything()) |> 
  filter_out(derep_bin == FALSE & coverage_bin == FALSE) |> 
  arrange(desc(derep_bin), desc(coverage_bin))

sum(test$derep_bin == TRUE & test$coverage_bin == TRUE)
#Hmmm looks like only 45 are included in both the abundance calculations and the dereplicated list. This is odd but I will just move forward
# with using the abundance list. 

# import tree data
tree_files <- list.files("data/gtdb_bin_tax/classify/", pattern = ".tree$", full.names = T) 
tree <- ape::read.tree(tree_files[9])
# Update tree tip labels
tree$tip.label <- tree$tip.label |> 
  str_remove("MAGScoT_cleanbin_000")

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

# Remove WM17-WM24

mt_mag$sample_table <- mt_mag$sample_table |> 
  filter_out(sample %in% c("WM17", "WM18", "WM19", "WM20", "WM21", "WM22", "WM23", "WM24"))


mt_mag
mt_mag$tidy_dataset()
mt_mag

# That leaves us with 116 MAGs 

# Next import the quality data ----------------------------------------------------------------------

all_bins_qual <- read_tsv("msi_downloads/all_bins_quality_report.tsv") |> 
  mutate(quality = ifelse(Completeness > 90 & Contamination < 5, "HIGH", "na")) |> 
  mutate(quality = ifelse(Completeness >= 50 & Contamination < 10 & quality != "HIGH", "MED", quality)) |>
  mutate(quality = ifelse(Completeness < 50 & Contamination < 10, "LOW", quality)) |> 
  mutate(quality = ifelse(Contamination >= 10, "contamination_flag", quality)) |> 
  select(Name, Completeness, Contamination, quality, everything()) |> 
  #rename the names to simple bins 
  mutate(bin = paste0("b__", str_remove(Name, "MA.*_000")))

mt_mag$tax_table <- mt_mag$tax_table |> 
  left_join(all_bins_qual, "bin")


# Calculations for later -----------------------------------------------------------
# Calculate relative abunance
mt_mag$cal_abund(select_cols = clade_cols)
head(mt_mag$taxa_abund)

# Calc Alpha Diversity
mt_mag$cal_alphadiv()

# Calc Bray Curtis Dissimilarity
mt_mag$cal_betadiv(method = c("bray","aitchison","robust.aitchison", "jaccard"), unifrac = T)
head(mt_mag$beta_diversity)

save(mt_mag, file = "output/data/mt_mag.RData")

remove <- ls() |> 
  str_subset("mt_.*", negate = T)

rm(list = remove)
rm(remove)
