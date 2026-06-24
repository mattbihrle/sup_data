# Import and clean sylph data

# Load Libraries
library(microeco)
library(tidyverse)

# # Import sample table -------------------------------------------------------------
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
    as.data.frame() |> 
  tidy_taxonomy()

# Create mt object
mt_sylph <-microtable$new(otu_table = otu_table, sample_table = sample_table, tax_table = tax_table)


mt_sylph
mt_sylph$tidy_dataset()
mt_sylph

# Filter out contamination -----------------------------------------------------------------------
load("output/data/contam_16s.RData")
contam_tax_16s <- contam_16s$tax_table
mt_sylph$tax_table <- mt_sylph$tax_table |> 
  # Match create a column based on exact matches and on fuzzy matches to look through by hand
  mutate(poss_16s_contaminant = ifelse(mt_sylph$tax_table$f %in% contam_tax_16s$f & mt_sylph$tax_table$g %in% contam_tax_16s$g, TRUE, FALSE),
        poss_fuzzy_contaminant = ifelse(ain(mt_sylph$tax_table$f,  contam_tax_16s$f, method = "jaccard", maxDist = 0.2) & ain(mt_sylph$tax_table$g, contam_tax_16s$g,method = "jaccard", maxDist = 0.2), TRUE, FALSE))
table(mt_sylph$tax_table$poss_16s_contaminant)
table(mt_sylph$tax_table$poss_fuzzy_contaminant)
# Filter out where the exact match and the fuzzy match differ
manual_check <- mt_sylph$tax_table |> 
  filter(poss_16s_contaminant != poss_fuzzy_contaminant)

view(manual_check)

# Listing out genomes that do match genus and family in contam_tax_16s to label as poss contaminants
extra_contam_genomes <- 
  c()
#START HERE NEXT TIME
# Set those genomes as TRUE in the possible contaminant column 
mag_df_clean <- mag_df_clean |> 
  mutate(poss_16s_contaminant = ifelse(genome %in% extra_contam_genomes, TRUE, poss_16s_contaminant))

mag_df_clean_small <- mag_df_clean |> 
  filter(poss_16s_contaminant == FALSE)
