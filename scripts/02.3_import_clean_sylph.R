# Import and clean sylph data

# Load Libraries
library(microeco)
library(stringdist)
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

# view(manual_check)

# Listing out genomes that do match genus and family in contam_tax_16s to label as poss contaminants
extra_contam_taxa <- 
  c("t__GCF_000702605.1", #all from exiguobacteraceae 
    "t__GCF_001939065.1", 
    "t__GCF_014524545.1", 
    "t__GCA_902363455.1",
    "t__GCA_945860265.1", # from spirosomaceae family
    "t__GCA_003259005.1",
    "t__GCA_030653545.1", # from nitrospira genus
    "t__GCA_021300015.1", # from planktomycetes, pirellula genus
    "t__GCA_021736575.1",
    "t__GCA_945900945.1", 
    "t__GCA_020718185.1", # removed from paracaedibacterium
    "t__GCF_014138435.1", # removed from rhizobiales
    "t__GCA_903842455.1",
    "t__GCA_002336985.1", # remove from sphingomonas
    "t__GCF_009768975.1", 
    "t__GCA_013823985.1",
    "t__GCF_029865925.1",
    "t__GCF_945890115.1", 
    "t__GCA_945903445.1", # removed from burkolderiaceae
    "t__GCA_945878355.1", 
    "t__GCA_009693505.1",
    "t__GCA_937897495.1",
    "t__GCA_945901045.1",
    "t__GCA_947447245.1",
    "t__GCF_900182955.1",
    "t__GCF_006974105.1",
    "t__GCA_004293725.1",
    "t__GCA_016463455.1",
    "t__GCF_006364715.1",
    "t__GCA_903914945.1",
    "t__GCA_947499625.1",
    "t__GCF_900103875.1",     # removed from psuedomonas
    "t__GCF_902498065.1",
    "t__GCF_900105495.1",
    "t__GCF_025397885.1",
    "t__GCF_900106065.1",
    "t__GCF_015461845.1",
    "t__GCF_021602155.1",
    "t__GCF_004683905.1",
    "t__GCF_007858255.1",
    "t__GCF_019145195.1",
    "t__GCF_002836515.1",
    "t__GCF_010095445.2",
    "t__GCF_900187425.1")

# Set those genomes as TRUE in the possible contaminant column 
mt_sylph$tax_table <- mt_sylph$tax_table |> 
  mutate(poss_16s_contaminant = ifelse(t %in% extra_contam_taxa, TRUE, poss_16s_contaminant))

table(mt_sylph$tax_table$poss_16s_contaminant)
mt_sylph$tax_table <- mt_sylph$tax_table |> 
  filter(poss_16s_contaminant == FALSE)

# Okay now remove the stations that are contamination (WM17-WM24)

mt_sylph$sample_table <- mt_sylph$sample_table |> 
  filter_out(str_detect(sample, "WM(1[7-9]{1}|2[0-9]{1})"))
#Tidy dataset 
mt_sylph
mt_sylph$tidy_dataset()
mt_sylph


# Calculations for later -----------------------------------------------------------
# Calculate relative abunance
mt_sylph$cal_abund(select_cols = clade_cols)
head(mt_sylph$taxa_abund)

# Calc Alpha Diversity
mt_sylph$cal_alphadiv()

# Calc Bray Curtis Dissimilarity
mt_sylph$cal_betadiv(method = c("bray","aitchison","robust.aitchison", "jaccard"))
head(mt_sylph$beta_diversity)

save(mt_sylph, file = "output/data/mt_sylph.RData")

remove <- ls() |> 
  str_subset("mt_.*", negate = T)

rm(list = remove)
rm(remove)
