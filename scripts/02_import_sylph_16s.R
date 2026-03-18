# Introductory Script

# Load Packages ------------------------------
  # For analysis
library(vegan) 
library(microeco)
library(tidyverse)
# For importing .biom file
library(phyloseq)
library(file2meco)
# Extra packages for plotting
library(ggalluvial)
library(ggnested) # Needs to be installed from github not CRAN
library(ggarrow)
library(Polychrome) 
library(forcats)
library(paletteer)
library(plotly)


# Import meta data ----------------------------------------------------------
sw_meta <- read_csv("output/data/metadata_supwinter.csv")
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



#Create long df 
long_b_df <- b_df %>% pivot_longer(cols = all_of(site_cols), 
  names_to = "sample", values_to = "abundance")

long_b_df <- long_b_df %>% 
  mutate(date = sw_meta$date[match(long_b_df$sample, sw_meta$sample)]
  )
long_b_df$date <- mdy(long_b_df$date)

# Turn the Sylph data into a microeco object ('mt_sylph')----------------------------------------------------------------------------------

otu_table <- b_df |> 
    rownames_to_column("otu") |> 
    select(all_of(site_cols), "otu") |> 
    column_to_rownames("otu") |> 
    as.data.frame()

tax_table <- b_df |> 
    rownames_to_column("otu") |> 
    select(all_of(clade_cols), "otu") |> 
    column_to_rownames("otu") |> 
    as.data.frame()

sample_table <- sw_meta |> 
    column_to_rownames("sample") |> 
    as.data.frame()

mt_sylph <-microtable$new(otu_table = otu_table, sample_table = sample_table, tax_table = tax_table)

mt_sylph$tidy_dataset()
 # Calculations for later
# Calculate relative abunance
mt_sylph$cal_abund()
mt_sylph$taxa_abund$p[1:5,1:5]

# Calc Alpha Diversity
mt_sylph$cal_alphadiv()

# Calc Bray Curtis Dissimilarity
mt_sylph$cal_betadiv()
mt_sylph$beta_diversity$bray

# Import 16S data----------------------------------------------------------------------------
  # First the sample names
  names_df <- readr::read_tsv("data/16s/lotus3_out/final_sample_map_4_R.txt")
  names_vec <- setNames(names_df$SampleID, names_df$RealNames)


  # Load microeco object
    meco_16s <- phyloseq::import_biom("data/16s/lotus3_out/OTU.biom") |> 
        phyloseq2meco()
      # Rename columns
      meco_16s$otu_table <- meco_16s$otu_table |> 
        rename(any_of(names_vec))
      # Reorder columns
      meco_16s$otu_table <- meco_16s$otu_table |> 
       select(order(colnames(meco_16s$otu_table)))
      colnames(meco_16s$otu_table)
      # rename rows
    meco_16s$sample_table <- meco_16s$sample_table |> 
      rownames_to_column("sample") |> 
      dplyr::left_join(names_df, by = join_by(sample == SampleID)) |> 
      column_to_rownames("RealNames") 
    # Create a column of sample names for later metadata
    meco_16s$sample_table <- meco_16s$sample_table|> 
      mutate(sample = str_extract(rownames(meco_16s$sample_table), ".{4}"))
      # Reorder rows alphabetically
      meco_16s$sample_table <- meco_16s$sample_table |> 
        slice(order(rownames(meco_16s$sample_table)))
    print(rownames(meco_16s$sample_table))
    # Rename Taxa
    taxa_names <- setNames( c("k", "p", "c", "o", "f", "g", "s"), colnames(meco_16s$tax_table))
    colnames(meco_16s$tax_table) <- taxa_names
    colnames(meco_16s$tax_table) 
  # Replace question marks with 'unknown'
    meco_16s |> tidy_taxonomy(pattern = "\\?", replacement = "unknown")
    meco_16s$tax_table[1:5,]
    
# Clean unneeded data
  # Remove 'Syn Mock' sample 
  meco_16s$sample_table <- meco_16s$sample_table |> 
    slice(
    str_which(meco_16s$sample_table$fastqFile, "Syn.*", negate = T)
      ) |> 
    rownames_to_column(var = "rownames") |> 
      # Add in the other metadata
    left_join(y = sw_meta) |> 
      # move rownames back
    column_to_rownames(var = "rownames")

  nrow(meco_16s$sample_table)
# Remove anything not Archaea or Bactera
   meco_16s$tax_table <- meco_16s$tax_table |> 
    dplyr::slice(
    stringr::str_which(meco_16s$tax_table$k, ".*Bacteria|.*Archaea")
   )
# Check to be sure we removed them
    meco_16s$tax_table$k |> 
      unique()
  # Remove mitochondria and chloroplasts
  print("removing mitochontria and chloroplast")
  meco_16s
  meco_16s$filter_pollution(taxa = c("mitochondria", "chloroplast", "metagenome"))
  meco_16s      

  # Tidy dataset
  print("Tidying dataset")
  meco_16s
  meco_16s$tidy_dataset()
  meco_16s


## Calculate relative abunance-------------------------------------------------------------------------------------

meco_16s$cal_abund()
meco_16s$taxa_abund$p[1:5,1:5]

# Calc Alpha Diversity
meco_16s$cal_alphadiv()

# Calc Bray Curtis Dissimilarity
meco_16s$cal_betadiv()
meco_16s$beta_diversity$bray
# 
# Remove unneeded data----------------------------------------------------
rm(names_vec)
rm(site_cols)
rm(taxa_names)
rm(names_df)
rm(otu_table)
rm(sample_table)
rm(tax_table)
rm(clade_cols)
