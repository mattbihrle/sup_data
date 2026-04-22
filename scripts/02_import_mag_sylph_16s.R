
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
library(tidyverse) #tidyverse last


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
mt_sylph$tidy_dataset()
 # Calculations for later
# Calculate relative abunance
mt_sylph$cal_abund()
mt_sylph$taxa_abund$p[1:5,1:5]

# Calc Alpha Diversity
mt_sylph$cal_alphadiv()

# Calc Bray Curtis Dissimilarity
mt_sylph$cal_betadiv(method = c("bray","aitchison","robust.aitchison"))
mt_sylph$beta_diversity$bray

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
       select(order(colnames(mt_16s$otu_table)))
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
    colnames(mt_16s$tax_table) 
  # Replace question marks with 'unknown'
    mt_16s |> tidy_taxonomy(pattern = "\\?", replacement = "unknown")
    mt_16s$tax_table[1:5,]
    
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
# Remove anything not Archaea or Bactera
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

  # Tidy dataset
  print("Tidying dataset")
  mt_16s
  mt_16s$tidy_dataset()
  mt_16s

mt_16s$sample_table$strat_season <- factor(mt_16s$sample_table$strat_season, 
  levels = c("summer", "fall", "winter", "spring"), 
  labels = c("Summer", "Fall", "Winter", "Spring"),
   ordered = T)
## Calculate relative abunance-------------------------------------------------------------------------------------

mt_16s$cal_abund()
mt_16s$taxa_abund$p[1:5,1:5]

# Calc Alpha Diversity
mt_16s$cal_alphadiv()

# Calc Bray Curtis Dissimilarity
mt_16s$cal_betadiv(method = c("bray","aitchison","robust.aitchison"))
mt_16s$beta_diversity$bray

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

# Import coverage/abundance data
de_coverage <- read_tsv("data/derep_bins/Drep_Bins_coverage.tsv")
  
de_coverage$genome <- str_replace_all(de_coverage$Genome, pattern = "_MAGScoT_cleanbin_000", replacement = "_")

colnames(de_coverage) <- str_remove_all(colnames(de_coverage), "001_val_.{2}fq.gz ")

# Recalculate relative abundances
   # First extract only forward reads and relative abundance to keep
col_keep <- colnames(de_coverage[str_which(colnames(de_coverage), 
  ".*1_Relative.*")])
col_keep
# Remove columns we dont want and the first 'unmapped' row
de_coverage_clean <- de_coverage |> 
  select(genome, all_of(col_keep)) |> 
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

mt_mag$sample_table$strat_season <- factor(mt_mag$sample_table$strat_season, 
  levels = c("summer", "fall", "winter", "spring"), 
  labels = c("Summer", "Fall", "Winter", "Spring"),
   ordered = T)
mt_mag$tidy_dataset()

 # Calculations for later
# Calculate relative abunance
mt_mag$cal_abund()
mt_mag$taxa_abund

# Calc Alpha Diversity
mt_mag$cal_alphadiv()

# Calc Bray Curtis Dissimilarity
mt_mag$cal_betadiv(method = c("bray","aitchison","robust.aitchison"), unifrac = T)
head(mt_mag$beta_diversity)

# Save the tables we care about
save(mt_mag, file = "output/data/mt_mag.RData")
save(mt_sylph, file = "output/data/mt_sylph.RData")
save(mt_16s, file = "output/data/mt_16s.RData")
# # Remove unneeded data----------------------------------------------------
remove <- ls() |> 
  str_subset("mt_.*", negate = T)

rm(list = remove)

rm(remove)
