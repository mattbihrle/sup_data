library(microeco)
library(file2meco)
library(tidyverse)
# Load in my little 16s data

load("output/data/mt_16s.RData")
# Import 16S data----------------------------------------------------------------------------
  # First the sample names
  names_df <- readr::read_tsv("data/16s/lotus3_out/final_sample_map_4_R.txt")
  names_vec <- setNames(names_df$SampleID, names_df$RealNames)


  # Load microeco object
    contam_mt_16s <- phyloseq::import_biom("data/16s/lotus3_out/OTU.biom") |> 
        phyloseq2meco()
      # Rename columns
      contam_mt_16s$otu_table <- contam_mt_16s$otu_table |> 
        rename(any_of(names_vec))
      # Reorder columns
      contam_mt_16s$otu_table <- contam_mt_16s$otu_table |> 
       dplyr::select(order(colnames(contam_mt_16s$otu_table)))
      colnames(contam_mt_16s$otu_table)
      # rename rows
    contam_mt_16s$sample_table <- contam_mt_16s$sample_table |> 
      rownames_to_column("sample") |> 
      dplyr::left_join(names_df, by = join_by(sample == SampleID)) |> 
      column_to_rownames("RealNames") 
    # Create a column of sample names for later metadata
    contam_mt_16s$sample_table <- contam_mt_16s$sample_table|> 
      mutate(sample = str_extract(rownames(contam_mt_16s$sample_table), ".{4}"))
      # Reorder rows alphabetically
      contam_mt_16s$sample_table <- contam_mt_16s$sample_table |> 
        slice(order(rownames(contam_mt_16s$sample_table)))
    print(rownames(contam_mt_16s$sample_table))
    # Rename Taxa
    taxa_names <- setNames( c("k", "p", "c", "o", "f", "g", "s"), colnames(contam_mt_16s$tax_table))
    colnames(contam_mt_16s$tax_table) <- taxa_names
    # Add an "otu" column for later use
    contam_mt_16s$tax_table <- contam_mt_16s$tax_table |> 
      mutate(otu = rownames(contam_mt_16s$tax_table))
    # check to make sure it worked
    rownames(contam_mt_16s$tax_table) == contam_mt_16s$tax_table$otu
    colnames(contam_mt_16s$tax_table) 
  # Replace question marks with 'unknown'
    contam_mt_16s |> tidy_taxonomy(pattern = "\\?", replacement = "unknown")
    contam_mt_16s$tax_table[1:5,]

  # Remove mitochondria and chloroplasts
  print("removing mitochontria and chloroplast")
  contam_mt_16s
  contam_mt_16s$filter_pollution(taxa = c("mitochondria", "chloroplast", "metagenome"))
  contam_mt_16s$tidy_dataset()
  contam_mt_16s     

  # Remove anything not Archaea or Bacteria
    contam_mt_16s$tax_table <-   contam_mt_16s$tax_table |> 
  dplyr::slice(
    stringr::str_which(  contam_mt_16s$tax_table$k, ".*Bacteria|.*Archaea")
  )
# Check to be sure we removed them
    contam_mt_16s$tax_table$k |> 
    unique()
  contam_mt_16s
  contam_mt_16s$tidy_dataset()
  contam_mt_16s
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
contam_tax_16s <-   contam_mt_16s$tax_table |> 
  filter_out(  contam_mt_16s$tax_table$otu %in% mt_16s$tax_table$otu)

# Add a possible_contaminant column to the df

mag_df_clean <- mag_df_clean |> 
  mutate(poss_contaminent = ifelse(mag_df_clean$f %in% contam_tax_16s$f & mag_df_clean$g %in% contam_tax_16s$g, TRUE, FALSE))
# Import coverage/abundance data
de_coverage <- read_tsv("msi_downloads/Drep_Bins_coverage_clean_6_11.tsv")
  
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

# Quick check to see what abundances look like once the unmapped row is removed. 
colSums(de_coverage_clean[2:ncol(de_coverage_clean)])

# import tree data
#  # Import and concatenate all tree files together - maybe work on this later, doesn't seem that important right now
# all_trees <- do.call(c, lapply(tree_files, ape::read.tree))
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

# Take this out before running in full script ---------------------------------

rownames<- paste0(mt_16s_contam$sample_table$sample)
sample_table <- mt_16s_contam$sample_table |> 
  mutate(rownames = rownames) |> 
  rownames_to_column("old_names") |> 
  column_to_rownames("rownames") |> 
  select(-old_names)
# -----------------------------------------------------
head(sample_table)
mt_mag <-microtable$new(otu_table = otu_table, sample_table =$sample_table, 
                        tax_table = tax_table, phylo_tree = tree)
mt_mag$sample_table$strat_season <- factor(mt_mag$sample_table$strat_season, 
  levels = c("summer", "fall", "winter", "spring"), 
  labels = c("Summer", "Fall", "Winter", "Spring"),
   ordered = T)
# mt_mag <- trans_norm$new(mt_mag)
# mt_mag <- mt_mag$norm(method = "rclr")
# Create a clone incase I mess up the original

mt_mag$tidy_dataset()
mt_mag
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

# Great, now put all those bins in a vector
con_vec <- unique(con_df$bin)

# Create a column to show which bins are potential contamination

mt_mag$tax_table <- mt_mag$tax_table |> 
  mutate(blank_contam = ifelse(mt_mag$tax_table$bin %in% paste0("b__", con_vec), TRUE, FALSE)
  )

view(mt_mag$tax_table)

mt_mag

mt_mag$cal_abund()

n <- "50"
rank <- "bin"

abund <-trans_abund$new(mt_mag, ntaxa = n, taxrank = rank, high_level = "poss_contaminent")$
  plot_bar(others_color = "grey70", ggnested = T, xtext_keep = TRUE, legend_text_italic = FALSE, 
facet = "strat_season") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=9)) +
  ggtitle("MAGs Data", subtitle = paste(n, "most abundant", rank, "green, possible contaminent from 16s data,", 
"brown, not possible contaminent")) +
  theme(legend.position = "none")
abund
plotly::plotly_build(abund)

abund <-trans_abund$new(mt_mag, ntaxa = n, taxrank = rank, high_level = "blank_contam")$
  plot_bar(others_color = "grey70", ggnested = T, xtext_keep = TRUE, legend_text_italic = FALSE, 
facet = "strat_season") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=9)) +
  ggtitle("MAGs Data", subtitle = paste(n, "most abundant", rank, "green, possible contaminent from blanks data,", 
"brown, not possible contaminent")) +
  theme(legend.position = "none")
abund
plotly::plotly_build(abund)

# Based on the information from the blanks, I'm going to remove the taxa that are present in the blank samples. 

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



load("output/data/mt_mag_large.RData")

# Expand foam_rollup function 

foam_rollup_df <- foam_rollup_df |> 
  # rename(func = `function`) |> 
  separate_wider_delim(col = func, delim = ";", names = c("abbrev", "long_name"))

foam_rollup_df <- foam_rollup_df |> 
  separate_wider_delim(abbrev, delim = ",", 
                       names = c("abbrev_1", "abbrev_2", "abbrev_3", "abbrev_4",
                                 "abbrev_5", "abbrev_6"), 
                       too_few = "align_start")
genes <- colnames(mt_mag_large$tax_table)[162:length(colnames(mt_mag_large$tax_table))]


foam_rollup_df <- foam_rollup_df |> 
  mutate(across(matches("abbrev.*"), ~ .x %in% genes, .names = "in_mt_{col}"))

foam_rollup_df <- foam_rollup_df |> 
  mutate(in_mt = if_any((matches("in_mt_abb.*"))))

