# Load required packages
library(tidyverse)
library(vegan)
library(Polychrome)
library(forcats)

# Read in the data
b_df <- readr::read_tsv("data/superior_sylph/Sylph_TaxAbund_out.tsv")

# Create vectors of column names
clade_cols <-
  stringr::str_extract_all(b_df$clade_name[8], "[a-z]{1}(?=_)")
clade_cols <- as.vector(clade_cols[[1]]) %>% 
  print()

site_cols <- stringr::str_extract(colnames(b_df), "(?<=/).{4}") %>% 
  na.omit() %>% 
  print()
colnames(b_df)[2:ncol(b_df)] <- site_cols

#Spread columns out, clean names, 

b_df <- b_df %>% 
  tidyr::separate_wider_delim(cols = clade_name, delim = "|", 
                       names = clade_cols, too_few = "align_start") %>% 
  # Drop all rows that don't get down to a taxa
  tidyr::drop_na("t") %>% 
  #remove additional characters and make data columns numeric
  mutate(across(clade_cols, ~ str_remove(., "[a-z]__")))

  # Quick check to see if abundances are all still 100
  colSums(b_df[9:ncol(b_df)])
  
#Create long df
long_b_df <- b_df %>% 
  pivot_longer(cols = all_of(site_cols),
               names_to = "sample", values_to = "abundance")

### Create bar chart

  #First create a unique palette based on how many catergories there are
palette <- length(unique(b_df$p)) %>% 
  glasbey.colors() %>% 
  unname()


long_b_df$p <- fct_reorder(long_b_df$p, long_b_df$abundance)

  # Plot
long_b_df %>% 
  ggplot(aes(x = sample, y = abundance, fill = p)) +
  geom_bar(stat = "identity") +
  # scale_fill_discrete() +
  scale_fill_manual(values = palette) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  NULL

## NMDS 
  b_dist <-
    b_df %>% 
    select(all_of(site_cols)) %>% 
    dist()
  
  b_dist
  nmds <- b_df %>% 
    select(all_of(site_cols)) %>% 
    t() %>% 
    # decostand("hellinger") %>% 
    metaMDS(distance="bray", autotransform=FALSE, binary=FALSE, 
                  noshare=TRUE, zerodist="ignore")
  
  nmds
  
  # Extract scores
  
  nmds_scores <- scores(nmds, display = "sites") %>%
    data.frame() %>% 
    rownames_to_column(var = "sample")
  
  # Basic NMDS plot
  nmds_basic_plot <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = sample)) +
    geom_point(size = 3, alpha = 0.8) +
    labs(title = "NMDS Ordination of WM Data",
         subtitle = paste("Stress =", round(nmds$stress, 3)),
         x = "NMDS1", 
         y = "NMDS2") +
    theme_minimal()
  
  nmds_basic_plot
  
  # Run stressplot
  stressplot(nmds)


