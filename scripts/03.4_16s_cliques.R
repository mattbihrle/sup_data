# Abby Code

#using but modifying the code from fullerton et al. BMS_16S_final_analysis so that i can use microeco

### Load required libraries
#library(microbiome) # data analysis and visualisation. this is discontinued but is an extension of phyloseq so may not be needed
library(microeco) # also the basis of data object. Data analysis and visualisation
library(RColorBrewer) # nice color options
library(dplyr) # data handling
library(network) # networks
library(intergraph)  # networks
# library(ggnetwork)   # network plotting with ggplot. I updated this from ggnet to ggnetwork, I hope the functions are similar.
library(igraph)  # networks. leaving this out bc its automatically loaded with microeco
library(ggplot2) # plotting library
library(gridExtra) # gridding plots
library(ape) # importing and handling phylogenetic trees
library(ggthemes) # additional themes fro ggplot2
library(magrittr) #
library(rioja) # plotting poackages for tabular bubbleplots
library(ggpubr)
library(ggtern) # ternary plots for geochemistry
library(plyr)
library(coda.base)
#library(tydiverse) is this supposed to be tidyverse??? this was in the og code i copied
library(tidyverse)
library(vegan) # Multivariate ecological analysis
#library(propr) #also not available for my version of R. leaving out for now.
library(missForest) # Imputing missing values in dataframes using Random Forests
library(VSURF) # Random Forests approach to variable importance identification
library(patchwork)


#they start with loading in datasets and tree. I'm going to load in my microtable and separate it out into the tax and counts and see what happens.
# mt_16s <- readRDS("/Users/abbysmason/SeaGrant2024/SeaGrant2024/outputs/estuary_microtableall.rds")

load("output/data/mt_16s.RData")

load("output/data/temp_df_long.RData")
#I'm just going to use my normal microtable that has everything and then use the contaminants list they have instead of just using cyanos
#if i only did cyanos i would revisit the contaminant list and see what's still applicable


###CONTAMINATION SCREENING###
# List of potential contaminant genera in subsurface 16S rRNA libraries after Sheik et al. 2018 Frontiers in Microbiology
# contaminants <- c("Afipia", "Aquabacterium", "Asticcacaulis", "Aurantimonas", "Beijerinckia", "Bosea", "Bradyrhizobium", "Brevundimonas", "Caulobacter", "Craurococcus", "Devosia", "Hoefleae", "Mesorhizobium", "Methylobacterium", "Novosphingobium", "Ochrobactrum", "Paracoccus", "Pedomicrobium", "Phyllobacterium", "Rhizobium", "Roseomonas", "Sphingobium", "Sphingomonas", "Sphingopyxis", "Acidovorax", "Azoarcus", "Azospira", "Burkholderia", "Comamonas", "Cupriavidus", "Curvibacter", "Delftiae", "Duganella", "Herbaspirillum", "Janthinobacterium", "Kingella", "Leptothrix", "Limnobacter", "Massilia", "Methylophilus", "Methyloversatilis", "Neisseria", "Oxalobacter", "Pelomonas", "Polaromonas", "Ralstonia", "Schlegelella", "Sulfuritalea", "Undibacterium", "Variovorax", "Acinetobactera", "Enhydrobacter", "Enterobacter", "Escherichia", "Nevskia", "Pasteurella", "Pseudomonas", "Pseudoxanthomonas", "Psychrobacter", "Stenotrophomonas", "Xanthomonas", "unclassified Acidobacteria Gp2", "Aeromicrobium", "Actinomyces", "Arthrobacter", "Beutenbergia", "Brevibacterium", "Corynebacterium", "Curtobacterium", "Dietzia", "Geodermatophilus", "Janibacter", "Kocuria", "Microbacterium", "Micrococcus", "Microlunatus", "Patulibacter", "Propionibacterium", "Rhodococcus", "Tsukamurella", "Chryseobacterium", "Dyadobacter", "Flavobacterium", "Hydrotalea", "Niastella", "Olivibacter", "Parabacteroides", "Pedobacter", "Prevotella", "Wautersiella", "Deinococcus", "Abiotrophia", "Bacillus", "Brevibacillus", "Brochothrix", "Facklamia", "Lactobacillus", "Paenibacillus", "Ruminococcus", "Staphylococcus", "Streptococcus", "Veillonella", "Fusobacterium")

#i'm just going to use microeco's built-in filtering of contaminants
# mt_16s$filter_pollution(taxa = contaminants)

#Remove mitochondria and chloroplasts. this is the default option for filter_polution so i don't need to put the taxa in there
# mt_16s$filter_pollution()

#stopped at line 227. going to skip to cliques
#filter by prevalence
mt_16s$filter_taxa(
  rel_abund = 0.0001,  # At least 0.1% in any sample
)

mt_16s
mt_16s$tidy_dataset()
mt_16s


# #separate out taxa and counts here
# abby_tax <- mt_16s$tax_table
# abby_counts <- t(mt_16s$otu_table)

#line 979:
#next they make the co-occurence matrix. they use igraph, this is part of the trans_diff function in microeco
# so i need to create a new trans_diff object.

# Retaining edges with Spearman correlation > 0.65
t_network <- trans_network$new(dataset = mt_16s,
                               cor_method = "spearman",
                               use_corr_p_adjust = FALSE)

#filtering correlations
#this step corresponds to `bac.cor[bac.cor < 0.7] = 0` in the original script
t_network$cal_network(COR_p_thres = 0.05, COR_optimization = TRUE)


# Calculate modules/cliques using the Louvain algorithm
# This replaces the `cluster_louvain(bac.cor.ig)` function call
t_network$cal_module(method = "cluster_louvain")

#success? 7 modules. idk if thats good
modules <- igraph::vertex_attr(t_network$res_network, "module")
head(modules)
table(modules)

test <- t_network$res_network[[1]]
#plotting network
V(t_network$res_network)$color <- modules

#make layout
layout <- layout_with_fr(t_network$res_network)

module_colors <- c("M1" = "red",
                   "M2" = "orange",
                   "M3" = "green",
                   "M4" = "purple",
                   "M5" = "blue",
                   "M6" = "yellow",
                   "M7" = "pink",
                   "M8" = "brown")

# map modules vector to colors
vertex_colors <- module_colors[modules]

# Maybe add sample names here?
plot(t_network$res_network,
     layout = layout,
     vertex.size = 5,
     vertex.color = vertex_colors,
     vertex.label = NA,
     edge.width = 0.5,
     edge.color = "gray80")


#giving names to module vector
modules <- igraph::vertex_attr(t_network$res_network, "module")
taxa_names <- igraph::vertex_attr(t_network$res_network, "name")

# Give names to the modules vector
names(modules) <- taxa_names
modules

#ok.... how about summing up the abundance of the modules and plotting them over time?
# #i need to figure out the best way to do relative abundance - there's no clear way with cal_abund
# 
# # First clean out the mt_16s table
# 
# otu_matrix <- mt_16s$otu_table  # Taxa as rows, samples as columns
# 
# # Convert to relative abundance (percentage)
# rel_abund_otu <- apply(otu_matrix, 2, function(x) x / sum(x) * 100)
# 
# # Now sum by module
# module_rel_abund <- data.frame(row.names = colnames(otu_matrix))
# 
# for(mod in unique(modules)) {
#   taxa_in_mod <- names(modules[modules == mod])
#   module_rel_abund[, mod] <- colSums(rel_abund_otu[taxa_in_mod, , drop = FALSE])
# }
# 
# head(module_rel_abund)
# 
# # Add time variable (replace 'Time' with your actual column name)
# module_rel_abund$date <- mt_16s$sample_table$date
# module_rel_abund$strat_season <- mt_16s$sample_table$strat_season
# # module_rel_abund$Site <- mt_16s$sample_table$Site
# 
# module_long <- pivot_longer(module_rel_abund,
#                             cols = c("M1", "M2", "M3", "M4", "M5", "M6", "M7"),
#                             names_to = "Module",
#                             values_to = "relabund")
# save(module_long, file = "output/data/module_long_16s.RData")

# Add modules to tax_table in 16s data

# first turn modules vector into df
modules_df <- data.frame(modules) |> 
  rownames_to_column("otu")
mt_16s$tax_table <- mt_16s$tax_table |>  
  rownames_to_column("otu_rows") |> 
  left_join(modules_df, by = join_by(otu_rows == otu)) |>
  column_to_rownames("otu_rows") |> 
  mutate(otu = rownames(mt_16s$tax_table)) 


mt_16s$tax_table <- mt_16s$tax_table |> 
  mutate(mod_num = str_pad(str_remove(modules, "M"), width = 2, side = "left", pad = "0")) |> 
  mutate(modules = ifelse(mod_num == "NA", NA, paste0("M", mod_num))) |> 
  select(-mod_num)

mt_16s$tax_table$modules <- factor(mt_16s$tax_table$modules)

mt_16s$cal_abund()
save(mt_16s, file = "output/data/mt_16s_cliques.RData")
save(t_network, file = "output/data/t_network_16s.RData")
load("output/data/mt_16s_cliques.RData")

mt_16s$tax_table <- mt_16s$tax_table |> 
  mutate(mod_num = str_pad(str_remove(modules, "M"), width = 2, side = "left", pad = "0")) |> 
  mutate(modules = ifelse(mod_num == "NA", NA, paste0("M", mod_num))) |> 
  select(-mod_num)

mt_16s$tax_table$modules <- factor(mt_16s$tax_table$modules)

mt_16s$tax_table$modules
mt_16s
# Look at specifically M1-----------------------------------------------------
mt_m1 <- clone(mt_16s)

mt_m1$tax_table <- mt_m1$tax_table |> 
  filter(modules == "M06")

mt_m1$tidy_dataset()$cal_abund()
# Plot abundances

rank = "g"
ntaxa = 30
# Get vector of 20 most abundant families

trans_abund$new(dataset = mt_m1, taxrank = rank, ntaxa = ntaxa)$
  plot_bar(facet = "strat_season") +
  ggtitle("16S", subtitle = paste("top", ntaxa, rank, "in M06"))


# Look at specifically M2
mt_m2 <- clone(mt_16s)

mt_m2$tax_table <- mt_m2$tax_table |> 
  filter(modules == "M01")

mt_m2$tidy_dataset()$cal_abund()
# Plot abundances

rank = "g"
ntaxa = 30
# Get vector of 20 most abundant families

trans_abund$new(dataset = mt_m2, taxrank = rank, ntaxa = ntaxa)$
  plot_bar(facet = "strat_season") +
  ggtitle("16S", subtitle = paste("top", ntaxa, rank, "in M02"))

# Look at specifically M3
mt_m3 <- clone(mt_16s)

mt_m3$tax_table <- mt_m3$tax_table |> 
  filter(modules == "M3")

mt_m3$tidy_dataset()$cal_abund()
# Plot abundances

rank = "g"
ntaxa = 30
# Get vector of 20 most abundant families

trans_abund$new(dataset = mt_m3, taxrank = rank, ntaxa = ntaxa)$
  plot_bar(facet = "strat_season") +
  ggtitle("16S", subtitle = paste("top", ntaxa, rank, "in M3"))

# Look at specifically M4
mt_m4 <- clone(mt_16s)

mt_m4$tax_table <- mt_m4$tax_table |> 
  filter(modules == "M4")

mt_m4$tidy_dataset()$cal_abund()
# Plot abundances

rank = "otu"
ntaxa = 30
# Get vector of 20 most abundant families

 trans_abund$new(dataset = mt_m4, taxrank = rank, ntaxa = ntaxa, high_level = "g")$
  plot_bar(facet = "strat_season", ggnested = TRUE) +
  ggtitle("16S", subtitle = paste("top", ntaxa, rank, "in M4"))



# Plot module abundances over time ------------------------------------------------------

t1 <- trans_abund$new(mt_16s, taxrank = "modules", ntaxa = 7)
plot_line(t1)
mt_16s$tax_table$modules |> 
  table()
t1$data_abund |> 
  filter(Taxonomy %in% unique(t1$data_abund$Taxonomy)[1:7]) |> 
  ggplot(aes(x = date, y = Abundance, color = Taxonomy, label = Taxonomy)) +
    # geom_rect(aes(xmin = sum_start, xmax = fall_start, ymax = max(relabund) +4 , ymin = max(relabund)+1),  fill = "grey", color = NULL) +
    # geom_rect(aes(xmin = fall_start, xmax = winter_start, ymax = max(relabund) +4 , ymin = max(relabund)+1),  fill = "green4",  color = NULL) +
    # geom_rect(aes(xmin = winter_start, xmax = spring_start, ymax = max(relabund) +4 , ymin = max(relabund)+1),  fill = "grey",  color = NULL) +
    # geom_rect(aes(xmin = spring_start, xmax = max(date), ymax = max(relabund) +4 , ymin = max(relabund)+1),  fill = "green4",  color = NULL) +
    geom_vline(aes(xintercept = date, color= strat_season), alpha = 0.3, linewidth = 8) +
    geom_line(size = 1) +
      geom_text() +
        geom_point(size = 2) +
  scale_color_brewer(palette = "Paired") +
  # scale_fill_discrete() +
  labs(title = "Module Abundance Over Time",
       y = "Relative Abundance") +
  theme_bw()

p2 <- ggplot(module_long, aes(x = date, y = relabund, color = module, group = module, label = rownames(module_long))) +
  # geom_rect(aes(xmin = sum_start, xmax = fall_start, ymax = max(relabund) +4 , ymin = max(relabund)+1),  fill = "grey", color = NULL) +
  # geom_rect(aes(xmin = fall_start, xmax = winter_start, ymax = max(relabund) +4 , ymin = max(relabund)+1),  fill = "green4",  color = NULL) +
  # geom_rect(aes(xmin = winter_start, xmax = spring_start, ymax = max(relabund) +4 , ymin = max(relabund)+1),  fill = "grey",  color = NULL) +
  # geom_rect(aes(xmin = spring_start, xmax = max(date), ymax = max(relabund) +4 , ymin = max(relabund)+1),  fill = "green4",  color = NULL) +
  # geom_line(size = 1) +
  geom_text() +
  geom_point(size = 2) +
  labs(title = "Module Abundance Over Time",
       y = "Relative Abundance") +
  theme_bw()
# Add annotation for seasons in ppt

p2


