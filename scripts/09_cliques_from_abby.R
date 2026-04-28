# Abby Code

#using but modifying the code from fullerton et al. BMS_16S_final_analysis so that i can use microeco

### Load required libraries
#library(microbiome) # data analysis and visualisation. this is discontinued but is an extension of phyloseq so may not be needed
library(microeco) # also the basis of data object. Data analysis and visualisation
library(RColorBrewer) # nice color options
library(dplyr) # data handling
library(network) # networks
library(intergraph)  # networks
library(ggnetwork)   # network plotting with ggplot. I updated this from ggnet to ggnetwork, I hope the functions are similar.
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


#they start with loading in datasets and tree. I'm going to load in my microtable and separate it out into the tax and counts and see what happens.
mt_16s <- readRDS("/Users/abbysmason/SeaGrant2024/SeaGrant2024/outputs/estuary_microtableall.rds")

load("output/data/mt_16s.RData")
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
  rel_abund = 0.0002,  # At least 0.2% in any sample
)

#separate out taxa and counts here
abby_tax <- mt_16s$tax_table
abby_counts <- t(mt_16s$otu_table)

#line 979:
#next they make the co-occurence matrix. they use igraph, this is part of the trans_diff function in microeco
# so i need to create a new trans_diff object.

# Retaining edges with Spearman correlation > 0.65
t_network <- trans_network$new(dataset = mt_16s,
                               cor_method = "spearman",
                               use_corr_p_adjust = FALSE)

#filtering correlations
#this step corresponds to `bac.cor[bac.cor < 0.7] = 0` in the original script
t_network$cal_network(COR_p_thres = 0.01, COR_cut = 0.7)


# Calculate modules/cliques using the Louvain algorithm
# This replaces the `cluster_louvain(bac.cor.ig)` function call
t_network$cal_module(method = "cluster_louvain")

#success? 5 modules. idk if thats good
modules <- igraph::vertex_attr(t_network$res_network, "module")
head(modules)
table(modules)

all_taxa <- rownames(mt_16s$otu_table)

# Extract taxa for each module
clique1_taxa <- all_taxa[modules == "M1"]
clique2_taxa <- all_taxa[modules == "M2"]
clique3_taxa <- all_taxa[modules == "M3"]
clique4_taxa <- all_taxa[modules == "M4"]
clique5_taxa <- all_taxa[modules == "M5"]
clique6_taxa <- all_taxa[modules == "M6"]

module1_taxonomy <- as.data.frame(mt_16s$tax_table[clique1_taxa, ])
module2_taxonomy <- as.data.frame(mt_16s$tax_table[clique2_taxa, ])
module3_taxonomy <- as.data.frame(mt_16s$tax_table[clique3_taxa, ])
module4_taxonomy <- as.data.frame(mt_16s$tax_table[clique4_taxa, ])
module5_taxonomy <- as.data.frame(mt_16s$tax_table[clique5_taxa, ])
module6_taxonomy <- as.data.frame(mt_16s$tax_table[clique6_taxa, ])


#plotting network
V(t_network$res_network)$color <- modules

#make layout
layout <- layout_with_fr(t_network$res_network)

module_colors <- c("M1" = "red",
                   "M2" = "blue",
                   "M3" = "green",
                   "M4" = "purple",
                   "M5" = "orange",
                   "M6" = "yellow")

#map modules vector to colors
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



#ok.... how about summing up the abundance of the modules and plotting them over time?
#i need to figure out the best way to do relative abundance - there's no clear way with cal_abund

otu_matrix <- mt_16s$otu_table  # Taxa as rows, samples as columns

# Convert to relative abundance (percentage)
rel_abund_otu <- apply(otu_matrix, 2, function(x) x / sum(x) * 100)

# Now sum by module
module_rel_abund <- data.frame(row.names = colnames(otu_matrix))

for(mod in unique(modules)) {
  taxa_in_mod <- names(modules[modules == mod])
  module_rel_abund[, mod] <- colSums(rel_abund_otu[taxa_in_mod, , drop = FALSE])
}

head(module_rel_abund)

# Add time variable (replace 'Time' with your actual column name)
module_rel_abund$date <- mt_16s$sample_table$date
module_rel_abund$strat_season <- mt_16s$sample_table$strat_season
# module_rel_abund$Site <- mt_16s$sample_table$Site

module_long <- pivot_longer(module_rel_abund,
                            cols = c("M1", "M2", "M3", "M4", "M5", "M6"),
                            names_to = "Module",
                            values_to = "relabund")


# Create dates for season starts
match("Summer", module_long$strat_season)
sum_start <- module_long$date[match("Summer", module_long$strat_season)]
fall_start <- module_long$date[match("Fall", module_long$strat_season)]
winter_start <- module_long$date[match("Winter", module_long$strat_season)]
spring_start <- module_long$date[match("Spring", module_long$strat_season)]
# Line plot - all regions
ggplot(module_long, aes(x = date, y = relabund, color = Module, group = Module, label = rownames(module_long))) +
  geom_rect(aes(xmin = sum_start, xmax = fall_start, ymax = max(relabund) +4 , ymin = max(relabund)+1),  fill = "grey", color = NULL) +
  geom_rect(aes(xmin = fall_start, xmax = winter_start, ymax = max(relabund) +4 , ymin = max(relabund)+1),  fill = "green4",  color = NULL) +
  geom_rect(aes(xmin = winter_start, xmax = spring_start, ymax = max(relabund) +4 , ymin = max(relabund)+1),  fill = "grey",  color = NULL) +
  geom_rect(aes(xmin = spring_start, xmax = max(date), ymax = max(relabund) +4 , ymin = max(relabund)+1),  fill = "green4",  color = NULL) +
  geom_line(size = 1) +
  geom_text() +
  geom_point(size = 2) +
  labs(title = "Module Abundance Over Time",
       y = "Relative Abundance") +
  theme_bw()
# Add annotation for seasons in ppt


#maybe just one region over time with bars

# ggplot(module_long, aes(x = as.factor(Time), y = relabund, fill = Sampling.Region)) +
#   geom_bar(stat = "identity", position = "stack") +
#   facet_wrap(~Module, scales = "free_y")


#the problem with this is that it shows all data points from each site which are all very different. I'm going to do boxplots with averages and see if that looks better

module_means <- module_long %>%
  group_by(Sampling.Region, Time, Site, Module) %>%
  dplyr::summarise(Mean_Abund = mean(relabund, na.rm = TRUE))

#so that each line only goes through averages, not every point
line_means <- module_long %>%
  group_by(Sampling.Region, Time, Module) %>%
  dplyr::summarise(
    Mean_Abund = median(relabund, na.rm = TRUE),
    .groups = 'drop'
  )

module_long_estuary <- module_long %>%
  filter(!Sampling.Region %in% c("Apostle Islands", "Jay Cooke State Park", "Pokegama River"))
 
line_means_estuary <- line_means %>%
  filter(!Sampling.Region %in% c("Apostle Islands", "Jay Cooke State Park", "Pokegama River"))

module_means_estuary <- module_means %>%
  filter(!Sampling.Region %in% c("Apostle Islands", "Jay Cooke State Park", "Pokegama River"))
 
ggplot(module_long_estuary, aes(x = factor(Time), y = relabund, fill = Module)) +
  #geom_boxplot(alpha = 1, outlier.size = 0.5) +
  # Line connecting the means
  geom_line(data = line_means_estuary,
            aes(x = factor(Time), y = Mean_Abund, group = Module, color = Module),
            size = 1.2, inherit.aes = FALSE) +
  # Points at the means
  #geom_point(data = module_means_estuary,
             #aes(x = factor(Time), y = Mean_Abund, color = Module),
             #size = 0.5, inherit.aes = FALSE) +
  facet_wrap(~Sampling.Region) +
  labs(y = "Relative Abundance (%)", x = "Time Point") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))