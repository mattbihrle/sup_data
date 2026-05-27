### Load required libraries
library(microbiome) # data analysis and visualisation
library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(RColorBrewer) # nice color options
library(dplyr) # data handling
library(network) # networks
library(intergraph)  # networks
# library(ggnetwork)   # network plotting with ggplot
library(igraph)  # networks
library(phyloseq) # ASV ecological analysis package
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
# library(tydiverse)
library(vegan) # Multivariate ecological analysis
library(propr)
library(missForest) # Imputing missing values in dataframes using Random Forests
library(VSURF) # Random Forests approach to variable importance identification
library(car) #for scatterplot
library(microeco)
library(tidyverse)
library(patchwork)


# Load in required data

load("output/data/mt_gene_mag.RData")

mt_gene_mag$filter_taxa(
  rel_abund = 0.00001
)

tax <- mt_gene_mag$tax_table

counts <- t(mt_gene_mag$otu_table)


# Retaining edges with Spearman correlation > 0.65
t_network <- trans_network$new(dataset = mt_gene_mag,
                               cor_method = "spearman",
                               use_corr_p_adjust = FALSE)
t_safe_network <- clone(t_network)
save(t_network, file = "output/data/t_network.RData")

# Start here 
load("output/data/t_network.RData")

#filtering correlations
#this step corresponds to `bac.cor[bac.cor < 0.7] = 0` in the original script
t_network$cal_network(COR_p_thres = 0.01, COR_cut = 0.5)

t_network$cal_module(method = "cluster_louvain")
save(t_network, file = "output/data/t_network.RData")
modules <- igraph::vertex_attr(t_network$res_network, "module")
view(modules)
table(modules)
all_taxa <- rownames(mt_gene_mag$otu_table)s
view(all_taxa)
# Extract taxa for each module
clique1_taxa <- all_taxa[modules == "M1"]
clique2_taxa <- all_taxa[modules == "M2"]
clique3_taxa <- all_taxa[modules == "M3"]
clique4_taxa <- all_taxa[modules == "M4"]
clique5_taxa <- all_taxa[modules == "M5"]

view(clique3_taxa)

V(t_network$res_network)$color <- modules

layout <- layout_with_fr(t_network$res_network)

module_colors <- c("M1" = "red",
                   "M2" = "orange",
                   "M3" = "green",
                   "M4" = "blue")

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

otu_matrix <- mt_gene_mag$otu_table  # Taxa as rows, samples as columns

# FIRST SEE IF THESE AREN'T ALREADY RELATIVE ABUNDANCE
# --------------------------------------------------------------
# Convert to relative abundance (percentage)
rel_abund_otu <- apply(otu_matrix, 2, function(x) x / sum(x) * 100)

# Now sum by module
module_rel_abund <- data.frame(row.names = colnames(otu_matrix))

for(mod in unique(modules)) {
  taxa_in_mod <- names(modules[modules == mod])
  module_rel_abund[, mod] <- colSums(rel_abund_otu[taxa_in_mod, , drop = FALSE])
}

head(module_rel_abund)
# ---------------------------------------------------------------------

# Add time variable (replace 'Time' with your actual column name)
module_rel_abund$date <- mt_gene_mag$sample_table$date
module_rel_abund$strat_season <- mt_gene_mag$sample_table$strat_season
# module_rel_abund$Site <- mt_gene_mag$sample_table$Site

module_long <- pivot_longer(module_rel_abund,
                            cols = c("M1", "M2", "M3", "M4"),
                            names_to = "Module",
                            values_to = "relabund")

# Add modules to tax_table in 16s data

  # first turn modules vector into df
modules_df <- data.frame(modules) |> 
  rownames_to_column("otu")
mt_gene_mag$tax_table <- mt_gene_mag$tax_table |> 
  rownames_to_column("otu") |> 
  left_join(modules_df) |>
  column_to_rownames("otu") |> 
  mutate(otu = rownames(mt_gene_mag$tax_table))

mt_gene_mag$cal_abund()

# Look at specifically M1
mt_m1 <- clone(mt_gene_mag)

mt_m1$tax_table <- mt_m1$tax_table |> 
  filter(modules == "M1")

mt_m1$tidy_dataset()$cal_abund()
# Plot abundances

rank = "gene"
ntaxa = 30
# Get vector of 20 most abundant families

trans_abund$new(dataset = mt_m1, taxrank = rank, ntaxa = ntaxa)$
  plot_bar(facet = "strat_season") +
     ggtitle("16S", subtitle = paste("top", ntaxa, rank, "in M1"))


# Look at specifically M2
mt_m2 <- clone(mt_gene_mag)

mt_m2$tax_table <- mt_m2$tax_table |> 
  filter(modules == "M2")

mt_m2$tidy_dataset()$cal_abund()
# Plot abundances

rank = "otu"
ntaxa = 30
# Get vector of 20 most abundant families

trans_abund$new(dataset = mt_m2, taxrank = rank, ntaxa = ntaxa)$
  plot_bar(facet = "strat_season") +
     ggtitle("16S", subtitle = paste("top", ntaxa, rank, "in M2"))

# Look at specifically M3
mt_m3 <- clone(mt_gene_mag)

mt_m3$tax_table <- mt_m3$tax_table |> 
  filter(modules == "M3")

mt_m3$tidy_dataset()$cal_abund()
# Plot abundances

rank = "otu"
ntaxa = 30
# Get vector of 20 most abundant families

trans_abund$new(dataset = mt_m3, taxrank = rank, ntaxa = ntaxa)$
  plot_bar(facet = "strat_season") +
     ggtitle("16S", subtitle = paste("top", ntaxa, rank, "in M3"))

# Look at specifically M4
mt_m4 <- clone(mt_gene_mag)

mt_m4$tax_table <- mt_m4$tax_table |> 
  filter(modules == "M4")

mt_m4$tidy_dataset()$cal_abund()
# Plot abundances

rank = "otu"
ntaxa = 30
# Get vector of 20 most abundant families

trans_abund$new(dataset = mt_m4, taxrank = rank, ntaxa = ntaxa)$
  plot_bar(facet = "strat_season") +
     ggtitle("16S", subtitle = paste("top", ntaxa, rank, "in M4"))

# Create dates for season starts
match("Summer", module_long$strat_season)
sum_start <- module_long$date[match("Summer", module_long$strat_season)]
fall_start <- module_long$date[match("Fall", module_long$strat_season)]
winter_start <- module_long$date[match("Winter", module_long$strat_season)]
spring_start <- module_long$date[match("Spring", module_long$strat_season)]
# Line plot - all regions
p2 <- ggplot(module_long, aes(x = date, y = relabund, color = Module, group = Module, label = rownames(module_long))) +
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
p2

# Add temp labels
temp_labels_1 <- seq(0, 17, by = 1)
temp_labels_2 <- seq(1, 18, by = 1)

temp_labels <- paste0(temp_labels_1, "-", temp_labels_2)

# Add samp dates


## Test plotting cliques over the temperature
load("output/data/temp_df_long_full.RData")
temp_df_long_full |> 
    ungroup() |> 
    # select(dttm, depth, temp) |> 
    drop_na(temp) |> 
    mutate(depth = as.numeric(depth)) |> 
    # filter(depth < 50) |> 
        dplyr::select(date, depth, temp, sample_date) |> 
        # drop_na(temp) |> 
        distinct() |> 
    ggplot(aes(x = date, y = depth, z = temp)) +
      # geom_contour_filled(breaks = breaks) +
      geom_contour_filled(binwidth = 1) +
        # scale_fill_brewer(palette = "RdBu", direction = -1) +
        scale_fill_viridis_d(option = "turbo") +
  labs(x = "Date", y = "Depth (m)", fill = "Temperature (°C)") +
    theme_classic() +
    theme(
      text = element_text(size = 15),
      panel.border = element_rect(colour = "black", fill = "red", linewidth = 3)
    ) +
  scale_y_reverse(expand = c(0,0), n.breaks = 10) +
  scale_x_date(expand = c(0,0), date_breaks = "1 month", date_label = "%b", limits = c(as_date("2024-07-20"), NA), ) +
    guides(fill = guide_legend(ncol = 1, reverse = T)) +
  # Clique Data
  geom_point(data = module_long, mapping = aes(x = date, y = 160 - relabund*1.5, z = NULL, color = Module, group = Module, label = rownames(module_long)), size = 2) +
  geom_line(data = module_long, mapping = aes(x = date, y = 160 - relabund*1.5, z = NULL, color = Module, group = Module, label = rownames(module_long), size = 1)) +
  scale_color_viridis_d(option = "turbo") +
    # geom_point(data = samp_dates, aes(x = sample_date, y = 38, z = NULL), color = 'white') +
    # geom_point(data = samp_dates, aes(x = sample_date, y = 22, z = NULL), color = 'white') +
    # geom_hline(yintercept = 22, color = 'white', linetype = "dashed" ) +
    # geom_vline(xintercept = sw_meta$date[which(sw_meta$strat_season == "summer")], color = "red") +
    # geom_vline(xintercept = sw_meta$date[which(sw_meta$strat_season == "fall")]) +
    # geom_vline(xintercept = sw_meta$date[which(sw_meta$strat_season == "winter")], , color = "red") +
    # geom_vline(xintercept = sw_meta$date[which(sw_meta$strat_season == "spring")]) +
    theme(text = element_text(size = 25), aspect.ratio = 0.7) +
  guides(size = "none") +
    NULL




p1 <- temp_df_long_full |> 
    ungroup() |> 
    # select(dttm, depth, temp) |> 
    drop_na(temp) |> 
    mutate(depth = as.numeric(depth)) |> 
    filter(depth < 50) |> 
        dplyr::select(date, depth, temp) |> 
        # drop_na(temp) |> 
        distinct() |> 
    ggplot(aes(x = date, y = depth, z = temp)) +
      # geom_contour_filled(breaks = breaks) +
      geom_contour_filled(binwidth = 1) +
        # scale_fill_brewer(palette = "RdBu", direction = -1) +
        scale_fill_viridis_d(option = "turbo", labels = temp_labels) +
  labs(x = "Date", y = "Depth (m)", fill = "Temperature (°C)") +
  geom_point(data = mt_gene_mag$sample_table, aes(x = date, y = 38, z = NULL), color = 'white') +
    theme_classic() +
    theme(
      text = element_text(size = 15),
      panel.border = element_rect(colour = "black", fill = "red", linewidth = 3)
    ) +
  scale_y_reverse(expand = c(0,0), n.breaks = 10) +
  scale_x_date(expand = c(0,0), date_breaks = "1 month", date_label = "%b", limits = c(as_date("2024-07-10"), as_date("2025-08-05"))) +
    guides(fill = guide_legend(ncol = 1, reverse = T)) 
p1
p2 <- ggplot(module_long, aes(x = date, y = relabund, color = Module, group = Module, label = rownames(module_long))) +
  # geom_rect(aes(xmin = sum_start, xmax = fall_start, ymax = max(relabund) +4 , ymin = max(relabund)+1),  fill = "grey", color = NULL) +
  # geom_rect(aes(xmin = fall_start, xmax = winter_start, ymax = max(relabund) +4 , ymin = max(relabund)+1),  fill = "green4",  color = NULL) +
  # geom_rect(aes(xmin = winter_start, xmax = spring_start, ymax = max(relabund) +4 , ymin = max(relabund)+1),  fill = "grey",  color = NULL) +
  # geom_rect(aes(xmin = spring_start, xmax = max(date), ymax = max(relabund) +4 , ymin = max(relabund)+1),  fill = "green4",  color = NULL) +
  geom_line(size = 1) +
  geom_text() +
  geom_point(size = 2) +
  labs(title = "Module Abundance Over Time",
       y = "Relative Abundance") +
  theme_bw() +
  scale_x_date(expand = c(0,0), date_breaks = "1 month", date_label = "%b", limits = c(as_date("2024-07-10"), as_date("2025-08-05")))

p2

p3 <- maestro_df_clean |> 
  select(dttm_local, matches(c("do_sat_12"))) |> 
  pivot_longer(cols = -dttm_local, names_to = "var", values_to = "value") |> 
    ggplot(aes(dttm_local, value)) +
  geom_point() + 
  facet_grid(var ~ ., scales = "free") +
  scale_x_date(expand = c(0,0), date_breaks = "1 month", date_label = "%b", 
               limits = c(as_date("2024-07-10"), as_date("2025-08-05"))) +
  theme_classic()
p3
p1 / p2 / p3 + plot_layout(heights = c(1, 2, 1))

load("output/data/maestro_df_clean.RData")

p3 <- maestro_df_clean |> 
  select(dttm_local, matches(c("pco2_12", "chla_12", "do_sat_12"))) |> 
  pivot_longer(cols = -dttm_local, names_to = "var", values_to = "value") |> 
    ggplot(aes(dttm_local, value)) +
  geom_point() + 
  facet_grid(var ~ ., scales = "free") +
  scale_x_date(expand = c(0,0), date_breaks = "1 month", date_label = "%b", 
               limits = c(as_date("2024-07-10"), as_date("2025-08-05"))) +
  theme_classic()
p3
