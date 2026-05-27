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
# t_network <- trans_network$new(dataset = mt_gene_mag,
#                                cor_method = "spearman",
#                                use_corr_p_adjust = FALSE)
t_safe_network <- clone(t_network)
# save(t_network, file = "output/data/t_network.RData")

# Start here 
load("output/data/t_network.RData")

#filtering correlations
#this step corresponds to `bac.cor[bac.cor < 0.7] = 0` in the original script
# Using optimized coefficient threshold based on RMT (Deng et al. 2012)
t_network$cal_network(COR_p_thres = 0.01, COR_optimization = TRUE)

t_network$cal_module(method = "cluster_louvain")

save(t_network, file = "output/data/t_network.RData")
load("output/data/t_network.RData")
# Open node table to look at abundances

# Create a summary table of the nodes and which module they are part of
node_df <- data.frame(module = vertex_attr(t_network$res_network, "module"),
                      gene = vertex_attr(t_network$res_network, "name"),
                      gene_rel_abund = vertex_attr(t_network$res_network, "RelativeAbundance")
                      ) 
table(node_df$module)

# Pad the modules so that they are in numberical order

node_df$module <- str_remove(node_df$module, "M") |> 
  str_pad(width = 2, side = "left", pad = "0") |> 
  str_pad(width = 3, side = "left", pad = "M")

# Double check it looks alright
table(node_df$module)

# Then add the module data into the tax_table for later plotting

mt_gene_mag$tax_table <- mt_gene_mag$tax_table |> 
  left_join(node_df)
# view(modules)

# view(mt_gene_mag$tax_table)
# Remove the genes that weren't included in the modules

mt_gene_mag$tax_table <- mt_gene_mag$tax_table |> 
  drop_na(module)

nrow(mt_gene_mag$tax_table)
save(mt_gene_mag, file = "output/data/mt_gene_mag_mods.RData")
# Great the rows match first add rownames then tidy

mt_gene_mag$tidy_dataset()
# Extract taxa for each module

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
     vertex.label = NA,
     edge.width = 0.5,
     edge.color = "gray80")


#ok.... how about summing up the abundance of the modules and plotting them over time?
#i need to figure out the best way to do relative abundance - there's no clear way with cal_abund

otu_matrix <- mt_gene_mag$otu_table  # Taxa as rows, samples as columns

# FIRST SEE IF THESE AREN'T ALREADY RELATIVE ABUNDANCE
# --------------------------------------------------------------
# Convert to relative abundance (percentage)
rel_abund_otu <- apply(otu_matrix, 2, function(x) x / sum(x) * 100)

# Okay looks like they all sum to 100
rel_abund_otu |> 
  colSums()

# Now sum by module 
module_rel_abund <- data.frame(row.names = colnames(otu_matrix))

for(mod in unique(node_df$module)) {
  taxa_in_mod <- node_df$gene[node_df$module == mod]
  module_rel_abund[, mod] <- colSums(rel_abund_otu[taxa_in_mod, , drop = FALSE])
}

head(module_rel_abund)

module_rel_abund |> 
  rowSums()
# ---------------------------------------------------------------------

# Add time variable (replace 'Time' with your actual column name)
module_rel_abund$date <- mt_gene_mag$sample_table$date
module_rel_abund$strat_season <- mt_gene_mag$sample_table$strat_season
# module_rel_abund$Site <- mt_gene_mag$sample_table$Site

module_long <- pivot_longer(module_rel_abund,
                            cols = c(node_df$module),
                            names_to = "Module",
                            values_to = "relabund")
save(module_long, file = "output/data/module_long_genes.RData")
# Calculate abundances of modules
mt_gene_mag$cal_abund()


# Look at gene patterns across the modules

l1_palette <- setNames(
  paletteer::paletteer_d("ggthemes::Classic_20"), # alphabet has 26 colors
  c(unique(mt_gene_mag$tax_table$l1))
)

load("output/data/mt_gene_mag_mods.RData") 
mt_gene_mag$tax_table <- mt_gene_mag$tax_table|> 
  mutate(id = gene) |> 
  column_to_rownames("id")
# Look at specifically M1
mt_m1 <- clone(mt_gene_mag)

mt_m1$tax_table <- mt_m1$tax_table |> 
  filter(module == "M01")

mt_m1$tidy_dataset()$cal_abund()
# Plot abundances

rank = "l2"
ntaxa = 20
# Get vector of 20 most abundant families

trans_abund$new(dataset = mt_m1, taxrank = rank, ntaxa = ntaxa, high_level = "l1")$
  plot_bar(facet = "strat_season", ggnested = TRUE) +
  # scale_fill_manual(values = l1_palette, drop = TRUE) +
     ggtitle("genes", subtitle = paste("top", ntaxa, rank, "in M1")) +
  # theme(legend.position = "none") +
  NULL


# Look at specifically M2
mt_m2 <- clone(mt_gene_mag)

mt_m2$tax_table <- mt_m2$tax_table |> 
  filter(module == "M02")

mt_m2$tidy_dataset()$cal_abund()
# Plot abundances

# Get vector of 20 most abundant families

trans_abund$new(dataset = mt_m2, taxrank = rank, ntaxa = ntaxa, high_level = "l1")$
  plot_bar(facet = "strat_season", ggnested = TRUE) +
     ggtitle("gene", subtitle = paste("top", ntaxa, rank, "in M2")) +
  # scale_fill_manual(values = l1_palette, drop = TRUE) + 
  NULL

# Look at specifically M3
mt_m3 <- clone(mt_gene_mag)

mt_m3$tax_table <- mt_m3$tax_table |> 
  filter(module == "M03")

mt_m3$tidy_dataset()$cal_abund()
# Plot abundances

# Get vector of 20 most abundant families

trans_abund$new(dataset = mt_m3, taxrank = rank, ntaxa = ntaxa, high_level = "l1")$
  plot_bar(facet = "strat_season", ggnested = TRUE) +
     ggtitle("16S", subtitle = paste("top", ntaxa, rank, "in M3")) +
  # scale_fill_manual(values = l1_palette, drop = TRUE) +
NULL

# Look at specifically M4
mt_m4 <- clone(mt_gene_mag)

mt_m4$tax_table <- mt_m4$tax_table |> 
  filter(module == "M04")

mt_m4$tidy_dataset()$cal_abund()
# Plot abundances
# Get vector of 20 most abundant families

trans_abund$new(dataset = mt_m4, taxrank = rank, ntaxa = ntaxa, high_level = "l1")$
  plot_bar(facet = "strat_season", ggnested = TRUE) +
     ggtitle("gene", subtitle = paste("top", ntaxa, rank, "in M4")) +
  # scale_fill_manual(values = l1_palette, drop = TRUE) + 
  NULL


# See notes.qmd for interpretations, now checking to test eigen values

t_network$cal_eigen()

# Try to create a correlation heat map-----------------------------------------------------------------------

t1 <- trans_env$new(dataset = mt_gene_mag, env_cols = c(1:22), standardize = F, complete_na = T)

#Check to see what values are missing and such
sum_stats <- skimr::skim(t1$data_env)
sum_stats

# Remove any variables that have NAs in them as per skimr
na_vars <- sum_stats$skim_variable[which(sum_stats$complete_rate < 1)]

t1$data_env <- t1$data_env |> 
  dplyr::select(-any_of(na_vars))

colnames(t1$data_env)

#Reskim

sum_stats <- skimr::skim(t1$data_env)
sum_stats

# Remove other variables I don't want
t1$data_env <- t1$data_env |> 
  dplyr::select(-matches(c("depth", "pres", "use_maestro", "date")))
colnames(t1$data_env)


t1$cal_cor(add_abund_table = t_network$res_eigen)

modules <- c(sort(unique(t_network$res_node_table$module)))
modules
t1$plot_cor(text_y_order = modules[1:4])

# Looks like temp, cdom, and do are primarily a negative relationship except for M3 which is lightly upregulated. 
#-----------------------------------------------------------------------------
module_colors <- setNames(pals::alphabet(14), levels(module_long$Module))

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
    scale_color_manual(values = module_colors) +
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
        dplyr::select(date, depth, temp) |> 
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
p2 <- module_long |> 
  filter(Module %in% c("M01", "M02", "M03", "M04")) |> 
  ggplot(aes(x = date, y = relabund, color = Module, group = Module)) +
  # geom_rect(aes(xmin = sum_start, xmax = fall_start, ymax = max(relabund) +4 , ymin = max(relabund)+1),  fill = "grey", color = NULL) +
  # geom_rect(aes(xmin = fall_start, xmax = winter_start, ymax = max(relabund) +4 , ymin = max(relabund)+1),  fill = "green4",  color = NULL) +
  # geom_rect(aes(xmin = winter_start, xmax = spring_start, ymax = max(relabund) +4 , ymin = max(relabund)+1),  fill = "grey",  color = NULL) +
  # geom_rect(aes(xmin = spring_start, xmax = max(date), ymax = max(relabund) +4 , ymin = max(relabund)+1),  fill = "green4",  color = NULL) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "Module Abundance Over Time",
       y = "Relative Abundance") +
  scale_color_manual(values = module_colors, drop = TRUE) +
  theme_bw() +
  scale_x_date(expand = c(0,0), date_breaks = "1 month", date_label = "%b", limits = c(as_date("2024-07-10"), as_date("2025-08-05")))

p2

# Environmental data
load("output/data/maestro_df_clean.RData")
p3 <- maestro_df_clean |> 
  select(dttm_local, matches(c("cdom_12"))) |> 
  pivot_longer(cols = -dttm_local, names_to = "var", values_to = "value") |> 
    ggplot(aes(dttm_local, value)) +
  geom_point() + 
  facet_grid(var ~ ., scales = "free") +
  scale_x_date(expand = c(0,0), date_breaks = "1 month", date_label = "%b", 
               limits = c(as_date("2024-07-10"), as_date("2025-08-05"))) +
  theme_classic()
p3

p1/p2 + plot_layout(heights = c(1, 2))
p1 / p2 / p3 + plot_layout(heights = c(1, 2, 1))

test <- mt_gene_mag$tax_table |> 
  filter(module == "M04")
