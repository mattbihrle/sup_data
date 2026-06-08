#' ---
#' title: 12_gene_cliques_sample
#' format: html
#' ---
#' 


library(microeco)

library(tidyverse)


load("output/data/mt_gene_sample.RData")


# Run this line only if you need to recreate the gene object, it will take a long time. 
t1_sample_ancomb <- trans_diff$new(dataset = mt_gene_sample, method = "ancombc2", group = "strat_season", fix_formula = "strat_season", alpha = 0.05, taxa_level = 'gene')
save(t1_sample_amcomb, file = "output/data/t1_sample_ancomb.RData")


view(t1$res_diff)

# t2 <- clone(t1)

# t2$res_diff$Factors <- t2$res_diff$Factors |> 
#     str_replace("\\((i|I)ntercept\\)", "summer") |> 
#     str_replace("strat_season.L", "fall") |> 
#     str_replace("strat_season.Q", "winter") |> 
#     str_replace("strat_season.C", "spring")

# sig_taxa <- t2$res_diff |> 
#     filter (P.adj <=0.05) |> 
#     dplyr::select(Taxa, Factors, Significance, lfc ) |> 
#     separate_wider_delim(Taxa, delim = "|", names = c("bin", "l1", "l2", "l3", "l4", "func", "id")) |> 
#     arrange(Factors, -abs(lfc), Significance)

# view(sig_taxa)

# t3 <- clone(t1)

# t3$res_diff$Factors <- t3$res_diff$Factors |> 
#     str_replace("\\((i|I)ntercept\\)", "summer") |> 
#     str_replace("strat_season.L", "fall") |> 
#     str_replace("strat_season.Q", "winter") |> 
#     str_replace("strat_season.C", "spring")

# t3$res_diff <- t3$res_diff |>  
#     dplyr::select(Taxa, Factors, Significance, lfc, P.adj ) |>
#     filter(P.adj <= 0.05) |> 
#     separate_wider_delim(Taxa, delim = "|", names = c("bin", "l1", "l2", "l3", "l4", "func", "id"))

# t3$plot_diff_bar(keep_full_name = TRUE, heatmap_cell = "P.adj", heatmap_sig = "Significance", heatmap_x = "Factors", heatmap_y = "id")

# fig_pair <- t2$res_diff |> 
#     filter(t2$res_diff$P.adj < 0.05) |> 
#     mutate(Taxa = str_remove(Taxa, "\\|WM.*")) |> 
#     mutate(Taxa = str_remove(Taxa, "d__Bacteria\\|")) |> 
#     mutate(Taxa = str_remove(Taxa, "\\|s__.*")) |> 
#     distinct(Taxa, .keep_all = T) |> 
#     ggplot(aes(x = Factors, y = Taxa, fill = lfc)) + 
#     geom_tile(color = "black") +
#     scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
#                          na.value = "white",
#                          name = NULL) +
#     geom_text(aes(Factors, Taxa, label = lfc), size = 5, color = "black") +
#     scale_color_identity(guide = FALSE) +
#     labs(x = NULL, y = NULL, title = "Gene Log fold changes compared to summer, (p < 0.1)") +
#     theme_minimal() +
#     theme(plot.title = element_text(hjust = 0.5), text = element_text(color = "black", size = 15))
# fig_pair


# Run this line only if you need to recreate the gene object, it will take a long time. 
t1_sample_alde <- trans_diff$new(dataset = mt_gene_sample, method = "ALDEx2_kw", group = "strat_season", fix_formula = "strat_season", alpha = 0.05, taxa_level = 'gene')
save(t1_sample_alde, file = "output/data/t1_sample_alde.RData")


# Load in all the crazy libraries

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
library(vegan) # Multivariate ecological analysis
library(propr)
library(missForest) # Imputing missing values in dataframes using Random Forests
library(VSURF) # Random Forests approach to variable importance identification
library(car) #for scatterplot
library(patchwork)



samp_genes <- clone(mt_gene_sample)

#Filter taxa with a certain relative abundane threshold
samp_genes$filter_taxa(
  rel_abund = 0.00001
)

# Retaining edges with Spearman correlation > 0.65
samp_network <- trans_network$new(dataset = samp_genes,
                               cor_method = "spearman",
                               use_corr_p_adjust = FALSE)
samp_safe_network <- clone(samp_network)
save(samp_network, file = "output/data/samp_network.RData")



#filtering correlations
#this step corresponds to `bac.cor[bac.cor < 0.7] = 0` in the original script
# Using optimized coefficient threshold based on RMT (Deng et al. 2012)
samp_network$cal_network(COR_p_thres = 0.01, COR_optimization = TRUE)

samp_network$cal_module(method = "cluster_louvain")

save(samp_network, file = "output/data/samp_network.RData")



# Create a summary table of the nodes and which module they are part of
node_df <- data.frame(module = vertex_attr(samp_network$res_network, "module"),
                      gene = vertex_attr(samp_network$res_network, "name"),
                      gene_rel_abund = vertex_attr(samp_network$res_network, "RelativeAbundance")
                      ) 
table(node_df$module)

# Pad the modules so that they are in numberical order

node_df$module <- str_remove(node_df$module, "M") |> 
  str_pad(width = 2, side = "left", pad = "0") |> 
  str_pad(width = 3, side = "left", pad = "M")

# Double check it looks alright
table(node_df$module)

# Then add the module data into the tax_table for later plotting

mt_gene_sample$tax_table <- mt_gene_sample$tax_table |> 
  left_join(node_df)
# view(modules)

# view(mt_gene_sample$tax_table)
# Remove the genes that weren't included in the modules

mt_gene_sample$tax_table <- mt_gene_sample$tax_table |> 
  drop_na(module)

nrow(mt_gene_sample$tax_table)
save(mt_gene_sample, file = "output/data/mt_gene_sample_mods.RData")
# Great the rows match first add rownames then tidy

mt_gene_sample$tidy_dataset()
# Extract taxa for each module

layout <- layout_with_fr(samp_network$res_network)

module_colors <- c("M1" = "red",
                   "M2" = "orange",
                   "M3" = "green",
                   "M4" = "blue")

vertex_colors <- module_colors[modules]
# Maybe add sample names here?
plot(samp_network$res_network,
     layout = layout,
     vertex.size = 5,
     vertex.label = NA,
     edge.width = 0.5,
     edge.color = "gray80")


#ok.... how about summing up the abundance of the modules and plotting them over time?
#i need to figure out the best way to do relative abundance - there's no clear way with cal_abund

otu_matrix <- mt_gene_sample$otu_table  # Taxa as rows, samples as columns

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
module_rel_abund$date <- mt_gene_sample$sample_table$date
module_rel_abund$strat_season <- mt_gene_sample$sample_table$strat_season
# module_rel_abund$Site <- mt_gene_sample$sample_table$Site

module_long <- pivot_longer(module_rel_abund,
                            cols = c(node_df$module),
                            names_to = "Module",
                            values_to = "relabund")
save(module_long, file = "output/data/module_long_genes_sample.RData")
# Calculate abundances of modules
mt_gene_sample$cal_abund()


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


## Test plotting cliques over the temperature
load("output/data/temp_df_long_full.RData")

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


