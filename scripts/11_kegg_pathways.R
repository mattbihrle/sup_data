# Install libraries that I may need

# BiocManager::install("pathview")

# BiocManager::install("KEGGREST")
BiocManager::install("gage")
#Load Libraries
library(tidyverse)
library(pathview)
library(KEGGREST)
library(gage)


# Read in the foam rollup data

bin_names <- list.files("msi_downloads/rollup/foam/") |> 
  str_extract("WM.*[0-9]{6}")
  # str_remove("MA.*_000")

length(bin_names)

# # Import core bin list
# core_bins <- read_tsv("msi_downloads/Derep_list.tsv", col_names = F) |> 
#   mutate(bin = str_remove(X1, "^all_bins/")) |> 
#   mutate(bin = str_remove(bin, ".fa$")) |> 
#   mutate(X1 = NULL)

# core_bins <- as.character(core_bins$bin)

# core_bins <- core_bins |> str_remove("MAGScoT_cleanbin_000")

for(i in 1:length(bin_names)) {

  bin_name <- bin_names[i]

  bin_name
  #Set file path here
  path <- file.path("msi_downloads", "rollup", "foam")
  path
  # Select the file you want with the regex here
int_df <- list.files(path, pattern = paste0(".*", bin_name, ".*", ".tsv"), full.names = T) |> 
  read_tsv(show_col_types = F) |> 
  mutate(bin = bin_name, .before = everything()) |> 
      rename_with(tolower)
  
    if(i == 1){
    gene_df <- int_df
  } else {
    gene_df <- bind_rows(gene_df, int_df)
  }
}
foam_rollup_df <- gene_df

foam_rollup_df$bin <- foam_rollup_df$bin |> 
  str_remove(pattern = "MAGScoT_cleanbin_000")

# Load in the mag_data

load("output/data/mt_mag.RData")

core_bins <- mt_mag$tax_table$bin |> 
  str_remove("b__")

head(core_bins)
head(foam_rollup_df$bin)
# Create a new column to filter out bins that I have abundance data for
foam_rollup_df <- foam_rollup_df |> 
  mutate(core_bin = foam_rollup_df$bin %in% core_bins)

# Now start following tutorial from pathview

foam_int <- foam_rollup_df |> 
  filter(core_bin == "TRUE") |> 
  group_by(bin, id) |> 
  tally(name = "count")

keggFind(database = "pathway", query = "Nitrogen Cycle")

n_cycle_id <- 01310

foam_path <- foam_int |> 
  select(bin, id, count) |> 
  pivot_wider(names_from = "bin", 
  values_from = "count", 
values_fill = NA) |> 
  column_to_rownames("id") |> 
  as.matrix()

pv_bin_0 <- pathview(
  gene.data = foam_path[, "WM01_S1_001"],
  pathway.id = n_cycle_id,
  species = "ko",
  out.suffix = "pv_bin_0"
)

mt_chitin <- clone(mt_mag)

mt_chitin$tax_table <- mt_chitin$tax_table |> 
  slice(str_which(mt_chitin$tax_table$f, ".*Chitin"))

mt_chitin$tidy_dataset()$cal_abund()
plot <- trans_abund$new(mt_chitin, taxrank = "bin", ntaxa = "30", high_level = "g")$
  plot_bar(ggnested = TRUE, xtext_angle = 30, facet = "strat_season")
plot
plotly::plotly_build(plot)

# Try something different, see if I can plot the abundance of gene pathways-------------------------------------------------


g1 <- clone(mt_gene_mag)

g1$tax_table <- g1$tax_table |> 
  mutate(l1_num = as.numeric(str_extract(g1$tax_table$l1, "^.{2}")), 
         l2 = ifelse(is.na(l2), "none", l2), 
         l3 = ifelse(is.na(l3), "none", l3))

g1$sample_table <- g1$sample_table |> 
  filter(date < "2025-03-30")

g1$tidy_dataset()

g1$cal_betadiv()

  tb2 <- trans_beta$new(g1, group = "strat_season_2", measure = "bray")
tb2$cal_ordination(method = "NMDS")
tb2$cal_manova(manova_all = TRUE)
print(tb2$res_manova)

test <- tb2$plot_ordination(plot_color = "strat_season_2", plot_shape = "strat_season", plot_type = c("point", "ellipse"), ellipse_level = 0.95) +
  #   geom_text(
  #   data = tb2$dataset$otu_table, # Filter data for the label layer
  #   aes(label = colnames(tb2$dataset$otu_table)),
  #   vjust = -1.1, # Adjust vertical position
  #   hjust = 0.5,  # Adjust horizontal position
  #   color = "black"
  # ) +
  ggarrow::geom_arrow_chain(colour = "black") + 
    ggtitle(paste0("bray-  MAGs PCoA")) +
  # scale_color_viridis_b() +
  NULL
test

g1$sample_table$strat_season_2

options <- data.frame(table(mt_gene_mag$tax_table$l1))
options$Var1

for(i in 1:length(options$Var1)){

  vars <- unique(g1$tax_table$l1_num) |> 
    sort()

  num <- vars[i]
g2 <- clone(g1)

g2$tax_table <- g2$tax_table |> 
  # filter_out(l1_num == num)
  filter(l1_num == num)

  if(nrow(g2$tax_table) == 0) next

g2$tidy_dataset()

g2$cal_abund()

  if(any(g2$tax_table$l3 == "none")){
    taxrank = "l2"
    ntaxa = "20"
    plot <- trans_abund$new(g2, taxrank = taxrank, ntaxa = ntaxa)$
      plot_bar(facet = "strat_season") +
    # theme(legend.position = "none") + 
    ggtitle(paste(options$Var1[i], "ntaxa =", ntaxa, "rank = ", taxrank)) +
    NULL
  } else {
taxrank = "l3"
ntaxa = "20"
plot <- trans_abund$new(g2, taxrank = taxrank, ntaxa = ntaxa, high_level = "l2")$
    plot_bar(facet = "strat_season", ggnested = TRUE) +
  # theme(legend.position = "none") + 
  ggtitle(paste(options$Var1[i], "ntaxa =", ntaxa, "rank = ", taxrank)) +
  NULL
  }
print(plot)
  }


# Okay after looking at those graphs I want to remove some l1s and look again at the different abundances through the year.

g1$tax_table <- g1$tax_table |> 
  filter_out(l1_num %in% c(7, 6, 21, 19))

g1$tidy_dataset()$cal_abund()
# Check to see if it worled
unique(g1$tax_table$l1) |> 
  sort()

# Okay now that it worked, go ahead and plot

taxrank = "l2"
ntaxa = "20"
plot <- trans_abund$new(g1, taxrank = taxrank, ntaxa = ntaxa, high_level = "l1")$
    plot_bar(facet = "strat_season", ggnested = TRUE) +
  # theme(legend.position = "none") + 
  ggtitle(paste("l1 modified, ntaxa =", ntaxa, "rank = ", taxrank)) +
  NULL
plot


plot <- trans_abund$new(g1, taxrank = taxrank, ntaxa = ntaxa)$
  plot_bar(facet = "strat_season") +
    # theme(legend.position = "none") + 
  ggtitle(paste("l1 modified ntaxa =", ntaxa, "rank = ", taxrank)) +
    NULL
plot
