library(microeco)
library(file2meco)
 
p2 <- ggplot(data, aes(x = Abundance, y = Growth_Rate, color = Species)) +
  geom_point(size = 4) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "black") +
  theme_minimal() +
  labs(title = "Growth Rate (r) vs. Abundance", x = "Abundance", y = "Growth Rate (r)")

# Print the plots
print(p1)
print(p2)

rollup <- read_tsv("msi_downloads/metacerb/all_files/uba_genomes/MetaCerb_WM15_S10_MAGScoT_cleanbin_000035/final/rollup/rollup_prodigal_WM15_S10_MAGScoT_cleanbin_000035-KOFam_all_KEGG.tsv") |> 
  arrange(ID) |> 
  select(ID) |> 
  distinct()

summary <- read_tsv("msi_downloads/metacerb/all_files/uba_genomes/MetaCerb_WM15_S10_MAGScoT_cleanbin_000035/final/annotations-prodigal_WM15_S10_MAGScoT_cleanbin_000035/final_annotation_summary.tsv") |> 
  drop_na(best_hit) |> 
  arrange(best_hit) |>
  select(best_hit) |> 
  distinct()
waldo::compare(rollup$ID, summary$best_hit)

kegg <- read_tsv("output/data/kegg_output/WM15_S10_035_kegg.tsv") |> 
  arrange(id) |> 
  distinct() |> 
  mutate(kegg = TRUE)

rollup <- kegg_table |> 
  filter(bin == "WM15_S10_035") |> 
  select(-n) |> 
  distinct() |> 
  arrange(id) |> 
  mutate(rollup = TRUE)
sum((rollup$id %in% kegg$id))

sum(kegg$id %in% rollup$id)

new <- rollup |> 
  full_join(kegg, by = c("id", "bin"))


library(tidyverse)
foam <- read_tsv("msi_downloads/metacerb/all_files/uba_genomes/MetaCerb_WM15_S10_MAGScoT_cleanbin_000035/final/annotations-prodigal_WM15_S10_MAGScoT_cleanbin_000035/annotation_summary_KOFam_all_FOAM.tsv") |> 
  drop_na(best_hit)
kegg <- read_tsv("msi_downloads/metacerb/all_files/uba_genomes/MetaCerb_WM15_S10_MAGScoT_cleanbin_000035/final/annotations-prodigal_WM15_S10_MAGScoT_cleanbin_000035/annotation_summary_KOFam_all_KEGG.tsv") |> 
  drop_na(best_hit)
sum <- read_tsv("msi_downloads/metacerb/all_files/uba_genomes/MetaCerb_WM15_S10_MAGScoT_cleanbin_000035/final/annotations-prodigal_WM15_S10_MAGScoT_cleanbin_000035/final_annotation_summary.tsv") |> 
  drop_na(best_hit) 

foam$best_hit |> 
  unique() |> 
  length()

kegg$best_hit |> 
  unique() |> 
  length()

sum$best_hit |> 
  unique() |> 
  length() %in%
  kegg$best_hit

sum(sum$best_hit %in% kegg$best_hit)
`%!in%` = Negate(`%in%`)
smol <- kegg |> 
  filter(kegg$best_hit %!in% sum$best_hit)

#looking at these it seems like summary KEGG is the way to go, try it again with 
# another bin just to see



foam <- read_tsv("msi_downloads/metacerb/all_files/uba_genomes/MetaCerb_WM12_S7_MAGScoT_cleanbin_000038/final/annotations-prodigal_WM12_S7_MAGScoT_cleanbin_000038/annotation_summary_KOFam_all_FOAM.tsv") |> 
  drop_na(best_hit)
kegg <- read_tsv("msi_downloads/metacerb/all_files/uba_genomes/MetaCerb_WM12_S7_MAGScoT_cleanbin_000038/final/annotations-prodigal_WM12_S7_MAGScoT_cleanbin_000038/annotation_summary_KOFam_all_KEGG.tsv") |> 
  drop_na(best_hit)
sum <- read_tsv("msi_downloads/metacerb/all_files/uba_genomes/MetaCerb_WM12_S7_MAGScoT_cleanbin_000038/final/annotations-prodigal_WM12_S7_MAGScoT_cleanbin_000038/final_annotation_summary.tsv") |> 
  drop_na(best_hit) 

foam$best_hit |> 
  unique() |> 
  length()

kegg$best_hit |> 
  unique() |> 
  length()

sum$best_hit |> 
  unique() |> 
  length()

sum(sum$best_hit %in% kegg$best_hit)
`%!in%` = Negate(`%in%`)
smol <- kegg |> 
  filter(kegg$best_hit %!in% sum$best_hit)


# Okay, one more check to see if I want to use rollups or kegg summary

bins <- list.files("msi_downloads/metacerb/all_files/uba_genomes/", full.names = T) |> 
  str_extract("WM.*_.*_.*_.*_[0-9]{6}") |> 
  na.omit()
 i <- 5
final_filename <- paste0("msi_downloads/metacerb/all_files/uba_genomes/MetaCerb_", bins[i], "/final/annotations-prodigal_", bins[i], "/final_annotation_summary.tsv")
rollup_filename <- paste0("msi_downloads/metacerb/all_files/uba_genomes/MetaCerb_", bins[i], "/final/rollup/rollup_prodigal_", bins[i], "-KOFam_all_KEGG.tsv")
summary_filename <- paste0("msi_downloads/metacerb/all_files/uba_genomes/MetaCerb_", bins[i], "/final/annotations-prodigal_", bins[i], "/final_annotation_summary.tsv")
final <- read_tsv(final_filename) |> 
  drop_na(best_hit)
rollup <- read_tsv(rollup_filename) |> 
  distinct(ID, .keep_all = T)

# Look at which any genes that are in the final table but not in the rollup table
sum(!final$best_hit %in% rollup$ID)
# Look at which genes are in the rollup table but not in the final table
sum(!rollup$ID %in% final$best_hit)
missing_from_final <- rollup |> 
  filter_out(rollup$ID %in% final$best_hit) |> 
  view()

if(sum(missing_from_final$ID %in% final$best_hit) == 0){
  paste("Rollup table has genes that are not in the final kegg table for bin", bins[i])
}
if(sum(final$best_hit %in% rollup$ID) == nrow(final)){
  paste("All genes from the 'final kegg table are present in the rollup table for bin", bins[i])
} else {
  paste("Not all genes from the 'final' kegg table are present in the rollup table for bin", str_remove(bins[i], "MAGScoT_cleanbin_000"))
}

sum <- read_tsv(summary_filename) |> 
  drop_na(best_hit)

sum |> 
  filter(missing_from_final$ID %in% sum$best_hit) |>
  view()

# Here I want to see everything in summary that is in 
sum |> 
  filter_out(sum$best_hit %in% final$best_hit %in% rollup$ID) |> 
  view()


id <- paste0("msi_downloads/metacerb/all_files/uba_genomes/MetaCerb_", bins[i],
             "/step_10-visualizeData/prodigal_", bins[i], "/KOFam_all_KEGG_level-id.tsv") |> 
  read_tsv()

final$best_hit %in% id$Id

sum(!rollup$ID %in% id$Id) 

sum(!id$Id %in% rollup$ID)
sum(!id$Id %in% final$best_hit)

# Okay so confirming here that I do want to use the rollup files. 

mt_mag$tax_table <- mt_mag$tax_table |> 
  mutate(bin = str_remove(bin, "b__"))
kegg_rollup <- read_tsv("output/data/kegg_rollup_all.tsv") |> 
  separate_wider_delim(cols = gene, delim = ",", 
                       names = c("gene1", "gene2", "gene3", "gene4", "gene5", "gene6", "gene7", "gene8", "gene9"), 
                       too_few = "align_start") |> 
  mutate(across(starts_with("gene"), ~str_trim(.x, side = "both") |> str_to_lower())) |> 
  left_join(mt_mag$tax_table, by = "bin")

mag_tax


uba_genomes <- mt_mag$tax_table |> 
  filter(f == "f__UBA2100") |> 
  select(bin) |> 
  mutate(bin = str_remove(bin, "b__")) |> 
  as_vector()

genes <- c("h4mptp")
gene <- "ftf"
kegg_rollup |> 
  filter(if_any(starts_with("gene"), ~ .x %in% genes)) |>
  # filter(str_detect(func, ".*methylglutamate.*")) |>
  # filter(bin %in% uba_genomes) |>
  nrow()
  # view()
  # distinct(func) |> 
  # print(n = 100)

kegg_rollup |> 
  filter(bin %in% uba_genomes) |> 
  filter(l1 == "metabolism") |> 
  view()
  
mt_mag$tax_table |> 
  filter(bin == "b__WM04_S4_024") |> 
  view()
mt_16s$sample_table |> 
  view()


mass <- clone(mt_16s)
mass$tax_table <- mass$tax_table |> 
  filter(g == "g__Massilia")
mass$cal_abund()
  trans_abund$new(mass, taxrank = "otu", ntaxa = 30)$
  plot_bar()

load("output/data/mt_16s.RData")
otus <- mt_16s$tax_table |> 
  filter(g == "g__Patulibacter") |> 
  select(otu) |> 
  as_vector() |> 
  str_remove("o__")

a <- mt_16s$otu_table |> 
  filter(rownames(mt_16s$otu_table) %in% otus) |> 
  select(!where(~sum(.x) == 0)) |> 
  colnames()


rose <- trans_abund$new(test, ntaxa = 200, taxrank = "g")$
  plot_bar()

rose_plot <- rose$data |> 
  filter(Taxonomy == "Roseomonas") |> 
  ggplot(aes(Sample, Abundance)) +
  geom_col()
plotly::plotly_build(rose_plot)

mt_16s$sample_table |> 
  select(full_id, ext_date) |> 
  view()

rosey <- clone(mt_16s)

rosey$tax_table <- rosey$tax_table |> 
  filter(g == "g__Roseomonas") 
rosey$cal_abund()

trans_abund$new(rosey, ntaxa = 15, taxrank = "otu")$plot_bar() |> 
  plotly::plotly_build()
smol <- clone(mt_16s)
smol$tax_table <- smol$tax_table |> 
  filter_out(g == "g__unknown")
trans_abund$new(smol, ntaxa = 1000, taxrank = "otu", high_level = "g")$
  plot_bar(ggnested = T) |> 
  plotly::plotly_build()


sed <- clone(mt_16s)

sed$tax_table <- sed$tax_table |> 
  filter(g == "g__Sphingomonas")
sed$cal_abund()
trans_abund$new(sed, ntaxa = 100, taxrank = "otu", high_level = "s")$
  plot_bar(ggnested = TRUE) |> 
  plotly::plotly_build()


full_t_df |> 
  filter(str_detect(xml_file, "Superior")) |> 
  filter(dttm_utc > "2024-08-15" & dttm_utc < "2025-09-01") |> 
  view()
load("output/data/mt_16s.RData")
load("output/data/sample_table.RData")

mt_new <- microtable$new(tax_table = mt_16s$tax_table, otu_table = mt_16s$otu_table, sample_table = sample_table)

mt_new$tidy_dataset()

mt_16s$tax_table |> 
  filter(str_detect(g, "Cyano")) |> 
  view()




# Function to test a bunch of different normalizations and plotting

plot_nmds <- function(mt_16s, norm = NULL, betadiv = NULL, group = "strat_season") {
# First normalize
    for(i in 1:length(norm)) {
  n <- trans_norm$new(mt_16s)
n <- n$norm(method = norm[i])
n$cal_betadiv(method = betadiv)
 # Then run the normalization with each of the beta diversity measures
 for(ii in 1:length(betadiv)){
tb2 <- trans_beta$new(n, group = "solar_season", measure = betadiv[ii])
tb2$cal_ordination(method = "NMDS")

plot <- tb2$plot_ordination(plot_color = group, plot_shape = "deployment", plot_type = c("point", "ellipse"), NMDS_stress_pos = NULL, ellipse_level = 0.95) +
  # geom_arrow_chain(colour = "black") + 
    ggtitle("NMDS Ordination of WM Data", subtitle = paste("Stress = ", round(tb2$res_ordination$model$stress, 2), 
  "Normalization:", norm[i], "Distance:", betadiv[ii])) + 
    # geom_label_repel(label = mt_16s$sample_table$full_id, color = "black") +
    labs(color = group) +
    guides(fill = "none") +
    theme_classic() +
    theme(
      text = element_text(size = 15),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 3)
    )
  print(plot)
    }
  }
}
norms <- c("hellinger", "rarefy", "SRS", "GMPR", "CSS", "DESeq2", "normalize", "standardize", "log")
betadiv = c("bray", "robust.aitchison")
plot_nmds(mt_16s, norm = norms, betadiv = betadiv, group = "strat_season")

library(mecoturn)
load("output/data/mt_16s.RData")
b <- betaturn$new(mt_16s, measure = "betaNRI", filter_thres = 0.0001)
save(b, "output/data/betaturn.RData")
b$cal_group_distance(group = "solar_season")
head(b$res_group_distance)
data(wheat_16S)
wheat_16S$sample_table

#Convert .nwk tree file to phylo format

tree <- ape::read.tree("data/16s/lotus3_out/OTUphylo.nwk")

mt_16s$phylo_tree <- tree
mt_16s
mt_16s$tidy_dataset()
mt_16s


# Functions I'd use?

# Load and install packages

load("output/data/mt_16s.RData")

s_mt <- clone(mt_16s)

for(i in 1:12){
  s_mt <- clone(mt_16s)
  s_mt$sample_table <- s_mt$sample_table |> 
  filter(lake == "superior") |> 
  mutate(month = month(date)) |> 
    filter(month == i)
  s_mt$tidy_dataset()

  a <- trans_abund$new(s_mt, taxrank = "g", ntaxa = 50)$plot_bar(facet = "deployment") +
  ggtitle(paste0("Month:", i))
  print(a)
}
load("output/data/mt_16s.RData")
mt_r <- clone(mt_16s) 

mt_r$rarefy_samples()

n <- trans_norm$new(mt_r)

n <- n$norm(method = "hellinger")


n$cal_betadiv(method = "bray")

tb2 <- trans_beta$new(n, group = "strat_season", measure = "bray")
tb2$cal_ordination(method = "NMDS")

plot <- tb2$plot_ordination(plot_color = "strat_season", plot_shape = "deployment", plot_type = c("point", "ellipse"), NMDS_stress_pos = NULL, ellipse_level = 0.95) +
  # geom_arrow_chain(colour = "black") + 
    ggtitle("NMDS Ordination of WM Data", subtitle = paste("Stress = ", round(tb2$res_ordination$model$stress, 2), 
  "Normalization:", "hellinger", "Distance:", "bray")) + 
    # geom_label_repel(label = mt_16s$sample_table$full_id, color = "black") +
    labs(color = "strat_season") +
    guides(fill = "none") +
    theme_classic() +
    theme(
      text = element_text(size = 15),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 3)
    )
  plot


mt_16s$sample_table |> 
  filter(deployment == "s_2122_WM") |> 
  select(full_id, strat_season, date, everything()) |> 
  view()

