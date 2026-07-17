library(microeco)
library(file2meco)
 |> |> |> |> |> |> |> |> |> |> |> |> |> |> |> |> |> |> |> |> |> |> |> |> |> ion and its growth rate (r).
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
  length()

# Starting from here, I should see if I can just remove the duplicates. 