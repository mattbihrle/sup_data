# Explore UBA genome and kegg genes

# Load Libraries
library(microeco)
library(tidyverse)

load("output/data/mt_mag.RData")

mt_mag$tax_table <- mt_mag$tax_table |> 
  mutate(bin = str_remove(bin, "b__"))
kegg_rollup <- read_tsv("output/data/kegg_rollup_all.tsv") |> 
  separate_wider_delim(cols = gene, delim = ",", 
                       names = c("gene1", "gene2", "gene3", "gene4", "gene5", "gene6", "gene7", "gene8", "gene9"), 
                       too_few = "align_start") |> 
  mutate(across(starts_with("gene"), ~str_trim(.x, side = "both") |> str_to_lower())) |> 
  left_join(mt_mag$tax_table, by = "bin")

uba_genomes <- mt_mag$tax_table |> 
  filter(f == "f__UBA2100") |> 
  select(bin) |> 
  mutate(bin = str_remove(bin, "b__")) |> 
  as_vector()
gene <- "pete"
kegg_rollup |> 
  filter(if_any(starts_with("gene"), ~ .x %in% gene)) |>
  # filter(str_detect(func, ".*methylglutamate.*")) |>
  # filter(bin %in% uba_genomes) |>
  # filter(kegg_path == "ko00195") |>
  # arrange(bin) |> 
  # nrow()
view()

# Now I want to look at other abundant taxa

trans_abund$new(mt_mag, ntaxa = 50, taxrank = "o")$plot_bar()


mt_mag$tax_table |> 
  # filter(g == "g__Bog-950") |> 
  # filter(f == "f__LD1") |> 
  view()

chitin <- clone(mt_mag)

chitin$tax_table <- chitin$tax_table |> 
  filter(o == "o__Chitinophagales")
chitin$cal_abund()

trans_abund$new(chitin, ntaxa = 20, taxrank = "f")$plot_bar() +
  ggtitle("relative abundance of families within order 'Chitinophagales'")

ld1 <- clone(mt_mag)
ld1$tax_table <- ld1$tax_table |> 
  filter(f == "f__LD1")
ld1$cal_abund()
trans_abund$new(ld1, ntaxa = 20, taxrank = "bin")$
  plot_bar() +
  ggtitle("relative abundance of bins within genus 'Bog-950'")

mt_mag$cal_abund()
bog <- trans_abund$new(mt_mag, ntaxa = 100, taxrank = "bin")


bog_bins <- mt_mag$tax_table |> 
  filter(g == "g__Bog-950") |> 
  select(bin) |> 
  as_vector()
bog$data_abund <- bog$data_abund |> 
  filter(Taxonomy %in% bog_bins)
bog$plot_bar()


# Read in the two bins we care about

WM12_003 <- read_tsv("output/data/kegg_output/WM12_S7_003_clean_kegg.tsv")
WM04_004 <- read_tsv("output/data/kegg_output/WM04_S4_004_clean_kegg.tsv")

join <- full_join(WM12_003, WM04_004, by = "id") |> 
  mutate(both = ifelse(is.na(bin.x) | is.na(bin.y), FALSE, TRUE)) |> 
  arrange(both)
table(join$both)

sum(is.na(join$bin.x))
sum(is.na(join$bin.y))

join |> 
  filter(id == "K01602") |> 
  nrow()

join |> 
  filter(both == FALSE) |> 
  mutate(bin = "both_bins") |> 
  select(bin, id) |> 
  write_tsv(file = "output/data/kegg_output/bog_950_genes_not_shared.tsv")
