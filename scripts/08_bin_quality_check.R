library(tidyverse)

all_bins_qual <- read_tsv("msi_downloads/all_bins_quality_report.tsv") |> 
  mutate(quality = ifelse(Completeness > 90 & Contamination < 5, "HIGH", "na")) |> 
  mutate(quality = ifelse(Completeness >= 50 & Contamination < 10 & quality != "HIGH", "MED", quality)) |>
  mutate(quality = ifelse(Completeness < 50 & Contamination < 10, "LOW", quality)) |> 
  mutate(quality = ifelse(Contamination >= 10, "contamination_flag", quality)) |> 
  select(Name, Completeness, Contamination, quality, everything())

# Pull out which are high quality

high_qual <- all_bins_qual |> 
  filter(quality == "HIGH")

med_qual <- all_bins_qual |> 
  filter(quality == "MED")

low_qual <- all_bins_qual |> 
  filter(quality == "LOW")

contam <- all_bins_qual |> 
  filter(quality == "contamination_flag")

bins_abundance <- read_tsv("msi_downloads/Drep_Bins_coverage.tsv") |> 
  slice(-1)

bins_derep <- read_tsv("msi_downloads/Derep_list.tsv", col_names = F) |> 
  mutate(bin = str_extract(X1, "WM.*000[0-9]{3}"))

all_bins_qual <- all_bins_qual |> 
  mutate(bin_abundance =  as.character(all_bins_qual$Name %in% bins_abundance$Genome), 
  bin_derep = all_bins_qual$Name %in% bins_derep$bin) |> 
  select(Name, bin_abundance, bin_derep, quality, everything())


bins_derep_2 <- read_tsv("msi_downloads/Derep_list_04_23.tsv", col_names = F) |> 
  mutate(bin = str_extract(X1, "WM.*000[0-9]{3}"))

bins_abundance_2 <- read_tsv("msi_downloads/Drep_Bins_coverage_04_23.tsv") |> 
  slice(-1)


  all_bins_qual <- all_bins_qual |> 
  mutate(bin_abundance_2 =  as.character(all_bins_qual$Name %in% bins_abundance_2$Genome), 
  bin_derep_2 = all_bins_qual$Name %in% bins_derep_2$bin) |> 
  select(Name, bin_abundance, bin_abundance_2, bin_derep, bin_derep_2, quality,  everything())


waldo::compare(all_bins_qual$bin_abundance, all_bins_qual$bin_abundance_2)

waldo::compare(all_bins_qual$bin_derep, all_bins_qual$bin_derep_2)

bins_abundance_2 <- read_tsv("msi_downloads/Drep_Bins_coverage_04_23.tsv")
bins_abundance <- read_tsv("msi_downloads/Drep_Bins_coverage.tsv") |> 
  select(!matches("Mean"))
colnames(bins_abundance) <- str_extract(colnames(bins_abundance), "WM.{2}_.*_R.{1}")


prefixes <- sprintf("WM%02d", 1:24)

df_final <- bind_cols(
  bins_abundance,
  map_dfc(prefixes, function(prefix) {
    bins_abundance %>%
      # Select columns that start with the current prefix
      select(starts_with(prefix)) %>%
      # Calculate row means and dynamically name the new column (e.g., WM01_mean)
      transmute(!!paste0(prefix, "_mean") := rowMeans(., na.rm = TRUE))
  })
)

df_final <- df_final |> 
  select(matches("_mean")) |> 
  # select(!WM06_mean:WM10_mean) |> 
  rename_with( ~str_extract(.x, pattern = "WM.{2}"))

bins_abundance_2 <- bins_abundance_2 |> 
  select(matches("Relative")) |> 
  rename_with( ~str_extract(.x, pattern = "WM.{2}"))


waldo::compare(bins_abundance_2, df_final)
