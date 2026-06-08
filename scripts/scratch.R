library(tidyverse)
library(microeco)

load("output/data/mt_16s.RData")


mt_16s$tax_table$p[mt_16s$tax_table$p |> 
  str_which(".*(C|c)yano")]

cyanos <- mt_16s$tax_table |> 
  filter(p == "p__Cyanobacteria")

mt_16s$tax_table <- mt_16s$tax_table |> 
  filter(p == "p__Cyanobacteria")

mt_16s$sample_table <- mt_16s$sample_table |> 
  filter(date < "2025-03-30")

mt_16s$tidy_dataset()

mt_16s$cal_abund()

mt_cyanos <- clone(mt_16s)
trans_abund$new(dataset = mt_16s, taxrank = "s", ntaxa = 20, high_level = "g")$
  plot_bar(facet = "strat_season", ggnested = TRUE, xtext_angle = 45)
save(mt_cyanos, file = "output/data/mt_cyanos.RData")

trans_abund$new(dataset = mt_16s, taxrank = "g", ntaxa = 20)$
  plot_bar(facet = "strat_season", xtext_angle = 45)

unique(mt_16s$tax_table$g)
mt_16s$tax_table |> 
  rownames_to_column("otu") |> 
  filter(otu == "OTU2") |> 
  print()


load("output/data/mt_mag.RData")

sed_tax <- mt_mag$tax_table |> 
  slice(str_which(mt_mag$tax_table$g, "(S|s)ediminibacterium"))


# Load in PPS info
library(tidyverse)
pps <- readxl::read_xlsx("pps_samples.xlsx", na = "NA") |> 
  mutate(pumped = as.logical(pumped), 
extracted = as.logical(extracted)) |> 
  select(full_id, everything()) |> 
  mutate(across(everything(), str_squish)) 
  mutate(across(c(6:8, 10), dmy_hms)) |> 
  separate_wider_delim(cols = inst_recovered, delim = " ", 
names = c("dow", "month", "day", "time", "year"), too_few = "align_start") |> 
  mutate(inst_recovered = mdy_hms(paste(month, day, year, time)), .after = inst_deployed) |> 
  select(-c(dow, month, day, year, time))


todo <- pps |> 
  filter(pumped == TRUE) |> 
  filter_out(extracted == TRUE) 
