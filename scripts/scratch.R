library(tidyverse)

mt_16s$tax_table |> 
  rownames_to_column("otu") |> 
  filter(otu == "OTU2") |> 
  print()


load("output/data/mt_mag.RData")

sed_tax <- mt_mag$tax_table |> 
  slice(str_which(mt_mag$tax_table$g, "(S|s)ediminibacterium"))


