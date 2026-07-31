# Load data
load("output/data/mt_16s.RData")


# Load Libraries 
packages <- c("tidyverse", "microeco")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
# Packages loading
invisible(lapply(packages, library, character.only = TRUE))
# Remove the list of packages and installed packages from the environment
rm(installed_packages, packages)

# -------------------------------------------------------------------------

clone_16s <- clone(mt_16s)
clone_16s$sample_table <- clone_16s$sample_table |> 
  filter(lake == "superior")
clone_16s$filter_taxa(
  rel_abund = 0.0001,  # At least 0.001% in any sample
)
clone_16s$tidy_dataset()
clone_16s

# t1 <- trans_diff$new(dataset = clone_16s, method = "ALDEx2_kw", group = "strat_season", alpha = 0.05, taxa_level = 'otu')
# save(t1, file = "output/data/16s_alde.RData")
t2 <- trans_diff$new(dataset = clone_16s, method = "ancombc2", group = "strat_season", alpha = 0.05, taxa_level = 'otu')
save(t2, file = "output/data/16s_ancombc2.RData")
print("All done!")
load("output/data/16s_ancombc2.RData")
load("output/data/16s_alde.RData")
sig_taxa <- t2$res_diff |> 
  filter(P.adj < 0.05) |> 
  rename(p_val = p) |> 
  separate_wider_delim(Taxa, delim = "|", names = colnames(mt_16s$tax_table[1:8]), too_few = "align_start") |> 
  mutate(across(any_of(colnames(mt_16s$tax_table[1:8])), ~str_remove(.x, "^[a-z]{1}__"))) |> 
  filter(!grepl(".*Intercept.*", Factors)) |> 
  arrange(desc(diff_robust), desc(passed_ss), lfc)


sig_taxa_alde <- t1$res_diff |> 
  filter(P.adj < 0.05) |> 
  # rename(p_val = p) |> 
  separate_wider_delim(Taxa, delim = "|", names = colnames(mt_16s$tax_table[1:8]), too_few = "align_start") |> 
  mutate(across(any_of(colnames(mt_16s$tax_table[1:8])), ~str_remove(.x, "^[a-z]{1}__"))) |> 
  arrange(P.adj)

sig_taxa <- sig_taxa|> 
  mutate(shared = otu %in% sig_taxa_alde$otu) |> 
  arrange(desc(shared))

sig_taxa_alde <- sig_taxa_alde|> 
  mutate(shared = otu %in% sig_taxa$otu) |> 
  arrange(desc(shared))

# Here I can pull out the taxa that are significant to both

mt_16s

mt_16s$taxa_abund |> 
  head()


mt_16s$taxa_abund$otu |> 
  mutate(otu = str_extract(rownames(mt_16s$taxa_abund$otu), "(?<=\\o__)OTU[0-9]*")) |> 
  filter(otu %in% sig_taxa$otu[sig_taxa$shared == TRUE]) |> 
  pivot_longer(cols = -otu, names_to = "full_id", values_to = "abundance") |>
  left_join(mt_16s$sample_table, by = "full_id") |>
  # now plot
  ggplot(aes(x = full_id, y = abundance, fill = otu)) +
    geom_bar(position = "stack", stat = "identity") +
  facet_wrap(~solar_season, scales = "free_x")

mt_16s$tax_table |> 
  filter(otu %in% paste0("o__", sig_taxa$otu[sig_taxa$shared == TRUE])) |>
  view()
