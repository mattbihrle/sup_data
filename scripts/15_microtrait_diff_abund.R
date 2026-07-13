# Microtrait differential abundance (July 13, 2026)


trans_abund$new(mt_microtrait, ntaxa = 50, taxrank = "l7")$plot_bar()

t1 <- trans_diff$new(dataset = mt_microtrait, method = "lefse", 
                     group = "strat_season_2", 
                     alpha = 0.05, taxa_level = 'all', lefse_subgroup = NULL)

sig_taxa <- t1$res_diff |>
  separate_wider_delim(cols = Taxa, delim = "|", names = c("l1", "l2", "l3", "l4",
                                                           "l5", "l6", "l7", "l8"),
                       too_few = "align_start") |>
  mutate(level = 8 - rowSums(is.na(across(matches("l[0-9]{1}"))))) |>
  select(level, Significance, everything()) |>
  arrange(l1, l2, l3, l4, l5, l6, l7, l8)

View(sig_taxa)


t2 <- trans_diff$new(dataset = mt_microtrait, method = "ancombc2", 
                     fix_formula = "strat_season", alpha = 0.05)
View(t2$res_diff)

sig_taxa_ancombc2 <- t2$res_diff |> 
  separate_wider_delim(cols = Taxa, delim = "|", 
                       names = c("l1", "l2", "l3", "l4", "l5", "l6", "l7", "l8", "ASV"),
                       too_few = "align_start") |>
  mutate(level = 8 - rowSums(is.na(across(matches("l[0-9]{1}"))))) |>
  select(level, Significance, everything()) |> 
  filter(P.adj < 0.05) |> 
  arrange(Significance, level) |> 
  View()


t2 <- trans_diff$new(dataset = mt_microtrait, method = "ancombc2", 
                     fix_formula = "strat_season_2", alpha = 0.05)
View(t2$res_diff)

sig_taxa_ancombc2_sum_winter <- t2$res_diff |> 
  separate_wider_delim(cols = Taxa, delim = "|", 
                       names = c("l1", "l2", "l3", "l4", "l5", "l6", "l7", "l8", "ASV"),
                       too_few = "align_start") |>
  mutate(level = 8 - rowSums(is.na(across(matches("l[0-9]{1}"))))) |>
  select(level, Significance, everything()) |> 
  filter(P.adj < 0.05) |> 
  filter_out(Factors == "(Intercept)") |> 
  arrange(Significance, level)

sig_taxa_ancombc2_sum_winter

mt_mag_large$tax_table |> 
  select(d, p, c, o, f, g, s, bin, matches("g3.*desaturase"))


trans_abund$new(mt_mag_large, ntaxa = 50, taxrank = "g", 
                high_level = "g3_stress_tolerance__specific__temperature__low__membrane_stabilization__fatty_acid_desaturase")$
  plot_bar(ggnested = T) +
  ggtitle("top 50 families split by whether they have the fatty acid desaturase gene")

fad <- clone(mt_microtrait)


fad$tax_table <- fad$tax_table |> 
  filter(l1 == "stress tolerance")

fad$cal_abund()

fad$tax_table
plot <- trans_abund$new(fad, ntaxa = 50, taxrank = "l7")$plot_bar()

plotly::plotly_build(plot)


# Now looking for hydrogenotrophy 


fad <- clone(mt_microtrait)


fad$tax_table <- fad$tax_table |> 
  filter(l1 == "resource use")

fad$cal_abund()

fad$tax_table
plot <- trans_abund$new(fad, ntaxa = 50, taxrank = "l7")$plot_bar()

plotly::plotly_build(plot)