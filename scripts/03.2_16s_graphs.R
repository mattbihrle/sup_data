# 16S Abundances

# Load data
load("output/data/mt_16s.RData")
load("output/data/full_16s.RData")
load("output/data/contam_16s.RData")
list.files("output/data/", pattern = ".RData")
# Load Libraries 
packages <- c("plotly", "ggalluvial", "ggarrow", "tidyverse", "microeco", "ggrepel", "ggh4x")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
# Packages loading
invisible(lapply(packages, library, character.only = TRUE))
# Remove the list of packages and installed packages from the environment
rm(installed_packages, packages)

# Start analysis ------------------------------------------------
load("output/data/mt_16s.RData")
mt_16s

# First I want to remove the lake erie samples
mt_16s$sample_table <- mt_16s$sample_table |> 
  filter_out(lake == "erie")
mt_16s$tidy_dataset()
mt_16s

mt_16s$rarefy_samples(method = "rarefy")

mt_16s$cal_abund()
mt_16s$taxa_abund$p[1:5,1:5]

# Calc Alpha Diversity
mt_16s$cal_alphadiv()

# Calc Beta Diversity
mt_16s$cal_betadiv(method = c("bray","aitchison","robust.aitchison", "jaccard"))
mt_16s$beta_diversity$bray
# 16S, Robust Aitchison, NMDS and stressplot
n <- trans_norm$new(mt_16s)
n <- n$norm(method = "clr")
n$cal_betadiv(method = "robust.aitchison")
tb2 <- trans_beta$new(n, group = "solar_season", measure = "robust.aitchison")
tb2$cal_ordination(method = "NMDS")
tb2$cal_manova(manova_all = TRUE)
print(tb2$res_manova)

stress <- vegan::stressplot(tb2$res_ordination$model)
stress

plot <- tb2$plot_ordination(plot_color = "solar_season", plot_shape = "deployment", plot_type = c("point", "ellipse"), NMDS_stress_pos = NULL, ellipse_level = 0.95) +
  # geom_arrow_chain(colour = "black") + 
    ggtitle("NMDS Ordination of WM Data", subtitle = paste("Stress = ", round(tb2$res_ordination$model$stress, 2))) + 
    # geom_label_repel(label = mt_16s$sample_table$full_id, color = "black") +
    labs(color = "Season") +
    guides(fill = "none") +
    theme_classic() +
    theme(
      text = element_text(size = 15),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 3)
    )
  plot



load("output/data/mt_16s.RData")

s_mt <- clone(mt_16s)

s_mt$sample_table <- s_mt$sample_table |> 
  filter(lake == "superior") |> 
  drop_na("solar_season")

s_mt$tidy_dataset()

s_mt$cal_betadiv(method = c("bray","aitchison","robust.aitchison", "jaccard"))
s_mt$beta_diversity$bray
# 16S, Robust Aitchison, NMDS and stressplot

tb2 <- trans_beta$new(s_mt, group = "solar_season", measure = "robust.aitchison")
tb2$cal_ordination(method = "NMDS")
tb2$cal_manova(manova_all = TRUE)
print(tb2$res_manova)

tb2 <- trans_beta$new(s_mt, group = "deployment", measure = "bray")
tb2$cal_ordination(method = "NMDS")
tb2$cal_manova(manova_all = TRUE)
print(tb2$res_manova)

stress <- vegan::stressplot(tb2$res_ordination$model)
stress

plot <- tb2$plot_ordination(plot_color = "solar_season", plot_shape = "deployment", plot_type = c("point", "ellipse"), NMDS_stress_pos = NULL, ellipse_level = 0.95) +
  # geom_arrow_chain(colour = "black") + 
    ggtitle("NMDS Ordination of WM Data", subtitle = paste("Stress = ", round(tb2$res_ordination$model$stress, 2))) + 
    # geom_label_repel(label = s_mt$sample_table$full_id, color = "black") +
    labs(color = "Season") +
    guides(fill = "none") +
    theme_classic() +
    theme(
      text = element_text(size = 15),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 3)
    )
  plot



  # Start analysis ------------------------------------------------
nmds <- function(mt_16s, measure = "bray", group = "solar_season",
 title = NULL){ 
  mt <- mt_16s$clone(deep =TRUE)
# First I want to remove the lake erie samples
mt$sample_table <- mt$sample_table |> 
  filter_out(lake == "erie")
mt$tidy_dataset()
# mt$cal_abund()
# Calc Alpha Diversity
# mt$cal_alphadiv()
# Calc Bray Curtis Dissimilarity
mt$cal_betadiv(method = c("bray","aitchison","robust.aitchison", "jaccard"))
# 16S, Robust Aitchison, NMDS and stressplot
tb2 <- trans_beta$new(mt, group = group, measure = measure)
tb2$cal_ordination(method = "NMDS")
# tb2$cal_manova(manova_all = TRUE)
# print(tb2$res_manova)
plot <- tb2$plot_ordination(plot_color = group, plot_shape = "deployment", plot_type = c("point", "ellipse"), NMDS_stress_pos = NULL, ellipse_level = 0.95) +
  # geom_arrow_chain(colour = "black") + 
    ggtitle("NMDS Ordination of WM Data", subtitle = paste("Stress = ", round(tb2$res_ordination$model$stress, 2))) + 
    # geom_label_repel(label = mt$sample_table$full_id, color = "black") +
    labs(color = "Stratification Season") +
    guides(fill = "none") +
    theme_classic() +
    theme(
      text = element_text(size = 15),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 3)
    ) + 
      ggtitle(title)
  return(plot)
  }

  nmds(full_mt_16s, measure = "robust.aitchison", group = "", title = "full")
  nmds(contam_16s,  measure = "robust.aitchison", group = "solar_season", title = "contam")
  nmds(mt_16s,  measure = "bray", group = "solar_season", title = "NMDS of 16S data, split by season (color) and year (shape)")


# Straight abundances-----------------------------------------------------------


load("output/data/mt_16s.RData")

s_mt <- clone(mt_16s)

s_mt$sample_table <- s_mt$sample_table |> 
  filter(lake == "superior") |> 
  mutate(month = month(date))

for(i in 1:12){
  s_mt <- clone(mt_16s)
  s_mt$sample_table <- s_mt$sample_table |> 
  filter(lake == "superior") |> 
  mutate(month = month(date)) |> 
    filter(month == i)
  s_mt$tidy_dataset()

  a <- trans_abund$new(s_mt, taxrank = "c", ntaxa = 20)$plot_bar(facet = "deployment") +
  ggtitle(paste0("Month:", i))
  print(a)
}

  s_mt <- clone(mt_16s)
  s_mt$sample_table <- s_mt$sample_table |> 
  filter(lake == "superior") |> 
  mutate(month = month(date)) 

t <- trans_abund$new(s_mt, taxrank = "c", ntaxa = 10) 
season_df <- 
  t$data_abund |> 
  filter(Taxonomy %in% t$data_taxanames) |> 
  group_by(deployment, strat_season) |> 
  mutate(date_max = max(date), 
          date_min = min(date)) |> 
  ungroup()
p1 <- t$data_abund |> 
  filter(Taxonomy %in% t$data_taxanames) |> 
  ggplot(aes(x = date, y = Taxonomy, fill = Abundance/100)) +
  geom_tile(aes(color = Abundance / 100)) + 
  scale_fill_viridis_c(transform = "log10", option = "magma", name = "Relative Abundance (%)", aesthetics = c("fill", "color")) +
  scale_x_date(
    date_breaks = "1 month",   # Interval spacing (e.g., "2 weeks", "1 year")
    date_labels = "%b %Y", 
    expand = c(0,0)
         # Strftime formatting (e.g., "Jan 2026")
  ) + 
  scale_y_discrete(expand = c(0,0)) +
    theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
axis.title.x = element_blank(),
    plot.margin = margin(t = 0, r = 5, b = 5, l = 5))

p2 <- season_df |> 
  ggplot(aes(date, y = "Season", fill = strat_season)) +
  geom_tile(aes(color = strat_season)) +
  scale_fill_brewer(palette = "Dark2", aesthetics = c("fill", "color"), name = "Season") +
scale_x_date(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_void() +
  theme(
    plot.margin = margin(t = 5, r = 5, b = 0, l = 5)
  )
p2

p2/p1 + patchwork::plot_layout(heights = c(1, 12))

# Above is a nice heat map
  plotly::plotly_build()

a
a <- trans_abund$new(s_mt, taxrank ="s", ntaxa = 100)$plot_bar(facet = "solar_season") |> 
  plotly::plotly_build()

a

a <- trans_abund$new(s_mt, taxrank ="s", ntaxa = 100)$plot_bar(facet = "deployment") |> 
  plotly::plotly_build()

a


s_mt$sample_table |> 
  view()


load("output/data/mt_16s.RData")

mt_16s$sample_table <- mt_16s$sample_table |> 
  filter_out(lake == "erie")

mt_16s$filter_taxa(
  rel_abund = 0.0001
)

mt_16s$tidy_dataset()

mt_16s$rarefy_samples()

mt_16s
v <- trans_venn$new(mt_16s)

v$plot_bar()


