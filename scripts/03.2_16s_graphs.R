# 16S Abundances

# Load data
load("output/data/mt_16s.RData")


# Load Libraries 
packages <- c("plotly", "ggalluvial", "ggarrow", "tidyverse", "microeco", "ggrepel")

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
mt_16s

# First I want to remove the lake erie samples
mt_16s$sample_table <- mt_16s$sample_table |> 
  filter_out(lake == "erie")
mt_16s$tidy_dataset()
mt_16s

mt_16s$cal_abund()
mt_16s$taxa_abund$p[1:5,1:5]

# Calc Alpha Diversity
mt_16s$cal_alphadiv()

# Calc Bray Curtis Dissimilarity
mt_16s$cal_betadiv(method = c("bray","aitchison","robust.aitchison", "jaccard"))
mt_16s$beta_diversity$bray
# 16S, Robust Aitchison, NMDS and stressplot

tb2 <- trans_beta$new(mt_16s, group = "strat_season", measure = "bray")
tb2$cal_ordination(method = "NMDS")
tb2$cal_manova(manova_all = TRUE)
print(tb2$res_manova)

stress <- vegan::stressplot(tb2$res_ordination$model)
stress

plot <- tb2$plot_ordination(plot_color = "strat_season", plot_shape = "deployment", plot_type = c("point", "ellipse"), NMDS_stress_pos = NULL, ellipse_level = 0.95) +
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
  drop_na(strat_season)

s_mt$tidy_dataset()

s_mt$cal_betadiv(method = c("bray","aitchison","robust.aitchison", "jaccard"))
s_mt$beta_diversity$bray
# 16S, Robust Aitchison, NMDS and stressplot

tb2 <- trans_beta$new(s_mt, group = "strat_season", measure = "robust.aitchison")
tb2$cal_ordination(method = "NMDS")
tb2$cal_manova(manova_all = TRUE)
print(tb2$res_manova)

tb2 <- trans_beta$new(s_mt, group = "deployment", measure = "bray")
tb2$cal_ordination(method = "NMDS")
tb2$cal_manova(manova_all = TRUE)
print(tb2$res_manova)

stress <- vegan::stressplot(tb2$res_ordination$model)
stress

plot <- tb2$plot_ordination(plot_color = "strat_season", plot_shape = "deployment", plot_type = c("point", "ellipse"), NMDS_stress_pos = NULL, ellipse_level = 0.95) +
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
