# Need temp labels and samp_dates from a previous script
# Also need to import the network data

library(tidyverse)
library(microeco)
library(patchwork)
load("output/data/temp_df_long_full.RData")
load("output/data/module_long_16s.RData")
load("output/data/maestro_df_clean.RData")


# Add temp labels
temp_labels_1 <- seq(0, 17, by = 1)
temp_labels_2 <- seq(1, 18, by = 1)

temp_labels <- paste0(temp_labels_1, "-", temp_labels_2)

p1 <- temp_df_long_full |> 
    ungroup() |> 
    # select(dttm, depth, temp) |> 
    drop_na(temp) |> 
    mutate(depth = as.numeric(depth)) |> 
    filter(depth < 50) |> 
        dplyr::select(date, depth, temp) |> 
        # drop_na(temp) |> 
        distinct() |> 
    ggplot(aes(x = date, y = depth, z = temp)) +
      # geom_contour_filled(breaks = breaks) +
      geom_contour_filled(binwidth = 1) +
        # scale_fill_brewer(palette = "RdBu", direction = -1) +
        scale_fill_viridis_d(option = "turbo", labels = temp_labels) +
  labs(x = "Date", y = "Depth (m)", fill = "Temperature (°C)") +
  geom_point(data = mt_mag$sample_table, aes(x = date, y = 38, z = NULL), color = 'white') +
    theme_classic() +
    theme(
      text = element_text(size = 15),
      panel.border = element_rect(colour = "black", fill = "red", linewidth = 3)
    ) +
  scale_y_reverse(expand = c(0,0), n.breaks = 10) +
  scale_x_date(expand = c(0,0), date_breaks = "1 month", date_label = "%b", limits = c(as_date("2024-07-20"), as_date("2025-08-05"))) +
    guides(fill = guide_legend(ncol = 1, reverse = T)) 
p1
p2 <- ggplot(module_long, aes(x = date, y = relabund, color = Module, group = Module, label = rownames(module_long))) +
  # geom_rect(aes(xmin = sum_start, xmax = fall_start, ymax = max(relabund) +4 , ymin = max(relabund)+1),  fill = "grey", color = NULL) +
  # geom_rect(aes(xmin = fall_start, xmax = winter_start, ymax = max(relabund) +4 , ymin = max(relabund)+1),  fill = "green4",  color = NULL) +
  # geom_rect(aes(xmin = winter_start, xmax = spring_start, ymax = max(relabund) +4 , ymin = max(relabund)+1),  fill = "grey",  color = NULL) +
  # geom_rect(aes(xmin = spring_start, xmax = max(date), ymax = max(relabund) +4 , ymin = max(relabund)+1),  fill = "green4",  color = NULL) +
  geom_line(size = 1) +
  geom_text() +
  geom_point(size = 2) +
  labs(title = "Module Abundance Over Time",
       y = "Relative Abundance") +
  theme_bw() +
  scale_x_date(expand = c(0,0), date_breaks = "1 month", date_label = "%b", limits = c(as_date("2024-07-20"), as_date("2025-08-05")))

p2

p3 <- maestro_df_clean |> 
  select(dttm_local, matches(c("cdom_12"))) |> 
  pivot_longer(cols = -dttm_local, names_to = "var", values_to = "value") |> 
    ggplot(aes(dttm_local, value)) +
  geom_point() + 
  facet_grid(var ~ ., scales = "free") +
  scale_x_date(expand = c(0,0), date_breaks = "1 month", date_label = "%b", 
               limits = c(as_date("2024-07-10"), as_date("2025-08-05"))) +
  theme_classic()
p1 / p2 / p3 + plot_layout(heights = c(1, 2, 0.5))

