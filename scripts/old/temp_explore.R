# Try to plot the temperature at the depth we need
# Load packages
library(tidyverse)
library(plotly)
# Run 01_intro_script.R to get sw_meta object

# Load Temp data ---------------------------------------------------------------------------
temp_df <- mat_df$T |> 
  as.data.frame()
# Set colnames and rownames
colnames(temp_df) <- mat_df$dep  
rownames(temp_df) <- mat_df$t 
# move dttm column out, convert to datetime and pivot long
# temp_df <- temp_df |> 
#   rownames_to_column(var = "dttm") |>
#   mutate(dttm = as_datetime((as.numeric(dttm) - 719529) * 86400)) |> 
#   pivot_longer(cols = !dttm, names_to = "depth", values_to = "temp")

# add 24 hr temp to metadata



## plot temp data ------------------------------------------------------------------------------
temp_plot <- temp_df_long |> 
  filter(depth < 50) |> 
  ggplot(aes(dttm, temp)) +
  geom_vline(xintercept = sw_meta$date) +
    geom_point(size = 0.5, color = "grey1")
    

temp_plot
plotly_build(temp_plot)


# Bill Plot ------------------------------------------------------------------
# temp plots
temp_df_plot <- temp_df_long %>% 
  filter(depth < 50) %>% 
  mutate(date = as_date(date)) %>% 
  select(date, depth, temp) %>% 
  mutate(week_date = floor_date(date, unit = "week"))%>% 
  # Group by that new weekly date column
  group_by(week_date, depth) %>% 
  summarize(mean_temp = mean(temp, na.rm = TRUE), .groups = "drop")


temp_plot <- temp_df_plot |>
  ggplot(aes(week_date, mean_temp, color = as.character(depth))) +
    geom_point(size = 0.5) +
  geom_line(linewidth = 0.5)+
    geom_vline(xintercept = sw_meta$date)
temp_plot
