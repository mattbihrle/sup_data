# just a little one to source
library(mecoturn)
library(tidyverse)
library(microeco)

load("output/data/mt_16s.RData")

mt_16s$phylo_tree <- ape::read.tree("data/16s/lotus3_out/OTUphylo.nwk")
mt_16s$sample_table <- mt_16s$sample_table |> 
  filter(lake != "erie") |> 
    # filter(deployment == "s_2122_WM") |> 
  mutate(month = month(date, label = TRUE))

  mt_16s
  mt_16s$tidy_dataset()
  # In this module, I want to try and caculate the days since stratification
  dates <- mt_16s$sample_table |> 
    filter(strat_season == "Summer") |> 
    mutate(year = year(date)) |> 
    group_by(year) |> 
    mutate(fall_date = as_date(max(date) + weeks(2))) |> 
    ungroup() |> 
    select(fall_date) |> 
    unique() 
  
  mt_16s$sample_table <- mt_16s$sample_table |> 
    rownames_to_column("rows") |> 
    left_join(dates, by = join_by(closest(date >= fall_date))) |> 
    mutate(days_since_fall = as.character(as.numeric(as_date(date) - as_date(fall_date)))) |> 
    drop_na(days_since_fall) |> 
    mutate(days_since_fall = str_pad(days_since_fall, side = "left", width = 3, pad = "0")) |> 
      column_to_rownames("rows")

    # split temperature into bins 
  mt_16s$sample_table <- mt_16s$sample_table |>
    mutate(temp_bin = cut(mt_16s$sample_table$temp_c_24h_med_pps, 
                          breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8)))
b <- betaturn$new(mt_16s, measure = "betaNTI", filter_thres = 0.001)
save(b, file = "output/data/betaturn.RData")

load("output/data/betaturn.RData")


# okay now try to calculate the distance based on what I learned last week

b$cal_group_distance(group = "date", within_group = FALSE)

b$res_group_distance |> head(15)
b$cal_group_distance_diff(method = "KW_dunn")
b$res_group_distance_diff
a <- b$plot_group_distance() +
  geom_jitter()
a
b$dataset$beta_diversity$betaNRI |> 
  view()




summer_start <- mt_16s$sample_table |> 
  mutate(year = year(date)) |> 
  group_by(year) |> 
  # filter(strat_season == "Summer") |> 
  mutate(first_strat = min(date[which(mt_16s$sample_table$strat_season == "Summer")])) |> 
  ungroup() |> 
  select(date, first_strat) 

# Make a plot of days between samples
plot_df <- b$res_group_distance |> 
  separate_wider_delim(cols = "date", delim = " vs ", names = c("start", "end")) |>
  mutate(across(start:end, ~as_datetime(.x))) |>
  mutate(start = as_date(start),
         end = as_date(end)) |> 
    left_join(mt_16s$sample_table, by = join_by(start == date), suffix =  c("_start", "")) |> 
    left_join(mt_16s$sample_table, by = join_by(end == date), suffix =  c("", "_end")) |> 
    mutate(date = start) |> 
    mutate(days_between = as.numeric(difftime(end, start, units = "days"))) |> 
  # Make a column that compares how different the samples are based on where they fall in the year
    mutate(fall_days_between = abs(as.numeric(days_since_fall) - as.numeric(days_since_fall_end)))


a <- plot_df |> 
  # mutate(temp_diff = abs(temp_c_24h_med_pps - temp_c_24h_med_pps_end)) |> 
  # mutate(days_from_strat =)
  # filter(abs(Value) >= 2) |> 
  filter(deployment == "s_2122_WM") |>
  filter(Value >= 2) |> 
  # mutate(is_sample = ifelse(full_id %in% c("s_2122_WM17", "s_2122_WM15", "s_2122_WM07"), TRUE, FALSE)) |> 
  ggplot(aes(x = fall_days_between, y = Value, color = strat_season, Value, z = start, w = end, p = strat_season, q = full_id, fill = full_id_end)) +
  geom_point(shape = 21, size = 4, stroke = 2) +
  # scale_color_discrete(palette = "Set2") +
  # scale_fill_discrete(palette = "Set2") +
  # scale_fill_viridis_d(option = "turbo") +
  # scale_color_viridis_d() +
  theme_classic() +
  NULL
a
ap <- plotly::plotly_build(a)
ap

b$res_group_distance

t <- clone(mt_16s)


t$sample_table <- t$sample_table |> 
  filter(deployment == "s_2122_WM") 
t$tidy_dataset()

trans_abund$new(t, ntaxa = 50, taxrank = "o")$plot_bar(facet = "strat_season")
# Make a plot of days between "days since fall" groups
trans_beta$new(t, group = "strat_season", measure = "robust.aitchison")$
  cal_ordination(method = "NMDS")$
  plot_ordination(add_sample_label = "full_id")

plot_df <- b$res_group_distance |> 
  separate_wider_delim(cols = "days_since_fall", delim = " vs ", names = c("start", "end")) |>
  mutate(start = as.numeric(start), 
        end = as.numeric(end)) |> 
  mutate(days_between = abs(end - start))
  # left_join(mt_16s$sample_table, by = join_by(start == date), suffix =  c("_start", "")) |> 
  # left_join(mt_16s$sample_table, by = join_by(end == date), suffix =  c("", "_end")) |> 
  # mutate(date = start)


a <- plot_df |> 
  # mutate(temp_diff = abs(temp_c_24h_med_pps - temp_c_24h_med_pps_end)) |> 
  # mutate(days_from_strat =)
  # filter(abs(Value) >= 2) |> 
  # filter(deployment == "s_2122_WM" & deployment_end == "s_2122_WM") |>
  ggplot(aes(days_between, Value, color = start, fill = end)) +
  geom_point(shape = 21, size = 4, stroke = 2) +
  # scale_fill_di(palette = "Set2") +
  # scale_color_viridis_c(option = "turbo") +
  theme_classic() +
  NULL
a
ap <- plotly::plotly_build(a)
ap
