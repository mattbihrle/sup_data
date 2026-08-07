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
b <- betaturn$new(mt_16s, measure = "betaNTI", filter_thres = 0.001)
save(b, file = "output/data/betaturn.RData")

load("output/data/betaturn.RData")

# okay now try to calculate the distance based on what I learned last week

b$cal_group_distance(group = "date", within_group = FALSE)
b$res_group_distance |> head()
b$cal_group_distance_diff(method = "KW_dunn")
b$res_group_distance_diff
b$plot_group_distance() +
  geom_jitter()
b$dataset$beta_diversity$betaNRI |> 
  view()

test <- clone(b)
test$res_group_distance <- test$res_group_distance |> 
  filter(deployment == "s_2122_WM")
test$plot_group_distance()
b$res_ses_betamntd


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
  mutate(days_between = as.numeric(difftime(end, start, units = "days"))) |>
  mutate(start = as_date(start),
         end = as_date(end)) |> 
  # mutate(dep_start = case_when(start %in% filter(mt_16s$sample_table, deployment == "s_2122_WM")$date ~ "s_2122_WM",
  #                              start %in% filter(mt_16s$sample_table, deployment == "s_2223_WM")$date ~ "s_2223_WM",
  #                              start %in% filter(mt_16s$sample_table, deployment == "s_2324_WM")$date ~ "s_2324_WM",
  #                              start %in% filter(mt_16s$sample_table, deployment == "s_2425_WM")$date ~ "s_2425_WM")) |> 
  # mutate(dep_end = case_when(end %in% filter(mt_16s$sample_table, deployment == "s_2122_WM")$date ~ "s_2122_WM",
  #                            end %in% filter(mt_16s$sample_table, deployment == "s_2223_WM")$date ~ "s_2223_WM",
  #                            end %in% filter(mt_16s$sample_table, deployment == "s_2324_WM")$date ~ "s_2324_WM",
  #                            end %in% filter(mt_16s$sample_table, deployment == "s_2425_WM")$date ~ "s_2425_WM")) |> 
  # mutate(strat_start = case_when(start %in% filter(mt_16s$sample_table, strat_season == "Summer")$date ~ "Summer",
  #                                start %in% filter(mt_16s$sample_table, strat_season == "Fall")$date ~ "Fall",
  #                                start %in% filter(mt_16s$sample_table, strat_season == "Winter")$date ~ "Winter",
  #                                start %in% filter(mt_16s$sample_table, strat_season == "Spring")$date ~ "Spring")) |> 
  # mutate(strat_end = case_when(end %in% filter(mt_16s$sample_table, strat_season == "Summer")$date ~ "Summer",
  #                              end %in% filter(mt_16s$sample_table, strat_season == "Fall")$date ~ "Fall",
  #                              end %in% filter(mt_16s$sample_table, strat_season == "Winter")$date ~ "Winter",
  #                              end %in% filter(mt_16s$sample_table, strat_season == "Spring")$date ~ "Spring")) |> 
  left_join(mt_16s$sample_table, by = join_by(start == date), suffix =  c("_start", "")) |> 
  left_join(mt_16s$sample_table, by = join_by(end == date), suffix =  c("", "_end")) |> 
  mutate(date = start) |> 
  left_join(summer_start, by = join_by(date, closest(start >= first_strat)))


a <- plot_df |> 
  mutate(temp_diff = abs(temp_c_24h_med_pps - temp_c_24h_med_pps_end)) |> 
  mutate(days_from_strat =)
  # filter(abs(Value) >= 2) |> 
  # filter(deployment == "s_2122_WM" & deployment_end == "s_2122_WM") |>
  ggplot(aes(days_between, Value, z = start, w = end, color = mixing, fill = mixing_end)) +
  geom_point(shape = 21, size = 4, stroke = 2) +
  scale_color_discrete(palette = "Set2") +
  scale_fill_discrete(palette = "Set2") +
  theme_classic() +
  NULL
a
ap <- plotly::plotly_build(a)
ap

b$res_group_distance
