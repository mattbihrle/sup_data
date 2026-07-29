# 01.1 Import PPS Data


# Load Libraries ----------------------------------------------
# List packages to load
packages <- c("tidyverse", "googlesheets4", "patchwork")
# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
# Packages loading
invisible(lapply(packages, library, character.only = TRUE))
# Remove the list of packages and installed packages from the environment
rm(installed_packages, packages)

# Import meta data from 2024-2025 year -----------------------------------------
sw_meta <- read_csv("output/data/metadata_supwinter.csv") |> 
  dplyr::select(!matches(c("12h_med", "round_date")))|> 
  distinct() |> 
  mutate(full_id = paste0("s_2425_", sample))

# Import master spreadsheet from google
# gs4_auth()
ext_sheet <- read_sheet("https://docs.google.com/spreadsheets/d/1Iq5AIxiH-EBa2lT83cKCUsLwYCvvCl0iHlCbkJw9rQo/edit?gid=1015683691#gid=1015683691") |> 
  # make sure all columns are the correct type 
  mutate(across(contains("dttm"), ~parse_date_time(.x, orders = c("dmyHMS", "mdyHMS"))),
         across(contains("date"), ~mdy(.x))
          )
# Create a sample table for import into the microtable
sample_table <- sw_meta |> 
  full_join(ext_sheet, by = "full_id") |> 
  mutate(rows = full_id) |> 
    column_to_rownames("rows") |> 
    as.data.frame() |> 
  mutate(date = as_date(dttm_sampled))

sample_table <- sample_table |> 
  mutate(strat_season = fct_reorder(strat_season, date, .na_rm = F), 
         strat_season_2 = fct_reorder(strat_season_2, date, .na_rm = F),
         solar_season = fct_reorder(solar_season, date, .na_rm = F), 
         mixing = fct_relevel(mixing, "stratified_std", "mixed", "stratified_inverse"))

# Now add additional data ------------------------------------------------

files <- list.files("data/metadata/pps/thermistor_data", pattern = ".asc", recursive = T, full.names = T)

full_t_df <- map(files, ~
read_lines(.x) |> 
  str_squish() |>
  str_subset("^[0-9]") |> 
  as.data.frame() |> 
  separate_wider_delim(cols = everything(), delim = ",", 
                       names = c("temp_c", "pres_db", "date_utc", 
                                  "time_utc"), 
                        too_many = "error", too_few = "error") |> 
  mutate(across(everything(), ~str_squish(.x)), 
        xml_file = str_replace(.x, ".asc", ".XML"),
        sn_inst = na.omit(str_extract(read_lines(.x), "(?<=SerialNumber: )[0-9]*"))
  )
) |> 
  list_rbind()

full_t_df <- full_t_df |> 
  mutate(across(c(temp_c, pres_db), ~as.numeric(.x)), 
          dttm_utc = round_date(dmy_hms(paste(date_utc, time_utc)), unit = "hour")) |> 
  select(-date_utc, -time_utc)

full_t_df |> 
  ggplot(aes(dttm_utc, temp_c)) +
  geom_point() +
  facet_wrap(~xml_file, scales = "free")
# Now that I have a master table of temperatures, I want to pull out the 24 hr median temp and pressure for each sample

t_df <- full_t_df |> 
  mutate(date = as_date(ceiling_date(dttm_utc, unit = "day"))) |> 
  group_by(xml_file, date) |> 
  mutate(temp_c_24h_med_pps = median(temp_c), 
          pres_db_24h_avg_pps = mean(pres_db)) |> 
  ungroup() |> 
  # select only the dates that match a sampling date
  filter(date %in% sample_table$date) |> 
  select(date, temp_c_24h_med_pps, pres_db_24h_avg_pps, sn_inst, xml_file) |> 
  distinct() |> 
  #add a column that pulls out the deployment code
  mutate(lake = str_to_lower(str_extract(xml_file, "(?<=Lake)(Superior|Erie)")), 
          loc = str_extract(xml_file, "(LOR|WM|ACB|Lor)")) |> 
  group_by(xml_file) |> 
  mutate(yr_start = min(year(date)),
          yr_end = yr_start + 1) |> 
  mutate(deployment = paste0(str_extract(lake, "^[a-z]{1}"), 
                              "_", 
                              str_extract(yr_start, "[0-9]{2}$"),
                              str_extract(yr_end, "[0-9]{2}$"), 
                              "_", 
                              str_to_upper(loc))) |> 
  ungroup() |> 
  mutate(deployment = str_replace(deployment, "NA", "e")) |> 
  select(-c(lake:yr_end))
stop()
# Now join to the other sheet

sample_table <- sample_table |> 
  mutate(deployment = str_extract(full_id, "^(e|s)_[0-9]{4}_[A-Z]*")) |> 
    left_join(t_df, by = c("date", "deployment")) |> 
  group_by(deployment) |> 
    mutate(mooring_depth = mean(pres_db_24h_avg_pps)) |> 
  ungroup()




# Quick graph to see how the years change 

p1 <- sample_table |> 
  drop_na(date) |> 
  group_by(deployment) |> 
  mutate(med_t = median(temp_c_24h_med_pps)) |> 
  ungroup() |> 
ggplot(aes(date, temp_c_24h_med_pps, color = deployment)) +
  geom_point() +
    scale_fill_viridis_d() +
    scale_x_date(date_breaks = "1 month") +
      stat_smooth(formula = y~1, aes(group = 1), method = "lm", se = FALSE) +
      # stat_summary(fun = median, geom = "hline", aes(yintercept = after_stat(y))) +
      facet_wrap(~deployment, scales = "free") +
geom_hline(aes(yintercept = med_t, group = deployment)) + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1))
# load("output/data/temp_df_long_full.RData")

# Add temp labels
temp_labels_1 <- seq(0, 17, by = 1)
temp_labels_2 <- seq(1, 18, by = 1)

temp_labels <- paste0(temp_labels_1, "-", temp_labels_2)

p2 <- temp_df_long_full |> 
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
  geom_point(data = filter(sample_table, deployment == "s_2425_WM"), aes(x = date, y = 38, z = NULL), color = 'white') +
    theme_classic() +
    theme(
      text = element_text(size = 15),
      panel.border = element_rect(colour = "black", fill = "red", linewidth = 3)
    ) +
  scale_y_reverse(expand = c(0,0), n.breaks = 10) +
  scale_x_date(expand = c(0,0), date_breaks = "1 month", date_label = "%b", limits = c(as_date("2024-07-20"), as_date("2025-08-15"))) +
    guides(fill = guide_legend(ncol = 1, reverse = T)) 

p1/p2

# Next up I could write some code for the superior samples to determine summer vs mixed vs winter. 

