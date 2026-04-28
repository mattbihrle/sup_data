# For loading .mat file
library(R.matlab)
library(tidyverse)
library(performance) # for outlier checks
library(skimr) # For summary stats
library(plotly) # For interactive plots
library(LakeMetabolizer) # For oxygen data
library(respR) # For working with oxygen data
 
# Load metadata-------------------------------------------------------------------------
sw_meta <- read_csv("data/sw_metadata.csv") |>
  mutate(date = mdy(date)) |> 
  select(!day:year)

# Create secondary stratification season based on combining fall and winter toegether
sw_meta <- sw_meta |> 
  mutate(strat_season_2 = ifelse(strat_season == "fall", "winter", strat_season))
sw_meta$strat_season <- factor(sw_meta$strat_season, levels = c("summer", "fall", "winter", "spring"), ordered = T)

sw_meta$strat_season_2 <- factor(sw_meta$strat_season_2, levels = c("summer", "winter", "spring"), ordered = T)
# Create another stratification based on the equinox

sw_meta <- sw_meta |> 
  mutate(solar_season = ifelse(date < "2024-12-21", "summer", "winter")) |> 
  mutate(solar_season = ifelse(solar_season == "summer" & date >= "2024-09-22", "fall", solar_season)) |> 
  mutate(solar_season = ifelse(solar_season == "winter" & date >= "2025-03-20", "spring", solar_season)) |> 
  mutate(solar_season = ifelse(date > "2025-06-20", "summer", solar_season))
#import temp file ----------------------------------------------------------------------------
mat_df <- readMat.default(con = "data/metadata/LSEO24h_temp.mat")

# Extract dataframe we want from 'T' for temp and convert matlab time to poxit time
temp_df <- mat_df$T |> 
  as.data.frame()
# Set colnames and rownames
colnames(temp_df) <- mat_df$dep  
rownames(temp_df) <- mat_df$t 

# Create date column 
temp_df <- temp_df |> 
    rownames_to_column(var = "dttm") |>
  mutate(dttm = as_datetime((as.numeric(dttm) - 719529) * 86400)) |> 
    # round to nearest hour to make things easier
  mutate(dttm = round_date(dttm, unit = "hour"))
  #   # Select only columns less than 50m
  # select(which(as.numeric(colnames(temp_df)) <= 50))
# move dttm column out, convert to datetime and pivot long
temp_df_long <- temp_df |> 
  # rownames_to_column(var = "dttm") |>
  # mutate(dttm = as_datetime((as.numeric(dttm) - 719529) * 86400)) |> 
  # # round to nearest hour to make things easier
  # mutate(dttm = round_date(dttm, unit = "hour")) |> 
  pivot_longer(cols = !dttm, names_to = "depth", values_to = "temp") |>
  # Create variable of local time
  mutate(dttm_local = with_tz(dttm, tzone = "US/Central")) |> 
  # Add averaged temperature column
  mutate(date = floor_date(dttm_local, unit = "day")) |> 
  group_by(date, depth) |> 
  mutate(temp_c_24h_med = median(temp)) |> 
  ungroup()

# First create temporary df of sample dates
  samp_dates <- sw_meta |> 
    select(sample) |> 
    mutate(sample_date = as_date(sw_meta$date))

## Create a 14 day averaged temperature---------------------------------
  
  # Then create 14 day average
temp_df_long <- temp_df_long |>
  # Force the date column to be a simple Date, dropping any hours/minutes/seconds
  mutate(date = as_date(date)) |>
  # Match every temp date to the upcoming sample date
  left_join(samp_dates, join_by(closest(date <= sample_date))) |>
  # Test to see how far apart sampling dates are
  mutate(test = as.numeric(sample_date - date)) |> 
  # NEW STEP: Force a strict 14-day window!
  # Keep only the rows where the temp date is within 14 days of the sample date
  filter(test <= 14) |>
  # Calculate the bin average
  group_by(sample_date, depth) |> 
  mutate(temp_c_14d_avg = mean(temp, na.rm = TRUE)) |> 
  ungroup()

# OLD VERSION
# locs <- seaprocess::find_near(vec = as_date(temp_df_long$date), vals = as_date(sw_meta$date))
# dates <- as_date(rep(NA, nrow(temp_df_long)))
# for(i in length(locs):1) {
# dates[1:locs[i]] <- as_date(sw_meta$date[i])
# }
# rm(i)
# # Add the 14ish day average
# temp_df_long<-temp_df_long |>
#   mutate(sample_date = dates) |> 
#   group_by(sample_date, depth) |> 
#   mutate(temp_c_14d_avg = mean(temp)) |> 
#   ungroup()
# Create a coefficient of variation to see if we should average or not

temp_df_long <- temp_df_long |> 
  group_by(date) |> 
  mutate(coeff_var = (sd(temp, na.rm = T)/mean(temp, na.rm = T)) * 100)

  # Plot coeff_var only on sample dates
temp_df_long |> 
  filter(sample_date == date, depth == 38) |> 
  mutate(coeff_test = ifelse(coeff_var > 10, "FALSE", "TRUE")) |> 
  ggplot(aes(dttm_local, temp, color = coeff_test)) +
  geom_point() +
  facet_wrap(~sample_date, scales = "free") +
  scale_x_datetime(date_labels = "%H")



# Plot temperature profiles on each sampling date
  # First only keep synoptic hours
hours_keep = c(00, 06, 12, 18)
# The create the plot
temp_df_long |> 
  filter(sample_date == date, minute(dttm_local) == 0, hour(dttm_local) %in% hours_keep) |> 
  ggplot(aes(x = temp, y = as.numeric(depth), color = as.character(hour(dttm_local)))) + 
  geom_point() + 
  scale_color_discrete() +
  scale_y_reverse() +
  facet_wrap(~ sample_date, scales = "free") +
  geom_hline(yintercept = 38) +
  NULL
  

### Add averaged temp column to metadata ---------------------
sw_meta_final <- temp_df_long |> 
  filter(depth == 38) |> 
  select(date, depth, temp_c_24h_med) |> 
  distinct(date, .keep_all = T) |> 
  right_join(y = sw_meta, by = "date") |> 
  relocate(sample)

# Import maestro data --------------------------------------------------------------------
# Recreate a dttm column
maestro_df <- readMat.default(con = "data/metadata/lseo_maestro_time.mat") |> 
  data.frame() |> 
  mutate(dttm = as_datetime((as.numeric(t) - 719529) * 86400), .before = CDOM) |>
  # round dttm to nearest minute
  mutate(dttm = round_date(dttm, unit = "minute")) |> 
  # Create a local dttm
  mutate(dttm_local = with_tz(dttm, tzone = "US/Central")) |> 
  # Create a date column rounded to nearest whole day for averaging later
  mutate(date = as_date(round_date(dttm_local, unit = "day"))) |> 
  # filter out the beginning of the deployment where pressure is less than 20m 
  filter(P > 20) |> 
  # add a sample date column for future averaging
  left_join(samp_dates, by = join_by(closest(date <= sample_date)))

  # Create summary stats to check for nas
sum_stats <- skim(maestro_df)
sum_stats

  # Create long df for plotting
long_maestro_df <- maestro_df |> 
  pivot_longer(cols = -c(dttm, date, sample, sample_date, dttm_local),
  names_to = "variable", values_to = "vals", values_drop_na = T)

  # Quick plot
long_maestro_df |> 
  ggplot(aes(dttm_local, vals, color = variable)) +
  geom_point()+
  facet_wrap(~variable, scales = "free")


# Create a vector of data variables to use later
data_vars <- maestro_df |> 
  select(-c(dttm, date, sample_date, sample, t, dttm_local)) |> 
  colnames() |> 
  as.character()
data_vars

## Start looking for outliers ---------------------------------------------------------
## Plot of outliers - multivariate ----------------------------------------------
# outliers <- maestro_df |> 
#   select(data_vars) |>  
#   drop_na() |> 
#   check_outliers(method = "mcd", verbose = T)
# outliers

# # Make a lil plot
# maestro_df |> 
#   select(c("dttm_local", data_vars)) |> 
#   drop_na() |> 
#   mutate(outlier = outliers) |> 
#   pivot_longer(-c(outlier, dttm_local), names_to = "variable", values_to = "vals") |> 
#   ggplot(aes(dttm_local, vals, color = outlier)) +
#   geom_point() +
#   facet_wrap(~variable, scales = "free")

# Based on looking at the individual variables (below) and conversations with Jay, I am going
# to remove all data after 2025-06-21 00:30:00

maestro_df <- maestro_df |> 
  filter(dttm_local < "2025-06-21 00:30:00")

##  Create otlier function ------------------------------------------------------------------

plot_outliers <- function(maestro_df, var, plotly = F) {
  outlier_df <- maestro_df |> 
    select(dttm_local,  var) |> 
    rename(data = var) |> 
    drop_na()
  # Create a new column
  outlier_df <- outlier_df |> 
    mutate(outlier = check_outliers(outlier_df, method = c("zscore", "iqr")))

  # Make a plot of outliers
  outlier_plot <- 
  outlier_df |> 
  ggplot(aes(dttm_local, data, color = outlier)) + 
  geom_point() +
  ggtitle(paste( var, "outliers"))
  
  outlier_plotly <- plotly_build(outlier_plot)
    if(plotly) {
    return(outlier_plotly)
    } else {
    return(outlier_plot)
    }
  }
### Clean CDOM ------------------------------------------------------------------

plot <- plot_outliers(maestro_df, "CDOM")
#plot

maestro_df_clean <- maestro_df
#Try outliers plot again
plot <- plot_outliers(maestro_df_clean, "CDOM")
#plot

# Looks good moving onto the next

### Clean CHL ------------------------------------------------------------------

plot <- plot_outliers(maestro_df, "CHL")
#plot

# Based on the plot there are some obvious seeming outliers at the bottom of the data
# Removing those

maestro_df_clean <- maestro_df_clean |> 
  mutate(CHL = ifelse(CHL > 0.5, CHL, as.numeric("NA")))

#Try outliers plot again
plot <- plot_outliers(maestro_df_clean, "CHL")
#plot

# Looks better moving on to the next

### Clean CND ------------------------------------------------------------------

plot <- plot_outliers(maestro_df, "CND")
#plot
# Conductivity looks good

### Clean DO ------------------------------------------------------------------

plot <- plot_outliers(maestro_df, "DO")
#plot

# DO looks like there is some craziness at the beginning but might just be diurnal cycles
# While the water column is stratified

# Plot just the beginning

maestro_df |> 
  filter(dttm_local < "2024-09-15") |> 
  mutate(day = as.logical(ifelse(hour(dttm_local) >= 7 & hour(dttm_local) <= 19, "TRUE", "FALSE"))) |> 
  ggplot(aes(dttm_local, DO, color = day)) +
  geom_point()

# Plot 3 days
maestro_df |> 
  filter(date(dttm_local) >= "2024-09-16" & date(dttm_local) <= "2024-09-20") |> 
  mutate(morning = as.logical(ifelse(hour(dttm) >= 12 & hour(dttm) <= 23, "TRUE", "FALSE"))) |> 
  ggplot(aes(dttm, DO)) +
  geom_line()

# No obvious pattern, going to calculate % saturation and compare with temp
maestro_df <- maestro_df |> 
  mutate(do_mgl = convert_DO(maestro_df$DO, from = "umol/L", to = "mg/L"))

  # Pull out a small df to use with lakemetabolizer do.sat function
do_df <- maestro_df |> 
  select(dttm_local, T) |> 
  rename("datetime" = dttm_local, "wtr" = T)

# Create another intermediate dataframe
do_sat <- o2.at.sat(ts.data = do_df, altitude = 183, salinity = 0) |> 
  mutate(o2_sat_percent = maestro_df$do_mgl/do.sat*100)

# Add the original column back into the maestro_df
maestro_df_clean <- maestro_df_clean |> 
  mutate(do_sat_percent = do_sat$o2_sat_percent)

# Add do_sat_percent to our 'data_vars list

data_vars <- c(data_vars, "do_sat_percent")

plot <- plot_outliers(maestro_df_clean, "do_sat_percent")
#plot

# Data look pretty good at this depth? Maybe later can get a better sense of trends


### Clean P ------------------------------------------------------------------

plot <- plot_outliers(maestro_df, "P")
#plot

# Pressure doesn't seem to have any crazy differences. I'm going to leave it as is

### Clean PCO2 ------------------------------------------------------------------

plot <- plot_outliers(maestro_df, "PCO2")
#plot

# pco2 goes fricking crazzzzy at the end of the deployment, plot a single day to see wahat is up
maestro_df |> 
  filter(date(dttm) > "2025-05-17" & date(dttm) < "2025-06-17") |> 
  mutate(morning = as.logical(ifelse(hour(dttm) >= 12 & hour(dttm) <= 23, "TRUE", "FALSE"))) |> 
  ggplot(aes(dttm, PCO2, color = morning)) +
  geom_point()

# Look at pco2 past april 
maestro_df |> 
  filter(date(dttm_local) > "2025-04-01")|> 
  mutate(morning = as.logical(ifelse(hour(dttm_local) >= 12 & hour(dttm_local) <= 23, "TRUE", "FALSE"))) |> 
  ggplot(aes(dttm_local, PCO2, color = morning)) +
  geom_point()

# Looks like there's some big shift around May. I may just remove everything since
maestro_df_clean <- maestro_df_clean |> 
  mutate(PCO2 = ifelse(date(dttm_local) <= "2025-05-01" , PCO2, as.numeric("NA")))

plot <- plot_outliers(maestro_df_clean, "PCO2")
#plot

# Looks okay but still there's a huge drop in late april
plot <- maestro_df |> 
  filter(date(dttm_local) > "2025-04-01")|> 
  mutate(morning = as.logical(ifelse(hour(dttm_local) >= 12 & hour(dttm_local) <= 23, "TRUE", "FALSE"))) |> 
  ggplot(aes(dttm_local, PCO2)) +
  geom_line()
#plotly_build(plot)

# It seems like there's not a clear cut-off. I think I'll keep it as is

### Clean PHYC ------------------------------------------------------------------

plot <- plot_outliers(maestro_df, "PHYC")
#plot

# Phycocyanin stays relatively normal and dull. I'm going to leave it as is. 



#Try outliers plot again
plot <- plot_outliers(maestro_df_clean, "PHYC")
#plot

# Looks better moving on to the next

### Clean TURB ------------------------------------------------------------------

plot <- plot_outliers(maestro_df, "TURB")
#plot

# Some pretty high turbidity levels in the middle but should keep those for now. Plus Turbidity
# isn't very interesting

#Try outliers plot again
plot <- plot_outliers(maestro_df_clean, "TURB")
#plot

# Looks better moving on to the next

### Clean pH ------------------------------------------------------------------

plot <- plot_outliers(maestro_df, "pH")
#plot

test <- maestro_df_clean |> 
  select(dttm_local, pH) |> 
  filter(dttm_local > "2024-12-01" & dttm_local < "2025-01-18") |>
  mutate(round_time = round_date(dttm_local, unit = "hour")) |> 
  group_by(round_time) |> 
  mutate(ph_mean = mean(pH)) |> 
  ungroup()

plot_outliers(test, "ph_mean")

plot <- plot_outliers(maestro_df_clean, "pH")
#plot

# pH has some outliers but they all seem to follow the trend. There is one instance of an acidic 
# pH which feels weird. The pH also goes pretty haywire at the end getting reall basic.
# I would assume that those aren't true but I will need to do more looking into it. It does seem like that 
# occures when pCO2 is highest which makes me think it could be a real-ish value. 

# Based on some literature, it looks like it would be unreasonable for it to be a real value.
# I am going to remove the data on 2024-12-08 which is a sampling date unfortunately. That one will have 
# to go without pH data. 
maestro_df_clean <- maestro_df_clean |> 
  mutate(pH = ifelse(date(dttm_local) == "2024-12-08", as.numeric("NA"), pH))

# Replot pH
plot <- plot_outliers(maestro_df_clean, "pH")
#plot
# Looks better with that day removed, now check out the end of the data

test <- maestro_df_clean |> 
  filter(date(dttm_local) > "2025-01-31") |>
  # Create an hourly average
    mutate(round_time = round_date(dttm_local, unit = "hour")) |> 
  group_by(round_time) |> 
  mutate(ph_mean = mean(pH)) |> 
  ungroup() |> 
  # Plot
  ggplot(aes(dttm_local, ph_mean)) +
  geom_line()
maestro_df_clean <- maestro_df_clean |>
  # Find the date after jan 1 where pH is first 8.6 (the previous max)
  # Then turn all pH after that date to "NA"
  mutate(pH = ifelse(dttm_local > maestro_df_clean$dttm_local[min(which(maestro_df_clean$dttm_local > "2025-01-01"
  & maestro_df_clean$pH > 8.6))], as.numeric('NA'), pH)) 
plot_outliers(maestro_df_clean, "pH")

# # Look at correlations
# test <- maestro_df_clean |> 
#   select(all_of(data_vars)) |> 
#   select(-c(CND, CDOM, TURB, T2, PHYC))
# colnames(test)
# GGally::ggpairs(test)
# test


plot_outliers(maestro_df_clean, "pH")
plot_outliers(maestro_df, "pH")
# pH looks better after removing that data at then end. It is hard to tell exactly where 
# things go wonky but
# I might just have to deal with it for now. 
# There is an interesting PCO2 and pH correlation

test <- maestro_df_clean |> 
  # filter(dttm_local > "2025-01-01" & dttm_local < "2025-05-15") |> 
  ggplot(aes(pH, T, color = dttm_local)) +
  geom_point() +
  scale_color_viridis_c()
plotly_build(test)

# It looks like it's not correlated with temperature at the end which seems like its not a temperature failue
# But I could possibly use the differing patterns to see when stuff changed in the future. For 
# now I will just leave it as is


## Recalculate multivariate outliers just for fun --------------------------------------
# #outliers <- maestro_df_clean |> 
#   select(data_vars) |>  
#   drop_na() |> 
#   check_outliers(method = "mcd", verbose = T)
# #outliers

# #maestro_df_clean |> 
#   select(data_vars, dttm) |> 
#   drop_na() |> 
#   mutate(outlier = outliers) |> 
#   pivot_longer(-c(dttm, outlier), names_to = "variable", values_to = "vals") |> 
#   ggplot(aes(dttm, vals, color = outlier)) +
#   geom_point() +
#   facet_wrap(~variable, scales = "free")


# Plot each variable looking specifically at the sampling dates
plot_samples <- function(maestro_df, var, plotly = T) {
  sample_df <- maestro_df |> 
    filter(date == sample_date) |> 
    select(dttm_local,  var, sample_date) |> 
    rename(data = var) |> 
    drop_na()
  # Create a new column
  sample_df <- sample_df |> 
    group_by(sample_date) |> 
    mutate(coeff_var = sd(data)/mean(data)*100) |>
    ungroup()

  # Make a plot of samples
  sample_plot <- 
  sample_df |> 
  ggplot(aes(dttm_local, data, color = coeff_var)) + 
  geom_point() +
  ggtitle(paste( var, "samples")) +
    facet_wrap(~sample_date, scales = "free")
  
  sample_plotly <- plotly_build(sample_plot)
    if(plotly) {
    return(sample_plotly)
    } else {
    return(sample_plot)
    }
  }

#for(i in 1:length(data_vars)) {
 # plot <- plot_samples(maestro_df_clean, data_vars[i], plotly = F)
#  print(plot)
#}

#maestro_df_clean_saved <- maestro_df_clean
 # Looking at all these plots it seems like turbidity is the only one with a crazy coefficient of variance.
 # I'm just going to take the mean value. 
## Average data ------------------------------------------------------------------------------
maestro_df_clean <- maestro_df_clean|>  
  group_by(date) |> 
  mutate(cdom_24h_avg = mean(CDOM, na.rm = TRUE),
         chla_24h_avg = mean(CHL, na.rm = TRUE),
         cond_24h_avg = mean(CND, na.rm = TRUE),
         do_24h_avg = mean(DO, na.rm = TRUE),
         do_sat_24h_avg = mean(do_sat_percent, na.rm = TRUE),
         pres_db_24h_avg = mean(P, na.rm = TRUE),
         pco2_24h_avg = mean(PCO2, na.rm = TRUE),
         phyc_24h_avg = mean(PHYC, na.rm = TRUE),
         mtemp_c_24h_avg = mean(T, na.rm = TRUE),
         mtemp_2_c_24h_avg = mean(T2, na.rm = TRUE),
         turb_24h_avg = mean(TURB, na.rm = TRUE),
         ph_24h_avg = mean(pH, na.rm = TRUE),
) |> 
  ungroup()

# Create noon and midnight mean values
maestro_df_clean <- maestro_df_clean|>  
  mutate(round_date = round_date(dttm_local, unit = "12 hours")) |> 
  group_by(round_date) |> 
  mutate(cdom_12h_med = median(CDOM, na.rm = TRUE),
         chla_12h_med = median(CHL, na.rm = TRUE),
         cond_12h_med = median(CND, na.rm = TRUE),
         do_12h_med = median(DO, na.rm = TRUE),
         do_sat_12h_med = median(do_sat_percent, na.rm = TRUE),
         pres_db_12h_med = median(P, na.rm = TRUE),
         pco2_12h_med = median(PCO2, na.rm = TRUE),
         phyc_12h_med = median(PHYC, na.rm = TRUE),
         mtemp_c_12h_med = median(T, na.rm = TRUE),
         mtemp_2_c_12h_med = median(T2, na.rm = TRUE),
         turb_12h_med = median(TURB, na.rm = TRUE),
         ph_12h_med = median(pH, na.rm = TRUE),
) |> 
  ungroup()


clean_long_m_df <- maestro_df_clean |> 
  select(all_of(data_vars), dttm, date, sample, sample_date) |> 
  pivot_longer(cols = -c(dttm, date, sample, sample_date),
  names_to = "variable", values_to = "vals", values_drop_na = T) |> 
  # Create a coefficient of variation to see if we should average or not 
  group_by(date, variable) |> 
  mutate(coeff_var = (sd(vals, na.rm = T)/mean(vals, na.rm = T)) * 100) |> 
  ungroup()


#   # Plot coeff_var only on sample dates
# clean_long_m_df |> 
#   mutate(coeff_test = ifelse(coeff_var > 20, "FALSE", "TRUE")) |> 
#   filter(sample_date == date, variable == "CHL") |> 
#   ggplot(aes(dttm, vals, color = coeff_var)) +
#   scale_color_viridis_c()+
#   geom_point() +
#   facet_wrap(~sample_date, scales = "free")


#### Add maestro data to broader metadata file --------------------------------------------------
sw_meta_output <- maestro_df_clean |> 
  # Filter so rounded date is equal to the sample date (the 24 hrs prior to 'sample_date' are 
  # averaged)
  filter(date == sample_date) |> 
  select(!any_of(c(data_vars, "dttm", "sample", "t", "dttm_local"))) |> 
  right_join(y = sw_meta_final, by = join_by(date)) |> 
# Relocate columns to make more sense
    select(sample, date, depth, strat_season, everything()) |> 
  # Put samples in alphabetical order
  arrange(sample) |> 
  distinct() |> 
  # Add a column for when maestro temp is within 1 deg of thermistor temp
  mutate(use_maestro = as.logical(ifelse(
    abs(mtemp_c_24h_avg - temp_c_24h_med) < 1, "TRUE", "FALSE")),
     .after = strat_season)

# maestro_df <- maestro_df |> 
#   left_join(avg_df, by = "date")
#   drop_na() |> 
#   check_outliers(ID = "dttm") |> 
#   which()
# maestro_plot <- maestro_df |>
#   slice(outliers) |> 
#   pivot_longer(cols = -dttm, names_to = "variable", values_to = "vals") |> 
#   ggplot(aes(x = dttm, y = vals, color = variable)) +
#   geom_point() +
#   facet_wrap(~variable, scales = "free")

# maestro_plotly <- plotly::plotly_build(maestro_plot) 

# maestro_plot
# maestro_plotly

# pres_plot <-
#   long_maestro_df |> 
#   filter(variable == "P") |> 
#   filter(vals < 12) |> 
#   ggplot(aes(x = dttm, y = vals, color = variable)) +
#   geom_point()

# pres_plot
# pres_plotly <- plotly::plotly_build(pres_plot)
# pres_plotly

stop()
# Write Files to output ----------------------------------------------------------------
sw_meta_output <- sw_meta_output |> 
  write_csv("output/data/metadata_supwinter.csv")
save(file = "output/data/maestro_df_clean.RData", maestro_df_clean)
save(file = "output/data/clean_long_m_df.RData", clean_long_m_df)
save(file = "output/data/temp_df_long.RData", temp_df_long)
# maestro_df_clean <- maestro_df_clean |> 
#   write_csv("output/data/maestro_data_clean.csv")
# # temp_df_long <- temp_df_long |> 
#   # write_csv("output/data/temp_data_long.csv")

obs <- ls()

rm(list = obs)
rm(obs)
