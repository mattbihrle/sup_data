# For loading .mat file
library(R.matlab)
library(tidyverse)
library(performance) # for outlier checks
library(skimr) # For summary stats
library(plotly) # For interactive plots
 
# Load metadata-------------------------------------------------------------------------
sw_meta <- read_csv("data/sw_metadata.csv") |>
  mutate(date = mdy(date)) |> 
  select(!day:year)
sw_meta$strat_season <- factor(sw_meta$strat_season, levels = c("summer", "fall", "winter", "spring"), ordered = T)

#import temp file ----------------------------------------------------------------------------
mat_df <- readMat.default(con = "data/metadata/LSEO24h_temp.mat")

# Extract dataframe we want from 'T' for temp and convert matlab time to poxit time
temp_df <- mat_df$T |> 
  as.data.frame()
# Set colnames and rownames
colnames(temp_df) <- mat_df$dep  
rownames(temp_df) <- mat_df$t 
# move dttm column out, convert to datetime and pivot long
temp_df_long <- temp_df |> 
  rownames_to_column(var = "dttm") |>
  mutate(dttm = as_datetime((as.numeric(dttm) - 719529) * 86400)) |> 
  pivot_longer(cols = !dttm, names_to = "depth", values_to = "temp") |> 
# Add averaged temperature column
  mutate(date = round_date(dttm, unit = "day")) |> 
  group_by(date, depth) |> 
  mutate(temp_c_24h_avg = mean(temp)) |> 
  ungroup()

## Create a 14 day averaged temperature---------------------------------
  # First create temporary df of sample dates
  samp_dates <- sw_meta |> 
    select(sample) |> 
    mutate(sample_date = as_date(sw_meta$date))
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

### Add averaged temp column to metadata ---------------------
sw_meta_final <- temp_df_long |> 
  filter(depth == 38) |> 
  select(date, depth, temp_c_24h_avg, temp_c_14d_avg) |> 
  distinct(date, .keep_all = T) |> 
  right_join(y = sw_meta, by = "date") |> 
  relocate(sample)

# Import maestro data --------------------------------------------------------------------
# Recreate a dttm column
maestro_df <- readMat.default(con = "data/metadata/lseo_maestro_time.mat") |> 
  data.frame() |> 
  mutate(dttm = as_datetime((as.numeric(t) - 719529) * 86400), .before = CDOM) |>
  # Create a date column rounded to nearest whole day for averaging later
  mutate(date = as_date(round_date(dttm, unit = "day"))) |> 
  # filter out the beginning of the deployment where pressure is less than 20m 
  filter(P > 20) |> 
  # add a sample date column for future averaging
  left_join(samp_dates, by = join_by(closest(date <= sample_date)))

  # Create summary stats to check for nas
sum_stats <- skim(maestro_df)
sum_stats

  # Create long df for plotting
long_maestro_df <- maestro_df |> 
  pivot_longer(cols = -c(dttm, date, sample, sample_date),
  names_to = "variable", values_to = "vals", values_drop_na = T)
  # Quick plot
long_maestro_df |> 
  ggplot(aes(dttm, vals, color = variable)) +
  geom_point()+
  facet_wrap(~variable, scales = "free")

# Create a vector of data variables to use later
data_vars <- maestro_df |> 
  select(-c(dttm, date, sample_date, sample, t)) |> 
  colnames() |> 
  as.character()
data_vars

## Start looking for outliers ---------------------------------------------------------
### Plot of outliers - multivariate ----------------------------------------------
# outliers <- maestro_df |> 
#   select(data_vars) |>  
#   drop_na() |> 
#   check_outliers(method = "mcd", verbose = T)
# outliers

# # Make a lil plot
# maestro_df |> 
#   mutate(outlier = outliers) |> 
#   pivot_longer(-c(dttm, outlier), names_to = "variable", values_to = "vals") |> 
#   ggplot(aes(dttm, vals, color = outlier)) +
#   geom_point() +
#   facet_wrap(~variable, scales = "free")

## Set outlier function ------------------------------------------------------------------

plot_outliers <- function(maestro_df, var, plotly = T) {
  outlier_df <- maestro_df |> 
    select(dttm,  var) |> 
    rename(data = var) |> 
    drop_na()
  # Create a new column
  outlier_df <- outlier_df |> 
    mutate(outlier = check_outliers(outlier_df, method = c("zscore", "iqr")))

  # Make a plot of outliers
  outlier_plot <- 
  outlier_df |> 
  ggplot(aes(dttm, data, color = outlier)) + 
  geom_point() +
  ggtitle(paste( var, "outliers"))
  
  outlier_plotly <- plotly_build(outlier_plot)
    if(plotly) {
    return(outlier_plotly)
    } else {
    return(outlier_plot)
    }
  }
### Outliers in CDOM ------------------------------------------------------------------

plot <- plot_outliers(maestro_df, "CDOM")
plot

# Based on the plot above the outliers at the beginning of the data seem to fit the trend. The outliers at the end however do not
# Remove all CDOM outliers at the end of the data

maestro_df_clean <- maestro_df |> 
  mutate(CDOM = ifelse(CDOM > 4.46e-02, CDOM, as.numeric("NA")))

#Try outliers plot again
plot <- plot_outliers(maestro_df_clean, "CDOM")
plot

# Looks good moving onto the next

### Outliers in CHL ------------------------------------------------------------------

plot <- plot_outliers(maestro_df, "CHL")
plot

# Based on the plot there are some obvious seeming outliers at the bottom of the data
# Removing those

maestro_df_clean <- maestro_df_clean |> 
  mutate(CHL = ifelse(CHL > 0.35, CHL, as.numeric("NA")))

#Try outliers plot again
plot <- plot_outliers(maestro_df_clean, "CHL")
plot

# Looks better moving on to the next
### Outliers in CND ------------------------------------------------------------------

plot <- plot_outliers(maestro_df, "CND")
plot
# Conductivity looks good apart from some gaps in the data, moving to the next

### Outliers in DO ------------------------------------------------------------------

plot <- plot_outliers(maestro_df, "DO")
plot

# DO looks like there is some craziness at the beginning but might just be diurnal cycles
# While the water column is stratified

# Plot just the beginning

maestro_df |> 
  filter(dttm < "2024-09-15") |> 
  mutate(morning = as.logical(ifelse(hour(dttm) < 5 | hour(dttm) > 17, "TRUE", "FALSE"))) |> 
  ggplot(aes(dttm, DO, color = morning)) +
  geom_point()

# Plot just a single day
maestro_df |> 
  filter(date(dttm) == "2024-09-16") |> 
  mutate(morning = as.logical(ifelse(hour(dttm) < 12, "TRUE", "FALSE"))) |> 
  ggplot(aes(dttm, DO, color = morning)) +
  geom_point()

# The single day does seem to be diurnal. I'll keep this data as is

### Outliers in P ------------------------------------------------------------------

plot <- plot_outliers(maestro_df, "P")
plot

# Pressure doesn't seem to have any crazy differences. I'm going to leave it as is

### Outliers in PCO2 ------------------------------------------------------------------

plot <- plot_outliers(maestro_df, "PCO2")
plot

# pco2 goes fricking crazzzzy at the end of the deployment, plot a single day to see wahat is up
maestro_df |> 
  filter(date(dttm) > "2025-05-17" & date(dttm) < "2025-06-17") |> 
  mutate(morning = as.logical(ifelse(hour(dttm) < 12, "TRUE", "FALSE"))) |> 
  ggplot(aes(dttm, PCO2, color = morning)) +
  geom_point()

maestro_df_clean <- maestro_df_clean |> 
  mutate(PCO2 = ifelse(PCO2 > 0.35, PCO2, as.numeric("NA")))

# PCO2 also seems to be diurnally driven. I'll keep this too

### Outliers in PHYC ------------------------------------------------------------------

plot <- plot_outliers(maestro_df, "PHYC")
plot

# Phycoerythrin goes way high and way low at the end of the deployment, going to zoom in on those
maestro_df |> 
  filter(date(dttm) > "2025-07-11") |> 
  mutate(morning = as.logical(ifelse(hour(dttm) < 12, "TRUE", "FALSE"))) |> 
  ggplot(aes(dttm, PHYC, color = morning)) +
  geom_point()

# It looks the the data is consistently bad after July 13, removing the data at the end

maestro_df_clean <- maestro_df_clean |> 
  mutate(PHYC = as.numeric(ifelse(date(dttm) >= "2025-07-13", "NA", PHYC)))

#Try outliers plot again
plot <- plot_outliers(maestro_df_clean, "PHYC")
plot

# Looks better moving on to the next

### Outliers in TURB ------------------------------------------------------------------

plot <- plot_outliers(maestro_df, "TURB")
plot

# Some pretty high turbidity levels in the middle but should keep those for now. Removing negative
# values
# Phycoerythrin goes way high and way low at the end of the deployment, going to zoom in on those
maestro_df |> 
  filter(date(dttm) > "2025-07-01") |> 
  mutate(morning = as.logical(ifelse(hour(dttm) < 12, "TRUE", "FALSE"))) |> 
  ggplot(aes(dttm, TURB, color = morning)) +
  geom_point()

 plot <- maestro_df |> 
  filter(date(dttm) > "2025-07-08") |> 
  mutate(morning = as.logical(ifelse(hour(dttm) < 12, "TRUE", "FALSE"))) |> 
  ggplot(aes(dttm, TURB, color = morning)) +
  geom_point()
 plot <- plotly_build(plot)
  plot
# It looks like the cutoff is around noon on july 09. Going to remove the data afterwards

maestro_df_clean <- maestro_df_clean |> 
  mutate(TURB = ifelse(dttm <= as_datetime("2025-07-09 12:00:00"), maestro_df$TURB, as.numeric("NA")))

#Try outliers plot again
plot <- plot_outliers(maestro_df_clean, "TURB")
plot

# Looks better moving on to the next

### Outliers in pH ------------------------------------------------------------------

plot <- plot_outliers(maestro_df, "pH")
plot

# pH has some outliers but they all seem to follow the trend. There is one instance of an acidic 
# pH which feels weird. The pH also goes pretty haywire at the end getting reall basic.
# I would assume that those aren't true but I will need to do more looking into it. 

# Leaving pH as is for now

## Recalculate multivariate outliers just for fun --------------------------------------
outliers <- maestro_df_clean |> 
  select(data_vars) |>  
  drop_na() |> 
  check_outliers(method = "mcd", verbose = T)
outliers

maestro_df_clean |> 
  select(data_vars, dttm) |> 
  drop_na() |> 
  mutate(outlier = outliers) |> 
  pivot_longer(-c(dttm, outlier), names_to = "variable", values_to = "vals") |> 
  ggplot(aes(dttm, vals, color = outlier)) +
  geom_point() +
  facet_wrap(~variable, scales = "free")

## Average data ------------------------------------------------------------------------------
maestro_df_clean <- maestro_df_clean |>  
  group_by(date) |> 
  mutate(cdom_24h_avg = mean(CDOM, na.rm = T),
         chla_24h_avg = mean(CND, na.rm = T),
         cond_24h_avg = mean(CND, na.rm = T),
         do_24h_avg = mean(DO, na.rm = T),
         pres_db_24h_avg = mean(P, na.rm = T),
         pco2_24h_avg = mean(PCO2, na.rm = T),
         phyc_24h_avg = mean(PHYC, na.rm = T),
         mtemp_c_24h_avg = mean(T, na.rm = T),
         mtemp_2_c_24h_avg = mean(T2, na.rm = T),
         turb_24h_avg = mean(TURB, na.rm = T),
         ph_24h_avg = mean(pH, na.rm = T),
) |> 
  ungroup()

#### Add maestro data to broader metadata file --------------------------------------------------
sw_meta_output <- maestro_df_clean |> 
  select(!any_of(c(data_vars, "dttm", "sample", "t"))) |> 
  right_join(y = sw_meta_final, by = join_by(date)) |> 
# Relocate columns to make more sense
    select(sample, date, depth, strat_season, everything()) |> 
  # Put samples in alphabetical order
  arrange(sample) |> 
  distinct() |> 
  # Add a column for when maestro temp is within 1 deg of thermistor temp
  mutate(use_maestro = as.logical(ifelse(
    abs(mtemp_c_24h_avg - temp_c_24h_avg) < 1, "TRUE", "FALSE")),
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


# Write Files to output ----------------------------------------------------------------
sw_meta_output <- sw_meta_output |> 
  write_csv("output/data/metadata_supwinter.csv")

# maestro_df_clean <- maestro_df_clean |> 
#   write_csv("output/data/maestro_data_clean.csv")
# # temp_df_long <- temp_df_long |> 
#   # write_csv("output/data/temp_data_long.csv")

rm(long_maestro_df, maestro_df, samp_dates, sum_stats, sw_meta, sw_meta_final,
temp_df, data_vars, mat_df, outliers, plot, plot_outliers)
