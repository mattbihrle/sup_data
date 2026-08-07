# Load required libraries
library(ggplot2)
library(dplyr)
library(scales)

# ==========================================
# 1. PHYSICAL CONSTANTS & SPM FORMULAS
# ==========================================
g <- 9.81              # Acceleration due to gravity (m/s^2)
kts_to_ms <- 0.514444  # Knots to m/s
nm_to_m   <- 1852      # Nautical miles to meters
m_to_ft   <- 3.28084   # Meters to feet

# Fetch-limited Wave Height (ft)
calc_H_fetch <- function(U_ms, X_m) {
  X_tilde <- (g * X_m) / (U_ms^2)
  H_tilde <- 0.283 * tanh(0.0125 * (X_tilde^0.42))
  H_m <- (H_tilde * U_ms^2) / g
  return(H_m * m_to_ft)
}

# Required Duration (hours) for Fetch X_m
calc_t_fetch <- function(U_ms, X_m) {
  X_tilde <- (g * X_m) / (U_ms^2)
  t_tilde <- 68.8 * (X_tilde^0.67)
  t_sec <- (t_tilde * U_ms) / g
  return(t_sec / 3600)
}

# Equivalent Fetch (m) for Duration t_hr
calc_X_eq <- function(U_ms, t_hr) {
  t_sec <- t_hr * 3600
  t_tilde <- (g * t_sec) / U_ms
  X_tilde <- (t_tilde / 68.8)^(1 / 0.67)
  X_m <- (X_tilde * U_ms^2) / g
  return(X_m)
}

# ==========================================
# 2. GENERATE DATASETS
# ==========================================
u_range_kts <- seq(5, 75, length.out = 150)

# A. Wind Speed Curves (Solid Black Lines)
wind_speeds_kts <- seq(10, 70, by = 5)
durations_hr <- exp(seq(log(0.1), log(72), length.out = 300))

wind_data <- do.call(rbind, lapply(wind_speeds_kts, function(u_kts) {
  u_ms <- u_kts * kts_to_ms
  data.frame(
    Duration_hr = durations_hr,
    Wind_kts = as.factor(u_kts),
    Wind_num = u_kts,
    Wave_ft = sapply(durations_hr, function(t) {
      X_eq <- calc_X_eq(u_ms, t)
      X_m <- min(X_eq, 500 * nm_to_m)
      calc_H_fetch(u_ms, X_m)
    })
  )
}))

# B. Standard Fetch Contours (Dashed Gray Lines)
fetches_nm <- c(1, 2, 5, 10, 20, 30, 50, 75, 100, 150, 200, 300, 400, 500)

fetch_data <- do.call(rbind, lapply(fetches_nm, function(x_nm) {
  x_m <- x_nm * nm_to_m
  data.frame(
    Fetch_nm = as.factor(x_nm),
    Fetch_num = x_nm,
    Duration_hr = sapply(u_range_kts, function(u) calc_t_fetch(u * kts_to_ms, x_m)),
    Wave_ft = sapply(u_range_kts, function(u) calc_H_fetch(u * kts_to_ms, x_m))
  )
})) %>% filter(Duration_hr >= 0.1 & Duration_hr <= 72 & Wave_ft >= 0.3 & Wave_ft <= 60)

# C. Lake Superior Overall Max Limits (Red Dotted Lines)
lake_superior_fetches <- c(139, 304)
lake_superior_names <- c("Lake Superior Max Width (139 NM)", "Lake Superior Max Length (304 NM)")

lake_data <- do.call(rbind, lapply(1:2, function(i) {
  x_nm <- lake_superior_fetches[i]
  x_m <- x_nm * nm_to_m
  data.frame(
    Label = lake_superior_names[i],
    Fetch_nm = x_nm,
    Duration_hr = sapply(u_range_kts, function(u) calc_t_fetch(u * kts_to_ms, x_m)),
    Wave_ft = sapply(u_range_kts, function(u) calc_H_fetch(u * kts_to_ms, x_m))
  )
})) %>% filter(Duration_hr >= 0.1 & Duration_hr <= 72 & Wave_ft >= 0.3 & Wave_ft <= 60)

# D. Duluth Entry Rays to Land (Green Dotted Lines - Only >= 1 NM)
duluth_directions <- c("ENE", "E", "ESE", "NE", "NNE", "SE")
duluth_fetches_nm <- c(225, 21.5, 10.0, 8.6, 3.6, 3.5)

duluth_data <- do.call(rbind, lapply(1:length(duluth_directions), function(i) {
  x_nm <- duluth_fetches_nm[i]
  x_m <- x_nm * nm_to_m
  data.frame(
    Direction = duluth_directions[i],
    Fetch_nm = x_nm,
    Duration_hr = sapply(u_range_kts, function(u) calc_t_fetch(u * kts_to_ms, x_m)),
    Wave_ft = sapply(u_range_kts, function(u) calc_H_fetch(u * kts_to_ms, x_m))
  )
})) %>% filter(Duration_hr >= 0.1 & Duration_hr <= 72 & Wave_ft >= 0.3 & Wave_ft <= 60)

# E. Wave Period Contours (Dotted Gray Lines)
periods_sec <- c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)
period_data <- do.call(rbind, lapply(periods_sec, function(t_sec) {
  u_vals <- seq(5, 75, length.out = 150)
  res <- data.frame()
  for (u in u_vals) {
    u_ms <- u * kts_to_ms
    T_tilde <- (g * t_sec) / u_ms
    arg <- T_tilde / (2.4 * pi)
    if (arg < 1) {
      X_tilde <- (atanh(arg) / 0.077)^4
      X_m <- min((X_tilde * u_ms^2) / g, 500 * nm_to_m)
      dur <- calc_t_fetch(u_ms, X_m)
      h_ft <- calc_H_fetch(u_ms, X_m)
      res <- rbind(res, data.frame(Period_sec = as.factor(t_sec), Duration_hr = dur, Wave_ft = h_ft))
    }
  }
  return(res)
})) %>% filter(Duration_hr >= 0.1 & Duration_hr <= 72 & Wave_ft >= 0.3 & Wave_ft <= 60)

# Line Label Data
wind_labels <- wind_data %>%
  group_by(Wind_kts) %>%
  filter(Duration_hr == max(Duration_hr))

fetch_labels <- fetch_data %>%
  group_by(Fetch_nm) %>%
  filter(Wave_ft == max(Wave_ft))

# ==========================================
# 3. BUILD PLOT
# ==========================================
p <- ggplot() +
  # Wave Period Contours
  geom_line(data = period_data, aes(x = Duration_hr, y = Wave_ft, group = Period_sec),
            linetype = "dotted", color = "gray50", linewidth = 0.35) +
  
  # Standard Fetch Contours
  geom_line(data = fetch_data, aes(x = Duration_hr, y = Wave_ft, group = Fetch_nm),
            linetype = "dashed", color = "gray30", linewidth = 0.45) +
  
  # Lake Superior Overall Limits (Red Dotted Lines)
  geom_line(data = lake_data, aes(x = Duration_hr, y = Wave_ft, group = Label),
            linetype = "dotted", color = "red", linewidth = 0.85) +
  
  # Duluth Entry Ray Fetches (Green Dotted Lines)
  geom_line(data = duluth_data, aes(x = Duration_hr, y = Wave_ft, group = Direction),
            linetype = "dotted", color = "forestgreen", linewidth = 0.85) +
  
  # Wind Speed Curves
  geom_line(data = wind_data, aes(x = Duration_hr, y = Wave_ft, group = Wind_kts),
            color = "black", linewidth = 0.75) +
  
  # Contour Labels
  geom_text(data = fetch_labels, aes(x = Duration_hr, y = Wave_ft, label = paste0(Fetch_nm, " NM")),
            vjust = -0.4, hjust = 1.05, size = 2.5, fontface = "italic", color = "gray20") +
  geom_text(data = wind_labels, aes(x = Duration_hr, y = Wave_ft, label = paste0(Wind_kts, " kt")),
            hjust = -0.15, size = 2.7, fontface = "bold") +
  
  # --- TOP LEGEND BOX: DULUTH ENTRY FETCH TO LAND (GREEN) ---
  annotate("rect", xmin = 31, xmax = 82, ymin = 0.72, ymax = 1.85,
           fill = "white", color = "forestgreen", linetype = "solid", linewidth = 0.35, alpha = 0.95) +
  annotate("text", x = 33, y = 1.68, label = "Green Dotted Lines:",
           color = "forestgreen", fontface = "bold", hjust = 0, size = 2.2) +
  annotate("text", x = 33, y = 1.18,
           label = "Duluth Entry Ray Fetch to Land:\n• NNE (22.5°): 3.6 NM (Lester Pt)\n• NE (45.0°): 8.6 NM (Stony Pt)\n• ENE (67.5°): 225 NM (Ontario)\n• E (90.0°): 21.5 NM (Port Wing)\n• ESE (112.5°): 10 NM (Amnicon R)\n• SE (135.0°): 3.5 NM (WI Point)",
           color = "forestgreen", hjust = 0, size = 1.85) +
  
  # --- BOTTOM LEGEND BOX: LAKE SUPERIOR MAX LIMITS (RED) ---
  annotate("rect", xmin = 31, xmax = 82, ymin = 0.315, ymax = 0.68,
           fill = "white", color = "red", linetype = "solid", linewidth = 0.35, alpha = 0.95) +
  annotate("text", x = 33, y = 0.56, label = "Red Dotted Lines:",
           color = "red", fontface = "bold", hjust = 0, size = 2.2) +
  annotate("text", x = 33, y = 0.40, label = "• Superior Length (304 NM)\n• Superior Width (139 NM)",
           color = "red", hjust = 0, size = 1.9) +
  
  # Axes & Log Scale Limits
  scale_x_log10(
    breaks = c(0.1, 0.25, 0.5, 1, 2, 3, 5, 8, 12, 18, 24, 36, 48, 72),
    labels = c("0.1", "0.25", "0.5", "1", "2", "3", "5", "8", "12", "18", "24", "36", "48", "72"),
    limits = c(0.1, 85)
  ) +
  scale_y_log10(
    breaks = c(0.3, 0.5, 1, 2, 3, 5, 8, 12, 16, 20, 25, 30, 40, 50, 60),
    labels = c("0.3", "0.5", "1", "2", "3", "5", "8", "12", "16", "20", "25", "30", "40", "50", "60"),
    limits = c(0.3, 65)
  ) +
  
  # Titles & Theme Formatting
  labs(
    title = "DEEP WATER WAVE FORECASTING CHART",
    subtitle = "Duluth Entry Ray Fetches to Land vs. Lake Superior Limits",
    x = "Duration (hours)",
    y = "Significant Wave Height (feet)"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle = element_text(hjust = 0.5, size = 9.5, color = "gray20"),
    panel.grid.major = element_line(color = "gray80", linewidth = 0.3),
    panel.grid.minor = element_line(color = "gray92", linewidth = 0.2),
    axis.title = element_text(face = "bold"),
    plot.margin = margin(15, 15, 15, 15)
  )
p
stop()

# ==========================================
# 4. EXPORT TO FILE
# ==========================================
ggsave("wave_fetch_superior.pdf", plot = p, width = 8.5, height = 11, units = "in")
ggsave("wave_fetch_superior.png", plot = p, width = 8.5, height = 11, units = "in", dpi = 300)
