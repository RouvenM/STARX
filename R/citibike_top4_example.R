library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)
library(openxlsx)


# 1. download data 

trips <- {
  url  <- "https://s3.amazonaws.com/tripdata/JC-202403-citibike-tripdata.csv.zip"
  dest <- "JC-202403-citibike-tripdata.csv.zip"
  if (!file.exists(dest))
    tryCatch(download.file(url, dest, mode = "wb", quiet = TRUE),
             error = function(e) cat("failed\n"))
  f <- unzip(dest, exdir = ".")
  read.csv(f[1], stringsAsFactors = FALSE)
}

# 1.5 Rough Map of Bike Stations

library(ggrepel)
library(ggspatial)

# 1. Summarize station locations and trip counts from raw data
station_summary <- trips %>%
  filter(!is.na(start_station_id), start_station_id != "") %>%
  group_by(start_station_id) %>%
  summarise(
    total_trips = n(),
    lat         = mean(start_lat, na.rm = TRUE),
    lng         = mean(start_lng, na.rm = TRUE),
    .groups     = "drop"
  )

# Version A: Minimalist Plot (No background map)
ggplot(station_summary, aes(x = lng, y = lat)) +
  # Station points scaled by trip volume
  geom_point(
    aes(size = total_trips),
    color = "#2b8cbe",
    alpha = 0.7
  ) +
  # Station ID labels
  geom_text_repel(
    aes(label = start_station_id),
    size = 2.8,
    max.overlaps = 50,
    box.padding = 0.2,
    point.padding = 0.1
  ) +
  scale_size_continuous(
    name = "Total Trips",
    labels = scales::label_comma()
  ) +
  coord_quickmap() +
  theme_minimal() +
  labs(
    title = "Citibike Stations Map",
    subtitle = "All stations with Station IDs",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom"
  )

# Version B: Plot with Basemap Background (OpenStreetMap)
ggplot() +
  # Render OpenStreetMap tile layer
  annotation_map_tile(type = "osm", zoom = 12) +
  # Station points
  geom_point(
    data = station_summary,
    aes(x = lng, y = lat, size = total_trips),
    color = "#2b8cbe",
    alpha = 0.85
  ) +
  # Labels with semi-transparent background box for readability over map features
  geom_label_repel(
    data = station_summary,
    aes(x = lng, y = lat, label = start_station_id),
    size = 2.5,
    fontface = "bold",
    fill = alpha("white", 0.75),
    label.padding = unit(0.1, "lines"),
    box.padding = 0.2,
    point.padding = 0.1,
    max.overlaps = 50
  ) +
  coord_sf(crs = 4326) +
  scale_size_continuous(
    name = "Total Trips",
    labels = scales::label_comma()
  ) +
  theme_minimal() +
  labs(
    title = "Citibike Stations Map",
    subtitle = "All stations with Station IDs (OSM)",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom"
  )

# 2. top 4 stations

stations <- trips %>%
  filter(!is.na(start_station_id), start_station_id != "") %>%
  group_by(start_station_id, start_station_name) %>%
  summarise(trips = n(),
            lat   = mean(start_lat, na.rm = TRUE),
            lng   = mean(start_lng, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(desc(trips)) %>%
  slice(1:4)


# 3. time series matrix (hourly)

trips$started_at <- as.POSIXct(trips$started_at, format = "%Y-%m-%d %H:%M:%S")
trips$period     <- floor_date(trips$started_at, unit = "1 hour")

periods <- seq(as.POSIXct("2024-03-01 00:00:00"),
               as.POSIXct("2024-03-31 23:00:00"), by = "1 hour")

counts <- trips %>%
  filter(start_station_id %in% stations$start_station_id) %>%
  group_by(period, start_station_id) %>%
  summarise(trips = n(), .groups = "drop")

ts_wide <- expand.grid(period           = periods,
                       start_station_id = stations$start_station_id,
                       stringsAsFactors = FALSE) %>%
  left_join(counts, by = c("period", "start_station_id")) %>%
  mutate(trips = replace_na(trips, 0)) %>%
  pivot_wider(names_from = start_station_id, values_from = trips) %>%
  arrange(period)

Y_raw           <- as.matrix(ts_wide[, -1])
colnames(Y_raw) <- paste0("area", 1:4)
safe_names      <- colnames(Y_raw)

cat("matrix:", nrow(Y_raw), "hours x", ncol(Y_raw), "stations\n")

as.data.frame(Y_raw) %>% 
  mutate(time = periods) %>%
  pivot_longer(-time, names_to = "station", values_to = "trips") %>%
  mutate(station = factor(station, levels = safe_names)) %>%
  ggplot(aes(time, trips)) +
  geom_line(color = "steelblue", linewidth = 0.3) +
  facet_wrap(~station, ncol = 1, scales = "free_y") +
  theme_bw() +
  labs(title = "Hourly trip counts — top 4 stations",
       subtitle = "Jersey City / Hoboken, MONTH YEAR", x = "", y = "Trips")

# 3.5 fit on all observations
fit_bike <- sparseVAR(Y = scale(Y_raw), VARpen = "LagFused", selection = "cv")

plot_cv(fit_bike)

fit_bike$Phihat
plot_coefficient_path(fit_bike)
plot_row_frobenius_path(fit_bike)
plot_frobenius_path(fit_bike)

# 4. Alternative: time series matrix (4-hourly aggregation)

# Round trip start timestamps down to 4-hour intervals (00:00, 04:00, 08:00, ...)
trips$period_4h <- floor_date(trips$started_at, unit = "4 hours")

# Full 4-hourly sequence across the entire month
periods_4h <- seq(as.POSIXct("2024-03-01 00:00:00"),
                  as.POSIXct("2024-03-31 20:00:00"), by = "4 hours")

counts_4h <- trips %>%
  filter(start_station_id %in% stations$start_station_id) %>%
  group_by(period_4h, start_station_id) %>%
  summarise(trips = n(), .groups = "drop")

ts_wide_4h <- expand.grid(period_4h        = periods_4h,
                          start_station_id = stations$start_station_id,
                          stringsAsFactors = FALSE) %>%
  left_join(counts_4h, by = c("period_4h" = "period_4h", "start_station_id" = "start_station_id")) %>%
  mutate(trips = replace_na(trips, 0)) %>%
  pivot_wider(names_from = start_station_id, values_from = trips) %>%
  arrange(period_4h)

Y_raw_4h           <- as.matrix(ts_wide_4h[, -1])
colnames(Y_raw_4h) <- paste0("area", 1:4)

cat("4-hourly matrix:", nrow(Y_raw_4h), "periods x", ncol(Y_raw_4h), "stations\n")

# Plot 4-hourly time series
as.data.frame(Y_raw_4h) %>% 
  mutate(time = periods_4h) %>%
  pivot_longer(-time, names_to = "station", values_to = "trips") %>%
  mutate(station = factor(station, levels = safe_names)) %>%
  ggplot(aes(time, trips)) +
  geom_line(color = "darkgreen", linewidth = 0.5) +
  facet_wrap(~station, ncol = 1, scales = "free_y") +
  theme_bw() +
  labs(title = "4-Hourly trip counts — top 4 stations",
       subtitle = "MONTH YEAR", x = "", y = "Trips")


# 4.5 fit on 4-hourly observations

fit_bike_4h <- sparseVAR(p=6, Y = scale(Y_raw_4h), VARpen = "LagFused", selection = "cv", VARgran = c(500, 30))

plot_cv(fit_bike_4h)

fit_bike_4h$Phihat
plot_coefficient_path(fit_bike_4h)
plot_row_frobenius_path(fit_bike_4h)
plot_frobenius_path(fit_bike_4h)

fit_bike_largelambda <- sparseVAR(p=6, Y = scale(Y_raw_4h), VARpen = "LagFused", VARlseq = 10000, eps = 1e-10)
fit_bike_largelambda$Phihat[, 1:4, 1] - fit_bike_largelambda$Phihat[, 5:8, 1]
fit_bike_largelambda$Phihat[, 5:8, 1] - fit_bike_largelambda$Phihat[, 9:12, 1]
fit_bike_largelambda$Phihat[, 9:12, 1] - fit_bike_largelambda$Phihat[, 13:16, 1]
fit_bike_largelambda$Phihat[, 13:16, 1] - fit_bike_largelambda$Phihat[, 17:20, 1]
fit_bike_largelambda$Phihat[, 17:20, 1] - fit_bike_largelambda$Phihat[, 21:24, 1]


# 5. initial fit 

window  <- 2 * 7 * 24
VAR_fit <- sparseVAR(Y = scale(Y_raw[1:window, ]), VARpen = "LagFused", selection = "cv")
cat("is.stable:", is.stable(VAR_fit), "\n")

Lhat <- lagmatrix(fit = VAR_fit, returnplot = FALSE)
lagmatrix(fit = VAR_fit, returnplot = TRUE)

print(diagnostics_plot(VAR_fit, variable = "area1"))


# 6. rolling window 

test_idx  <- (window + 1):nrow(Y_raw)
n_test    <- length(test_idx)
k         <- ncol(Y_raw)
forecasts <- matrix(NA, nrow = n_test, ncol = k, dimnames = list(NULL, safe_names))
actuals   <- matrix(NA, nrow = n_test, ncol = k, dimnames = list(NULL, safe_names))
naive_sc  <- matrix(NA, nrow = n_test, ncol = k, dimnames = list(NULL, safe_names))

for (i in seq_along(test_idx)) {
  t     <- test_idx[i]
  Y_win <- Y_raw[(t - window):(t - 1), ]  # rows t-336 to t-1 only
  Y_sc  <- scale(Y_win)
  mu_w  <- attr(Y_sc, "scaled:center")
  sd_w  <- attr(Y_sc, "scaled:scale")
  sd_w[sd_w == 0] <- 1
  
  fit  <- sparseVAR(Y = Y_sc, VARpen = "LagFused", selection = "cv")
  fcst <- tryCatch(as.numeric(directforecast(fit, h = 1)),
                   error = function(e) Y_win[nrow(Y_win), ])
  
  forecasts[i, ] <- fcst
  actuals[i, ]   <- (Y_raw[t, ] - mu_w) / sd_w      # t not in Y_win
  naive_sc[i, ]  <- (Y_win[nrow(Y_win), ] - mu_w) / sd_w
  
  if (i %% 24 == 0) cat("progress:", i, "/", n_test, "\n")
}

cat("done!\n")


# 7. accuracy 

errors     <- forecasts - actuals
errors_sq  <- errors^2
errors_abs <- abs(errors)
msfe_var   <- colMeans(errors_sq)
naive_msfe <- colMeans((naive_sc - actuals)^2)

comparison <- data.frame(
  station     = safe_names,
  msfe_var    = round(msfe_var, 4),
  msfe_naive  = round(naive_msfe, 4),
  improvement = paste0(round((1 - msfe_var / naive_msfe) * 100, 1), "%")
) %>% arrange(msfe_var)

print(comparison)
cat("\noverall VAR MSFE:  ", round(mean(msfe_var), 4),
    "\noverall naive MSFE:", round(mean(naive_msfe), 4),
    "\nimprovement:       ", round((1 - mean(msfe_var) / mean(naive_msfe)) * 100, 1), "%\n")

results_per_window <- data.frame(window_idx = seq_along(test_idx),
                                 hour_idx   = test_idx,
                                 timestamp  = periods[test_idx])
for (s in safe_names) {
  results_per_window[[paste0("forecast_",  s)]] <- forecasts[, s]
  results_per_window[[paste0("actual_",    s)]] <- actuals[, s]
  results_per_window[[paste0("error_",     s)]] <- errors[, s]
  results_per_window[[paste0("sq_error_",  s)]] <- errors_sq[, s]
  results_per_window[[paste0("abs_error_", s)]] <- errors_abs[, s]
}

cumulative <- data.frame(window_idx = seq_along(test_idx))
for (s in safe_names) {
  cumulative[[paste0("cum_msfe_", s)]] <- cumsum(errors_sq[, s]) / seq_along(test_idx)
  cumulative[[paste0("cum_mae_",  s)]] <- cumsum(errors_abs[, s]) / seq_along(test_idx)
}


# 8. plots 

data.frame(index = test_idx, 
           Forecast = forecasts[, "area1"], 
           Actual   = actuals[, "area1"]) %>% 
  pivot_longer(-index, names_to = "type", values_to = "value") %>%
  mutate(type = factor(type, levels = c("Actual", "Forecast"))) %>%
  ggplot(aes(index, value, color = type)) +
  geom_line(linewidth = 0.4) +
  scale_color_manual(values = c("Actual" = "black", "Forecast" = "steelblue")) +
  theme_bw() +
  labs(title    = "Forecast vs. Actual — area1",
       subtitle = "Rolling window | 1-step-ahead hourly | CV",
       x = "Hour index", y = "Trips (standardized)", color = "")

