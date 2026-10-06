library(bigtime)
library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)
library(openxlsx)


# 1. download data 

trips <- {
  url  <- "https://s3.amazonaws.com/tripdata/JC-202401-citibike-tripdata.csv.zip"
  dest <- "JC-202401-citibike-tripdata.csv.zip"
  if (!file.exists(dest))
    tryCatch(download.file(url, dest, mode = "wb", quiet = TRUE),
             error = function(e) cat("failed\n"))
  f <- unzip(dest, exdir = ".")
  read.csv(f[1], stringsAsFactors = FALSE)
}


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

periods <- seq(as.POSIXct("2024-01-01 00:00:00"),
               as.POSIXct("2024-01-31 23:00:00"), by = "1 hour")

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

p_series <- as.data.frame(Y_raw) %>%
  mutate(time = periods) %>%
  pivot_longer(-time, names_to = "station", values_to = "trips") %>%
  mutate(station = factor(station, levels = safe_names)) %>%
  ggplot(aes(time, trips)) +
  geom_line(color = "steelblue", linewidth = 0.3) +
  facet_wrap(~station, ncol = 1, scales = "free_y") +
  theme_bw() +
  labs(title = "Hourly trip counts — top 4 stations",
       subtitle = "Jersey City / Hoboken, January 2024", x = "", y = "Trips")
ggsave("plot_series.png", plot = p_series, width = 10, height = 8, dpi = 150)


# 4. initial fit 

window  <- 2 * 7 * 24
VAR_fit <- sparseVAR(Y = scale(Y_raw[1:window, ]), VARpen = "HLag", selection = "bic")
cat("is.stable:", is.stable(VAR_fit), "\n")

Lhat <- lagmatrix(fit = VAR_fit, returnplot = FALSE)
png("plot_lagmatrix.png", width = 800, height = 700, res = 150)
lagmatrix(fit = VAR_fit, returnplot = TRUE)
dev.off()

png("plot_diagnostics.png", width = 900, height = 600, res = 150)
print(diagnostics_plot(VAR_fit, variable = "area1"))
dev.off()


# 5. rolling window 

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

  fit  <- sparseVAR(Y = Y_sc, VARpen = "HLag", selection = "bic")
  fcst <- tryCatch(as.numeric(directforecast(fit, h = 1)),
                   error = function(e) Y_win[nrow(Y_win), ])

  forecasts[i, ] <- fcst
  actuals[i, ]   <- (Y_raw[t, ] - mu_w) / sd_w      # t not in Y_win
  naive_sc[i, ]  <- (Y_win[nrow(Y_win), ] - mu_w) / sd_w

  if (i %% 24 == 0) cat("progress:", i, "/", n_test, "\n")
}

cat("done!\n")


# 6. accuracy 

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


# 7. plots 

p_fcst <- data.frame(index    = test_idx,
                     Forecast = forecasts[, "area1"],
                     Actual   = actuals[, "area1"]) %>%
  pivot_longer(-index, names_to = "type", values_to = "value") %>%
  mutate(type = factor(type, levels = c("Actual", "Forecast"))) %>%
  ggplot(aes(index, value, color = type)) +
  geom_line(linewidth = 0.4) +
  scale_color_manual(values = c("Actual" = "black", "Forecast" = "steelblue")) +
  theme_bw() +
  labs(title    = "Forecast vs. Actual — area1",
       subtitle = "Rolling window | 1-step-ahead hourly | BIC",
       x = "Hour index", y = "Trips (standardized)", color = "")
ggsave("plot_forecast.png", plot = p_fcst, width = 12, height = 4, dpi = 150)
