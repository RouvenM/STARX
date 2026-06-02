library(bigtime)
library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)
library(openxlsx)


# 1. download data 

months <- paste0("2024", sprintf("%02d", 1:12))

trips <- lapply(months, function(m) {
  url  <- paste0("https://s3.amazonaws.com/tripdata/JC-", m, "-citibike-tripdata.csv.zip")
  dest <- paste0("JC-", m, "-citibike-tripdata.csv.zip")
  if (!file.exists(dest))
    tryCatch(download.file(url, dest, mode = "wb", quiet = TRUE),
             error = function(e) cat("failed:", m, "\n"))
  if (file.exists(dest)) {
    f <- unzip(dest, exdir = ".")
    read.csv(f[1], stringsAsFactors = FALSE)
  }
}) %>% bind_rows()


# 2. top stations 

stations <- trips %>%
  filter(!is.na(start_station_id), start_station_id != "") %>%
  group_by(start_station_id, start_station_name) %>%
  summarise(trips = n(),
            lat   = mean(start_lat, na.rm = TRUE),
            lng   = mean(start_lng, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(desc(trips)) %>%
  slice(1:50)


# 3. build time series matrix 

trips$started_at <- as.POSIXct(trips$started_at, format = "%Y-%m-%d %H:%M:%S")
trips$period     <- floor_date(trips$started_at, unit = "30 minutes")

periods <- seq(
  floor_date(min(trips$period, na.rm = TRUE), "day"),
  floor_date(max(trips$period, na.rm = TRUE), "day") + days(1) - minutes(30),
  by = "30 min"
)

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

Y_raw <- as.matrix(ts_wide[, -1])
name_map        <- setNames(stations$start_station_name, stations$start_station_id)
colnames(Y_raw) <- substr(name_map[colnames(Y_raw)], 1, 20)


# 4. plot raw data 

plot_series <- function(Y, title = "") {
  as.data.frame(Y) %>%
    mutate(Time = 1:n()) %>%
    pivot_longer(-Time, names_to = "Series", values_to = "vals") %>%
    mutate(Series = factor(Series, levels = colnames(Y))) %>%
    ggplot() +
    geom_line(aes(Time, vals)) +
    facet_wrap(~Series, ncol = 1, scales = "free_y") +
    ylab("") + theme_bw() + ggtitle(title)
}

plot(plot_series(Y_raw[, 1:6], "Trips per half-hour — top 6 stations (2024)"))

ggplot(stations, aes(lng, lat)) +
  geom_point(aes(size = trips), color = "steelblue", alpha = 0.7) +
  geom_text(aes(label = substr(start_station_name, 1, 18)),
            size = 2.2, vjust = -0.8, check_overlap = TRUE) +
  theme_bw() +
  labs(title = "Citi Bike — Jersey City / Hoboken 2024",
       x = "longitude", y = "latitude", size = "trips")


# 5. initial fit (4 weeks) 

n_init  <- 4 * 7 * 48
VAR_fit <- sparseVAR(Y = scale(Y_raw[1:n_init, ]), VARpen = "HLag", selection = "cv")

plot_cv(VAR_fit)
Lhat <- lagmatrix(fit = VAR_fit, returnplot = TRUE)
plot(diagnostics_plot(VAR_fit, variable = colnames(Y_raw)[1]))
cat("is.stable:", is.stable(VAR_fit), "\n")


# 6. direct + recursive forecast 

n_train   <- n_init - 1
VAR_train <- sparseVAR(Y = scale(Y_raw[1:n_train, ]), VARpen = "HLag", selection = "cv")

directforecast(VAR_train, h = 1)

if (is.stable(VAR_train)) {
  rec <- recursiveforecast(VAR_train, h = 10)
  plot(rec, series = colnames(Y_raw)[1], last_n = 50)
}


# 7. rolling window evaluation

window  <- 4 * 7 * 48
idx     <- (window + 1):nrow(Y_raw)
k       <- ncol(Y_raw)

forecasts <- matrix(NA, nrow = length(idx), ncol = k, dimnames = list(NULL, colnames(Y_raw)))
actuals   <- matrix(NA, nrow = length(idx), ncol = k, dimnames = list(NULL, colnames(Y_raw)))

for (i in seq_along(idx)) {
  t     <- idx[i]
  Y_win <- Y_raw[(t - window):(t - 1), ]
  Y_sc  <- scale(Y_win)
  mu    <- attr(Y_sc, "scaled:center")
  sd    <- attr(Y_sc, "scaled:scale")
  sd[sd == 0] <- 1  # avoid division by zero

  fit  <- sparseVAR(Y = Y_sc, VARpen = "HLag", selection = "bic")
  fcst <- tryCatch(directforecast(fit, h = 1),
                   error = function(e) Y_win[nrow(Y_win), ])

  forecasts[i, ] <- as.numeric(fcst) * sd + mu
  actuals[i, ]   <- Y_raw[t, ]

  if (i %% 200 == 0) cat("progress:", i, "/", length(idx), "\n")
}


# 8. accuracy 

errors     <- forecasts - actuals
msfe_var   <- colMeans(errors^2)
naive_msfe <- colMeans((Y_raw[idx - 1, ] - actuals)^2)

comparison <- data.frame(
  station     = colnames(Y_raw),
  msfe_var    = round(msfe_var, 4),
  msfe_naive  = round(naive_msfe, 4),
  improvement = paste0(round((1 - msfe_var / naive_msfe) * 100, 1), "%")
) %>% arrange(msfe_var)

cat("\noverall VAR MSFE:  ", round(mean(msfe_var), 4),
    "\noverall naive MSFE:", round(mean(naive_msfe), 4),
    "\nimprovement:       ", round((1 - mean(msfe_var) / mean(naive_msfe)) * 100, 1), "%\n")


# 9. plots 

data.frame(station = colnames(Y_raw),
           msfe    = msfe_var,
           better  = msfe_var < naive_msfe) %>%
  arrange(msfe) %>%
  mutate(station = factor(station, levels = station)) %>%
  ggplot(aes(station, msfe, fill = better)) +
  geom_col() +
  geom_hline(yintercept = mean(naive_msfe), linetype = "dashed",
             color = "red", linewidth = 0.8) +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "#E57373"),
                    labels = c("TRUE" = "better than naive", "FALSE" = "worse than naive")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) +
  labs(title = "Forecast accuracy by station — 2024",
       x = "", y = "MSFE", fill = "")

data.frame(index    = idx,
           forecast = forecasts[, 1],
           actual   = actuals[, 1]) %>%
  slice(1:(14 * 48)) %>%
  pivot_longer(-index, names_to = "type", values_to = "value") %>%
  ggplot(aes(index, value, color = type)) +
  geom_line(linewidth = 0.4) +
  scale_color_manual(values = c("actual" = "black", "forecast" = "steelblue")) +
  theme_bw() +
  labs(title    = paste("forecast vs actual —", colnames(Y_raw)[1]),
       subtitle = "first 2 weeks of evaluation period",
       x = "half-hour index", y = "trips", color = "")
