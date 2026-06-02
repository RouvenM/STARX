library(bigtime)
library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)
library(openxlsx)


# 1. download data (January 2024)

trips <- {
  url  <- "https://s3.amazonaws.com/tripdata/JC-202401-citibike-tripdata.csv.zip"
  dest <- "JC-202401-citibike-tripdata.csv.zip"
  if (!file.exists(dest))
    tryCatch(download.file(url, dest, mode = "wb", quiet = TRUE),
             error = function(e) cat("failed\n"))
  f <- unzip(dest, exdir = ".")
  read.csv(f[1], stringsAsFactors = FALSE)
}


# 2. top 20 stations 

stations <- trips %>%
  filter(!is.na(start_station_id), start_station_id != "") %>%
  group_by(start_station_id, start_station_name) %>%
  summarise(trips = n(),
            lat   = mean(start_lat, na.rm = TRUE),
            lng   = mean(start_lng, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(desc(trips)) %>%
  slice(1:20)


# 3. build time series matrix

trips$started_at <- as.POSIXct(trips$started_at, format = "%Y-%m-%d %H:%M:%S")
trips$period     <- floor_date(trips$started_at, unit = "30 minutes")

periods <- seq(
  as.POSIXct("2024-01-01 00:00:00"),
  as.POSIXct("2024-01-31 23:30:00"),
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

Y_raw      <- as.matrix(ts_wide[, -1])
name_map   <- setNames(stations$start_station_name, stations$start_station_id)
long_names <- substr(name_map[colnames(Y_raw)], 1, 20)  # for plots/output
safe_names <- paste0("S", seq_len(ncol(Y_raw)))          # for bigtime
colnames(Y_raw) <- safe_names

cat("matrix:", nrow(Y_raw), "half-hours x", ncol(Y_raw), "stations\n")


# 4. train/test split: 3 weeks train, 1 week test

n_train <- 3 * 7 * 48
n_test  <- 1 * 7 * 48

Y_train <- Y_raw[1:n_train, ]
Y_test  <- Y_raw[(n_train + 1):(n_train + n_test), ]

cat("train:", nrow(Y_train), "half-hours (3 weeks)\n")
cat("test: ", nrow(Y_test),  "half-hours (1 week)\n")


# 5. fit model 

Y_train_sc <- scale(Y_train)
mu <- attr(Y_train_sc, "scaled:center")
sd <- attr(Y_train_sc, "scaled:scale")
sd[sd == 0] <- 1

VAR_fit <- sparseVAR(Y = Y_train_sc, VARpen = "HLag", selection = "bic")

Lhat <- lagmatrix(fit = VAR_fit, returnplot = TRUE)
cat("is.stable:", is.stable(VAR_fit), "\n")


# 6. forecast over test week 

if (is.stable(VAR_fit)) {
  rec <- recursiveforecast(VAR_fit, h = n_test)
  plot(rec, series = safe_names[1], last_n = 50)
  # rec$fcst is 336 x 20 — backtransform columnwise
  fcst_matrix <- sweep(rec$fcst, 2, sd, "*")
  fcst_matrix <- sweep(fcst_matrix, 2, mu, "+")
} else {
  cat("model not stable — falling back to naive\n")
  fcst_matrix <- matrix(rep(Y_train[nrow(Y_train), ], n_test),
                        nrow = n_test, byrow = TRUE)
}


# 7. accuracy vs naive

msfe_var   <- colMeans((fcst_matrix - Y_test)^2)
naive_msfe <- colMeans((matrix(rep(Y_train[nrow(Y_train), ], n_test),
                               nrow = n_test, byrow = TRUE) - Y_test)^2)

comparison <- data.frame(
  station     = long_names,
  msfe_var    = round(msfe_var, 4),
  msfe_naive  = round(naive_msfe, 4),
  improvement = paste0(round((1 - msfe_var / naive_msfe) * 100, 1), "%")
) %>% arrange(msfe_var)

print(comparison)
cat("\noverall VAR MSFE:  ", round(mean(msfe_var), 4),
    "\noverall naive MSFE:", round(mean(naive_msfe), 4),
    "\nimprovement:       ", round((1 - mean(msfe_var) / mean(naive_msfe)) * 100, 1), "%\n")


# 8. plots

data.frame(station = long_names,
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
  labs(title = "forecast accuracy by station — test week (Jan 2024)",
       x = "", y = "MSFE", fill = "")

data.frame(index    = 1:n_test,
           forecast = fcst_matrix[, 1],
           actual   = Y_test[, 1]) %>%
  pivot_longer(-index, names_to = "type", values_to = "value") %>%
  ggplot(aes(index, value, color = type)) +
  geom_line(linewidth = 0.4) +
  scale_color_manual(values = c("actual" = "black", "forecast" = "steelblue")) +
  theme_bw() +
  labs(title    = paste("forecast vs actual —", long_names[1]),
       subtitle = "test week (week 4, January 2024)",
       x = "half-hour index", y = "trips", color = "")
