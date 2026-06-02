# STARX
Spatio-temporal autogression with exogenous variables

Sparse VAR baseline on Citi Bike station-level demand data (Jersey City / Hoboken).  
This repository documents the data pipeline, model setup, and forecast evaluation as a foundation for the STARX extension with exogenous variables.

---

## Data

**Source:** [Citi Bike System Data](https://s3.amazonaws.com/tripdata/index.html)  
**Coverage:** Jersey City / Hoboken, January 2024  
**Raw trips:** ~50,600 individual rides  
**Stations:** Top 20 by trip volume  
**Resolution:** 30-minute intervals → 1,488 half-hours × 20 stations

Data files are not included in this repository. The script downloads them automatically on first run.

---

## Method

### Time Series Construction

Each ride is assigned to its 30-minute interval via `floor_date()`. Trip counts are aggregated per station and interval; missing combinations (zero rides) are filled with 0 to produce a complete rectangular matrix.

### Model

Sparse VAR with hierarchical lag penalty (HLag) from the [`bigtime`](https://github.com/ineswilms/bigtime) package:

```
Y_t = A_1 Y_{t-1} + A_2 Y_{t-2} + ... + A_p Y_{t-p} + ε_t
```

- Penalty: `HLag` — encourages whole lags to drop out before individual coefficients
- Selection: BIC
- Standardization: `scale()` fit on training data only, applied before model estimation

### Train / Test Split

| Set | Period | Half-hours |
|-----|--------|-----------|
| Train | Weeks 1–3 (Jan 1–21) | 1,008 |
| Test | Week 4 (Jan 22–31) | 336 |

Forecast method: `recursiveforecast(h = 336)` — one model fit, rolled forward over the full test week.

### Baseline

Naive forecast: last observed value carried forward.

---

## Results

**Overall improvement over naive: 25.3%**  
VAR beats naive on 19 out of 20 stations.

| Station | MSFE (VAR) | MSFE (Naive) | Improvement |
|---------|-----------|-------------|-------------|
| Newport Pkwy | 0.71 | 0.98 | 27.9% |
| Marshall St & 2 St | 0.80 | 1.08 | 25.5% |
| Adams St & 2 St | 0.81 | 1.10 | 26.4% |
| Hoboken Ave at Monmouth | 0.83 | 1.10 | 25.1% |
| Clinton St & 7 St | 0.86 | 1.22 | 29.5% |
| Marin Light Rail | 0.87 | 1.21 | 28.1% |
| 11 St & Washington St | 1.08 | 1.20 | 9.8% |
| Columbus Park - Clinton | 1.12 | 1.61 | 30.6% |
| Madison St & 1 St | 1.16 | 1.60 | 27.1% |
| 8 St & Washington St | 1.32 | 1.33 | 1.0% |
| Hamilton Park | 1.45 | 1.96 | 26.4% |
| South Waterfront Walkway | 1.52 | 1.94 | 21.8% |
| City Hall - Washington St | 1.62 | 2.27 | 28.5% |
| Exchange Pl | 1.67 | 2.18 | 23.4% |
| Bergen Ave & Sip Ave | 2.09 | 2.98 | 29.9% |
| Newport PATH | 2.09 | 2.85 | 26.5% |
| River St & 1 St | 2.47 | 3.21 | 23.0% |
| Hoboken Terminal - Hudson St | 2.73 | 2.72 | -0.5% |
| Grove St PATH | 8.07 | 11.01 | 26.7% |
| Hoboken Terminal - River St | 8.45 | 12.32 | 31.3% |

MSFE values are in original scale (trips²). The two highest-MSFE stations are the busiest hubs in the network; their absolute error is larger but relative improvement is consistent with the rest.

---

## Repository Structure

```
STARX/
├── README.md
├── Literature/
└── R/
    └── citibike_sparseVAR_split.R
```

---

## Dependencies

```r
install.packages(c("bigtime", "ggplot2", "dplyr", "tidyr", "lubridate", "openxlsx"))
```

R version used: 4.4.1


