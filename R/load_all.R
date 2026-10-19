# Load all required R and C++ functions 

library(Rcpp)
library(devtools)

base_url <- "https://raw.githubusercontent.com/RouvenM/STARX/main/"

r_files <- c(
  "R/diagnostics.R",
  "R/hvar.R", 
  "R/hvaralgorithms.R", 
  "R/plot_cv.R", 
  "R/simVAR.R"
)

message("Load R scripts...")
for (file in r_files) {
  source_url(paste0(base_url, file))
}

cpp_files <- c(
  "src/simVAR.cpp", 
  "src/hvar.cpp"
)

message("Compile C++ scripts...")
for (file in cpp_files) {
  tmp_cpp <- tempfile(fileext = ".cpp")
  download.file(paste0(base_url, file), tmp_cpp, quiet = TRUE)
  sourceCpp(tmp_cpp)
}

message("All functions loaded successfully.")