source("renv/activate.R")
# dcurver binary on aarch64-apple-darwin20 links against OpenMP (libiomp5),
# which is unavailable on macOS CI (fopenmp stripped). Force source install.
options(renv.config.binary.exclude = "dcurver")
