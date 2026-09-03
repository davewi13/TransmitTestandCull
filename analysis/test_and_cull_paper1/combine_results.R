library(dplyr)
library(purrr)

# Path to results
res_path <- "../Results_onlyR0_quarterlatent"

# List all chunk files
files <- list.files(res_path, pattern = "^chunk\\d+\\.rds$", full.names = TRUE)

# Read them all in and row-bind
all_results <- files %>%
  map(readRDS) %>%
  bind_rows()

# Save combined results
saveRDS(all_results, file = file.path(res_path, "all_results_onlyR0_quarterlatent.rds"))
