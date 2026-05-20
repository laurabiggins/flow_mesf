library(tidyverse)

wide_data <- read_csv("data/raw_data/mesf_beads_wide.csv")
tidy_flow <- wide_data |>
  pivot_longer(cols = `APC-Cy7`:SB780, names_to = "fluoro")
write_rds(tidy_flow, file = "mesf_tidied.rds")

