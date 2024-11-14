# data-raw/zheng_data.R

# Read the original data file
zheng_data <- readRDS("data-raw/Zheng10x_Immune+Jurkat.rda")

# Save it properly as package data
usethis::use_data(zheng_data, overwrite = TRUE)