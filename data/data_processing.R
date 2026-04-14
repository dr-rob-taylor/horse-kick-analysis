library(tidyverse)

# ------------------------------------------------------------------------------
# Horse kick data: deaths per year, summed across 14 cavalry corps
# Source: von Bortkiewicz (1898)
# ------------------------------------------------------------------------------

data_path <- list.files( pattern = ".csv", recursive = T)
data_raw  <- read.csv( data_path, row.names = 1)
data_long <- data_raw |> 
  pivot_longer( everything(), names_to = "Corp", values_to = "Deaths" )

deaths_per_year <- rowSums( data_raw )
deaths_per_corp <- colSums( data_raw )

save(data_long, deaths_per_year, deaths_per_corp, file = "data/prussian_data.RData")
