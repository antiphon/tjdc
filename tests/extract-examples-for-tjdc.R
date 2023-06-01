# synthetic
library(stars)
library(dplyr)

tfi <- "~/work/Projektit/MORE/manuscript-1/codes/data/case9.2.z1.tif"

c92 <- read_stars(tfi)
c92i <- readRDS(gsub(".tif", "_info.rds", tfi))

# make proper stars:
names(c92)  <-  "values"

c92b <- st_set_dimensions(c92, "band", c92i$timesteps, names = "z")

stars::st_get_dimension_values(c92b, "z")

plot(c92b)

case_medium_zoom1         <- c92b
case_medium_zoom1_details <- c92i

usethis::use_data(case_medium_zoom1, case_medium_zoom1_details)
