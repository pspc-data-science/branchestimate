library(choroplethrMaps)
library(choroplethr)
data("df_pop_state")
data("df_pop_county")

usethis::use_data(df_pop_state, overwrite = TRUE)
usethis::use_data(df_pop_county, overwrite = TRUE)
