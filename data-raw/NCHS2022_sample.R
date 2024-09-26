## code to prepare `NCHS2022_sample` dataset goes here

save_wd <- "C:/Users/KIW58/OneDrive/Documents/GitHub/COMMA/data-raw/"
NCHS2022_sample <- read.csv(paste0(save_wd, "Nat2022US_sample_COMMA.csv"))

usethis::use_data(NCHS2022_sample, overwrite = TRUE)
