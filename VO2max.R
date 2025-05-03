library(readxl)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(whippr)
library(lubridate)
library(segmented)
library(car)
library(tibble)

# data clean --------------------------------------------------------------
CPET_RAMP <- read_data(".xlsx", metabolic_cart = "nspire", time_column = "CPET Results")
CPET_RAMP <- CPET_RAMP |> rename_with(~ as.character(1:16), everything()) |> 
  dplyr::select(1,4,12) |> 
  rename(time = "1",
         vo2 = "4",
         po  = "12") |> 
  mutate(vo2 = as.numeric(as.character(vo2))) |> 
  slice(which(po == 0) [1]:n())
start_index <- which(CPET_RAMP$po == 50) [1]
end_index <- which(CPET_RAMP$po == 0 & seq_along(CPET_RAMP$po) > start_index) [1] -1
CPET_RAMP <- CPET_RAMP |> slice(start_index:end_index) 


CPET_RAMP <- CPET_RAMP |> 
  mutate(time = period_to_seconds(hms(time)))
intercept <- tail(CPET_RAMP$time[CPET_RAMP$po == 50],1)
CPET_RAMP$time <- CPET_RAMP$time - intercept #wanneeer het altijd 50 baseline waarden zijn zo doen, anders met seq along



CPET_RAMP <- CPET_RAMP |> rename("CPET Results" = time)

# normalize ramp test with whippr -----------------------------------------
ramp_normalized <- CPET_RAMP |> 
  incremental_normalize(
    .data = CPET_RAMP,
    incremental_type = "ramp",
    has_baseline = TRUE,
    baseline_length = 0, ## 4-min baseline
    work_rate_magic = TRUE,
    baseline_intensity = 50, ## baseline was performed at 50 W
    ramp_increase = 10 ## 10 W/min ramp
  )

# detect and remove outliers ----------------------------------------------
data_ramp_with_outliers <- 
  detect_outliers(
    .data = ramp_normalized,
    test_type = "incremental",
    vo2_column = "vo2",
    cleaning_level = 0.90,
    method_incremental = "anomaly", ## changed to anomaly detection
    verbose = TRUE) |>   
  plot_outliers() 

# Remove outliers
data_ramp_without_outliers <- data_ramp_with_outliers[["data"]]|> 
  dplyr::filter(outlier != "yes") 


data_ramp_without_outliers <- data_ramp_without_outliers|> rename("time" = `CPET Results`)


VO2max <-  mean(tail(data_ramp_without_outliers$vo2, 6), na.rm = TRUE)


# plot data ---------------------------------------------------------------
ggplot(data_ramp_without_outliers, aes(x = time, y = vo2 )) +
  geom_point(shape = 21, size = 3, fill = "white") + 
  labs(x = "time (s)", y = "vo2 (L/min)", title = "vo2/time") +
  theme_whippr() +
  scale_y_continuous(limits = c(0,5)) +
  scale_y_continuous(breaks = seq(0, 5, by = 0.5) )


table <-  tibble(
  Term = c("VO2max"),
  estimate = c(round(VO2max,5)))


