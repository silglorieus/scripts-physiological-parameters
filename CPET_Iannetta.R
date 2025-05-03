#load packages
library(readxl)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(whippr)
library(lubridate)
library(segmented)
library(car)
library(tibble)

##monoexponential fit
# data clean --------------------------------------------------------------
CPET_CWR <- read_data(".xlsx", metabolic_cart = "nspire", time_column = "CPET Results")

CPET_CWR <- CPET_CWR |> 
  rename_with(~ as.character(seq_along(.)), everything()) |> 
  dplyr::select("1","4", "12") |> 
  rename(time = '1',
         vo2 = '4',
         po  = '12') |> 
  mutate(vo2 = as.numeric(as.character(vo2)))

first_80_index <- which(CPET_CWR$po == 80)[1]
start_index <- max(1, first_80_index - 1)  
end_index <- which(CPET_CWR$po == 0 & seq_along(CPET_CWR$po) > which(CPET_CWR$po == 80) [1])[1] - 1
CPET_CWR <- CPET_CWR |> slice(start_index:end_index) 

CPET_CWR$time_corrected <- seq(from = 0, by = 5, length.out = nrow(CPET_CWR))


# plot data ---------------------------------------------------------------
ggplot(CPET_CWR, aes(x = time_corrected, y = vo2 )) +
  geom_point(shape = 21, size = 3, fill = "white") + 
  labs(x = "time (s)", y = "vo2 (L/min)", title = "vo2/time") +
  theme_whippr() +
  scale_y_continuous(limits = c(0,4)) +
  scale_y_continuous(breaks = seq(0, 3, by = 0.5) )

# convert time ------------------------------------------------------------
CPET_CWR <- CPET_CWR |> 
  mutate(time = period_to_seconds(hms(time)))

CPET_CWR <- CPET_CWR |> 
  rename("CPET Results" = time_corrected) |> 
  dplyr::select(4,2,3)

# Regression analysis ------------------------------------------------------
Results <- vo2_kinetics(
  .data = CPET_CWR,
  intensity_domain = "moderate",
  vo2_column = vo2,
  protocol_n_transitions = 1,
  protocol_baseline_length = 5,
  protocol_transition_length = 355,
  cleaning_level = 0.90,
  cleaning_baseline_fit = "linear",
  fit_level = 0.90,
  fit_bin_average = 5,
  fit_phase_1_length = 20,
  fit_baseline_length = 5,
  fit_transition_length = 355,
  verbose = TRUE
)


Results$plot_model[[1]]
Results$model_summary

# remove outliers
Results <- Results[[1]][[1]] |> 
  dplyr::filter(outlier != "yes")

last_120_seconds <- mean(tail(Results$VO2,12, na.rm = TRUE))
 


##ramp fit
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
CPET_RAMP$time <- CPET_RAMP$time - intercept 

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
    method_incremental = "anomaly", 
    verbose = TRUE) |>   
  plot_outliers() 

# remove outliers
data_ramp_without_outliers <- data_ramp_with_outliers[["data"]]|> 
  dplyr::filter(outlier != "yes") 


data_ramp_without_outliers <- data_ramp_without_outliers|>
  rename("time" = `CPET Results`) 
data_ramp_without_outliers$po <- as.numeric(data_ramp_without_outliers$po)  

#fit2
df_filtered <- data_ramp_without_outliers |> 
  filter(time >= 70, vo2 <= 2.6)


fit2 <- lm(vo2 ~ po, data = df_filtered)

Intercept_treshold <- (fit2[["coefficients"]][["(Intercept)"]] - 2.6)/(0 - fit2[["coefficients"]][["po"]])

filtered_data_fit2 <- data_ramp_without_outliers |> 
  filter(time >= 70, po <= Intercept_treshold)

fit2 <- lm(vo2 ~ po, data = filtered_data_fit2)




#time delay
time_delay_po <- (fit2[["coefficients"]][["(Intercept)"]] - last_120_seconds)/(0 - fit2[["coefficients"]][["po"]])

# Plot 
ggplot(data_ramp_without_outliers, aes(x = po, y = vo2 )) +
  geom_point(shape = 21, size = 3, fill = "white") + 
  labs(x = "power output (W)", y = "VO2 (L/min)", title = "vo2/po") +
  scale_x_continuous(breaks = seq(10, 400, by = 20)) +
  theme_whippr() +
  
  
  geom_abline(intercept = last_120_seconds , slope = 0, color = "red", linetype = "solid", size = 1) +
  
  
  geom_abline(intercept = coef(fit2)[["(Intercept)"]], slope = fit2[["coefficients"]][["po"]], 
              color = "red", linetype = "solid", size = 1) +
  geom_vline(xintercept = Intercept_treshold, linetype = "dashed", color = "red")


time_delay_x <- (time_delay_po - 80) * 6

slope_fit2 <- coefficients(fit2)[["po"]]
se_slope2 <- summary(fit2)$coefficients["po", "Std. Error"]
t_value_slope2 <- summary(fit2)$coefficients["po", "t value"]
p_value_slope2 <- summary(fit2)$coefficients["po", "Pr(>|t|)"]
conf_intervals_fit2 <- confint(fit2, level = 0.95)
ci_low_slope2 <- conf_intervals_fit2[2,1]
ci_high_slope2 <- conf_intervals_fit2[2,2]



table <- tibble(
  term = c("VO2_CWR", "slope 2", "MRT"),
  estimate = c(round(last_120_seconds,5), round(slope_fit2,5), round(time_delay_x,5)),
  std.error = c("", formatC(se_slope2, format = "E", digits = 5), ""),
  t.value = c("", round(t_value_slope2,5), ""),
  p.value = c("", formatC(p_value_slope2, format = "E", digits = 5), ""),
  conf.low = c("", round(ci_low_slope2,5),""),
  conf.high = c("", round(ci_high_slope2,5),""),
)  



