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
    cleaning_level = 0.95,
    method_incremental = "anomaly", ## changed to anomaly detection
    verbose = TRUE) |>   
  plot_outliers() 

# remove outliers
data_ramp_without_outliers <- data_ramp_with_outliers[["data"]]|> 
  dplyr::filter(outlier != "yes") 


data_ramp_without_outliers <- data_ramp_without_outliers|>
  rename("time" = `CPET Results`) 
data_ramp_without_outliers$po <- as.numeric(data_ramp_without_outliers$po)  


#fit2
df_filtered <- data_ramp_without_outliers[
  which(data_ramp_without_outliers$vo2 >= 2.6 & data_ramp_without_outliers$vo2 <= 3.1),
]


fit2 <- lm(vo2 ~ po, data = df_filtered)
Intercept_GET <- (fit2[["coefficients"]][["(Intercept)"]] - 2.6)/(0 - fit2[["coefficients"]][["po"]])
Intercept_RCP <- (fit2[["coefficients"]][["(Intercept)"]] - 3.1)/(0 - fit2[["coefficients"]][["po"]])

filtered_data_fit2 <- data_ramp_without_outliers |> 
  filter(po >= Intercept_GET, po <= Intercept_RCP)

fit3 <- lm(vo2 ~ po, data = filtered_data_fit2)



#plot
ggplot(data_ramp_without_outliers, aes(x = po, y = vo2)) +
  geom_point(shape = 21, size = 3, fill = "white") +  
  geom_abline(intercept = coef(fit3)[["(Intercept)"]], slope = fit3[["coefficients"]][["po"]], 
              color = "red", linetype = "solid", size = 1) +  
  labs(title = "VO2/po",
       x = "Power output (W)",
       y = "VO2 (L/min)") +
  scale_x_continuous(breaks = seq(10, 400, by = 20))  +
  geom_vline(xintercept = Intercept_GET, linetype = "dashed", color = "red") +  # Eerste verticale lijn
  geom_vline(xintercept = Intercept_RCP, linetype = "dashed", color = "red") +
  theme_whippr()
  


slope_fit3 <- coefficients(fit3)[["po"]]
se_slope3 <- summary(fit3)$coefficients["po", "Std. Error"]
t_value_slope3 <- summary(fit3)$coefficients["po", "t value"]
p_value_slope3 <- summary(fit3)$coefficients["po", "Pr(>|t|)"]
conf_intervals_fit3 <- confint(fit3, level = 0.95)
ci_low_slope3 <- conf_intervals_fit3[2,1]
ci_high_slope3 <- conf_intervals_fit3[2,2]

table <- tibble(
  term = "slope GET-RCP",
  estimate = round(slope_fit3,5),
  std.error = formatC(se_slope3, format = "E", digits = 5),
  t.value = round(t_value_slope3,5),
  p.value = formatC(p_value_slope3, format = "E", digits = 5),
  conf.low = round(ci_low_slope3,5),
  conf.high = round(ci_high_slope3,5),
)  


