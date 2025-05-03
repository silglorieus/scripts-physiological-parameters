# load packages -----------------------------------------------------------
install.packages("whippr")
install.packages("readxl")
install.packages("ggplot2")
install.packages("tidyverse")
install.packages("dplyr")
install.packages("lubridate")
install.packages("segmented")
install.packages("tibble")
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
CPET_RAMP <- read_data("JDK 1/CPET_JDK_01_Colosio_Alessandro_10W_5sec.xlsx", metabolic_cart = "nspire", time_column = "CPET Results")
CPET_RAMP <- CPET_RAMP |> rename_with(~ as.character(1:16), everything()) |> 
    dplyr::select(1,4,12) |> 
    rename(time = "1",
           vo2 = "4",
           po  = "12") |> 
    mutate(vo2 = as.numeric(as.character(vo2))) |> 
    slice(which(po == 0) [1]:n())
start_index <- which(CPET_RAMP$po == 50) [1]
end_index <- nrow(CPET_RAMP)

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
      baseline_length = 0, 
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


data_ramp_without_outliers <- data_ramp_without_outliers|> rename("time" = `CPET Results`)




# regression analysis -----------------------------------------------------
intercept_index <- tail(which(CPET_RAMP$po == 50), 1)


# fit1
fit1 <- lm(vo2 ~ 1, data = data_ramp_without_outliers[1:intercept_index,])
conf_intervals_fit1 <- confint(fit1, level = 0.95)
error_margin_slope1 <- conf_intervals_fit1[1, 2] - conf_intervals_fit1[1, 1] 

#fit2
filtered_data_fit2 <- data_ramp_without_outliers |> 
  filter(time >= 70, vo2 <= 2.6)

fit2 <- lm(vo2 ~ time, data = filtered_data_fit2)

Intercept_GET <- (fit2[["coefficients"]][["(Intercept)"]] - 2.6)/(0 - fit2[["coefficients"]][["time"]])

filtered_data_fit2 <- data_ramp_without_outliers |> 
  filter(time >= 70, time <= Intercept_GET)

fit2 <- lm(vo2 ~ time, data = filtered_data_fit2)


conf_intervals_fit2 <- confint(fit2, level = 0.95)
error_margin_slope2 <- conf_intervals_fit2[2, 2] - conf_intervals_fit2[2, 1] 


# Time delay ------------------------------------------------------
Time_delay_x <- ( fit2[["coefficients"]][["(Intercept)"]] - fit1[["coefficients"]][["(Intercept)"]])/ (0 - fit2[["coefficients"]][["time"]]) 
Time_delay_y <- 0 * Time_delay_x + fit1[["coefficients"]][["(Intercept)"]]

# Plot 
ggplot(data_ramp_without_outliers, aes(x = time, y = vo2 )) +
  geom_point(shape = 21, size = 3, fill = "white") + 
  labs(x = "time (s)", y = "vo2 (L/min)", title = "vo2/time") +
  theme_whippr() +
  

  geom_abline(intercept = fit1[["coefficients"]][["(Intercept)"]], slope = 0, color = "red", linetype = "solid", size = 1) +
  
  
  geom_abline(intercept = coef(fit2)[["(Intercept)"]], slope = coef(fit2)[["time"]], 
              color = "red", linetype = "solid", size = 1)




# Extract regression parameters
intercept_fit1 <- fit1[["coefficients"]][["(Intercept)"]]  # Intercept van fit1 (vast)
slope_fit1 <- 0              # Slope van fit1 (altijd 0)
intercept_fit2 <- coef(fit2)[["(Intercept)"]]
slope_fit2 <- coef(fit2)[["time"]]


# Make textstring
annotation_text <- paste0(
  "i1: ", round(intercept_fit1, 5), " ± ", round(error_margin_slope1,5), "\n",
  "s1: ", round(slope_fit1, 2), "\n",
  "TD: ", round(Time_delay_x, 5), "\n",
  "s2: ", round(slope_fit2, 5), " ± ", round(error_margin_slope2, 5)
  )


# Plot 
ggplot(data_ramp_without_outliers, aes(x = time, y = vo2 )) +
  geom_point(shape = 21, size = 3, fill = "white") + 
  labs(x = "time (s)", y = "vo2 (L/min)", title = "vo2/time") +
  theme_whippr() +
  
geom_abline(intercept = fit1[["coefficients"]][["(Intercept)"]], slope = 0, color = "red", linetype = "solid", size = 1) +
  

geom_abline(intercept = coef(fit2)[["(Intercept)"]], slope = coef(fit2)[["time"]], 
              color = "red", linetype = "solid", size = 1) +
  

geom_vline(xintercept = Time_delay_x, color = "black", linetype = "dashed", size = 1) + 

geom_vline(xintercept = Intercept_GET, linetype = "dashed", color = "red")


# add text
annotate("text", x = 200, y = max(data_ramp_without_outliers$vo2), 
         label = annotation_text, hjust = 0, vjust = 1, 
         color = "black", size = 3, fontface = "bold", 
         family = "sans", lineheight = 1.2, parse = FALSE) 





#parameters
se_intercept1 <- summary(fit1)$coefficients["(Intercept)", "Std. Error"]
se_intercept2 <- summary(fit2)$coefficients["(Intercept)", "Std. Error"]
se_slope2 <- summary(fit2)$coefficients["time", "Std. Error"]
t_value_intercept1 <- summary(fit1)$coefficients[, "t value"]
p_value_intercept1 <- summary(fit1)$coefficients[, "Pr(>|t|)"]
t_value_slope2 <- summary(fit2)$coefficients["time", "t value"]
p_value_slope2 <- summary(fit2)$coefficients["time", "Pr(>|t|)"]
ci_low_intercept1 <- conf_intervals_fit1[,1]
ci_high_intercept1 <- conf_intervals_fit1[,2]
ci_low_slope2 <- conf_intervals_fit2[2,1]
ci_high_slope2 <- conf_intervals_fit2[2,2]
  
  

table <-  tibble(
  Term = c("intercept", "slope 2", "MRT"),
  estimate = c(round(intercept_fit1,5), round(slope_fit2,5), round(Time_delay_x,5)),
  std.error = c(round(se_intercept1,5), formatC(se_slope2, format = "E", digits = 5), ""),
  t.value = c(round(t_value_intercept1,5), round(t_value_slope2,5), ""),
  p.value = c(formatC(p_value_intercept1, format = "E", digits = 5), formatC(p_value_slope2, format = "E", digits = 5), ""),
  conf.low = c(round(ci_low_intercept1,5), round(ci_low_slope2,5),""),
  conf.high = c(round(ci_high_intercept1,5), round(ci_high_slope2,5),""),
  
)


