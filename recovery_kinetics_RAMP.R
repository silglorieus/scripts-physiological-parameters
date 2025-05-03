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


# data clean --------------------------------------------------------------
CPET_RAMP_1 <- read_data(".xlsx", metabolic_cart = "nspire", time_column = "CPET Results")
CPET_RAMP <- CPET_RAMP_1 |> rename_with(~ as.character(1:16), everything()) |> 
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
ramp_normalized <- CPET_RAMP_1 |> 
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


Results <- data_ramp_without_outliers|> rename("time" = `CPET Results`)
Results$po <- as.numeric(Results$po)


last_60_seconds_start <- mean(tail(Results$vo2, 12), na.rm = TRUE)


#recovery kinetics
Results <- CPET_RAMP_1 |> rename_with(~ as.character(1:16), everything()) |> 
  dplyr::select(1,4,12) |> 
  rename(time = "1",
         vo2 = "4",
         po  = "12") |> 
  mutate(vo2 = as.numeric(as.character(vo2))) |> 
  slice(which(po == 0) [1]:n())
start_index <- which(Results$po == 50) [1]
end_index <- nrow(Results)
Results <- Results |> slice(start_index:end_index) 


Results <- Results |> 
  mutate(time = period_to_seconds(hms(time)))
intercept <- tail(Results$time[Results$po == 50],1)
Results$time <- Results$time - intercept

Results$po <- as.numeric(Results$po)

last_po_max_index <- tail(which(Results$po == max(Results$po, na.rm = TRUE)), 1) 
end_index <- nrow(Results)


Results <- Results |> 
  slice(last_po_max_index:end_index)


last_60_seconds_end <- mean(tail(Results$vo2,12, na.rm = TRUE))



a_fixed <- last_60_seconds_start - last_60_seconds_end
c_fixed <- last_60_seconds_end

Results$t0 <- Results$time - min(Results$time)




model <- nls(vo2 ~ a_fixed * exp(-b * t0) + c,
             start = list(b = 0.01, c = c_fixed),
             data = Results,
              )

summary(model)

plot(Results$time, Results$vo2, pch = 1, col = "black", bg = "lightyellow", xlab = "t", ylab = "VO2") +
  lines(Results$time, predict(model), col = "red", lwd = 2) +
  theme_whippr()



model_summary <- summary(model)
ci <- confint(model)
b_fixed     <- model_summary$parameters["b", "Estimate"]
std.error   <- model_summary$parameters["b", "Std. Error"]
t.value     <- model_summary$parameters["b", "t value"]
p.value     <- model_summary$parameters["b", "Pr(>|t|)"]
c_fixed     <- model_summary$parameters["c", "Estimate"]
c_std.error   <- model_summary$parameters["c", "Std. Error"]
c_t.value     <- model_summary$parameters["c", "t value"]
c_p.value     <- model_summary$parameters["c", "Pr(>|t|)"]
a_fixed     <- last_60_seconds_start - c_fixed
b_low <- ci["b", 1]
b_high <- ci["b", 2]
c_low <- ci["c", 1]
c_high <- ci["c", 2]
tau <- 1 / b_fixed
T_half <- 0.693 / b_fixed

table <-  tibble(
  Term = c("a", "b", "c", "tau", "T1/2"),
  estimate = c(round(a_fixed,5), round(b_fixed,5), round(c_fixed,5), round(tau,5), round(T_half ,5)),
  std.error = c("", formatC(std.error, format = "E", digits = 5), formatC(c_std.error, format = "E", digits = 5),"",""),
  t.value = c("", round(t.value,5), round(c_t.value,5),"",""),
  p.value = c("", formatC(p.value, format = "E", digits = 5), formatC(c_p.value, format = "E", digits = 5),"",""),
  conf.low = c("", round(b_low,5),round(c_low,5),"",""),
  conf.high = c("", round(b_high,5), round(c_high,5),"",""),
)

plot(CPET_RAMP$`CPET Results`, CPET_RAMP$vo2, pch = 1, col = "black", bg = "lightyellow", xlab = "t", ylab = "VO2")


# Plot 
ggplot(Results, aes(x = time, y = vo2 )) +
  geom_point(shape = 21, size = 3, fill = "white") + 
  labs(x = "time (s)", y = "vo2 (L/min)", title = "vo2/time") +
  theme_whippr() +
  geom_line(aes(y = predict(model)), color = "red")




