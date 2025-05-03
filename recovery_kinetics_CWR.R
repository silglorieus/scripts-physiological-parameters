#load packages
library(readxl)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(whippr)
library(lubridate)
library(segmented)
library(purrr)
library(lubridate)

# data clean --------------------------------------------------------------
# Function to read all files
read_all_excels <- function(directory) {
  
  files <- list.files(directory, pattern = "\\.xlsx$", full.names = TRUE)
  
  
  data_list <- lapply(files, function(file) {
    read_data(file, metabolic_cart = "nspire", time_column = "CPET Results")
  })
  
  names(data_list) <- basename(files)  
  
  return(data_list)
}

# perform function for the map "JDK ..."
CPET_CWR_list <- read_all_excels("JDK 1")

#Function to analyse one file
process_CPET_data <- function(data) {
  data <- data |> 
    rename_with(~ as.character(seq_along(.)), everything()) |> 
    dplyr::select(1, 4, 12) |> 
    rename(time = "1", vo2 = "4", po = "12") |> 
    dplyr::mutate(vo2 = as.numeric(as.character(vo2)))
  
  

  first_80_index <- which(data$po == 80)[1]
  start_index <- max(1, first_80_index - 1)  
  
  first_50_index <- which(data$po == 50)[1]
  end_index <- tail(which(data$po == 0 & seq_along(data$po) < first_50_index), 1)
    if (!is.na(start_index) && !is.na(end_index) && start_index < end_index) {
    data <- data |> slice(start_index:end_index)
  }
  data$time_corrected <- seq(from = 0, by = 5, length.out = nrow(data))
  
  return(data)
}

# apply function on all files in the list
processed_CPET_CWR_list <- lapply(CPET_CWR_list, process_CPET_data)

#merge datalist
merged <- reduce(processed_CPET_CWR_list,left_join, by = "time_corrected")
CPET_CWR <- as_tibble(merged) |>  
  dplyr::select(1,4,2,6,9,12,13)

#average
CPET_CWR <- CPET_CWR |> 
  mutate(mean = rowMeans(CPET_CWR[,3:6], na.rm = TRUE))

# plot data ---------------------------------------------------------------
ggplot(CPET_CWR, aes(x = time_corrected, y = mean )) +
  geom_point(shape = 21, size = 3, fill = "white") + 
  labs(x = "time (s)", y = "vo2 (L/min)", title = "vo2/time") +
  theme_whippr() +
  scale_y_continuous(limits = c(0, 4)) +
  scale_y_continuous(breaks = seq(0, 3, by = 0.5) )

CPET_CWR <- CPET_CWR |> 
  rename("CPET Results" = time_corrected) |> 
  rename("vo2" = mean) |> 
  rename("po" = po.y.y)

CPET_CWR <- CPET_CWR[, -c(1, 3, 4, 5, 6)]




# Regression analysis ------------------------------------------------------
Results <- vo2_kinetics(
  .data = CPET_CWR,
  intensity_domain = "moderate",
  vo2_column = vo2,
  protocol_n_transitions = 1,
  protocol_baseline_length = 360,
  protocol_transition_length = 360,
  cleaning_level = 0.90,
  cleaning_baseline_fit = c("exponential"),
  fit_level = 0.90,
  fit_bin_average = 5,
  fit_phase_1_length = 20,
  fit_baseline_length = 5,
  fit_transition_length = 360,
  verbose = TRUE
)


Results$plot_model[[1]]
Results$model_summary

# Remove outliers
Results <- Results[[1]][[1]] |> 
  dplyr::filter(outlier != "yes")



last_60_seconds_start <- mean(tail(Results$VO2[Results$po == 80], 12), na.rm = TRUE)

last_60_seconds_end <- mean(tail(Results$VO2,12, na.rm = TRUE))


last_po_80_index <- tail(which(Results$po == 80),1 ) 
Results <- Results |> 
  slice(last_po_80_index:nrow(Results))

a_fixed <- last_60_seconds_start - last_60_seconds_end
c_fixed <- last_60_seconds_end

Results$t0 <- Results$t - min(Results$t)


model <- nls(VO2 ~ a_fixed * exp(-b * t0) + c,
             start = list(b = 0.01, c = c_fixed),
             data = Results)

summary(model)

plot(Results$t, Results$VO2, pch = 1, col = "black", bg = "lightyellow", xlab = "t", ylab = "VO2") +
  lines(Results$t, predict(model), col = "red", lwd = 2) +
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



# Plot 
ggplot(Results, aes(x = t, y = VO2 )) +
  geom_point(shape = 21, size = 3, fill = "white") + 
  labs(x = "time (s)", y = "vo2 (L/min)", title = "vo2/time") +
  theme_whippr() +
  geom_line(aes(y = predict(model)), color = "red")
  