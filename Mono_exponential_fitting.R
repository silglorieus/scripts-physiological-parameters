# load packages -----------------------------------------------------------
install.packages("whippr")
install.packages("readxl")
install.packages("ggplot2")
install.packages("tidyverse")
install.packages("dplyr")
install.packages("lubridate")
install.packages("segmented")
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
# Function to read all excel files in a map
read_all_excels <- function(directory) {
  
  files <- list.files(directory, pattern = "\\.xlsx$", full.names = TRUE)
  

  data_list <- lapply(files, function(file) {
    read_data(file, metabolic_cart = "nspire", time_column = "CPET Results")
  })
  
  names(data_list) <- basename(files)  
  
  return(data_list)
}

# perform the function for "JDK ..."
CPET_CWR_list <- read_all_excels("JDK 1")

#Function to process one file
process_CPET_data <- function(data) {
  data <- data |> 
    rename_with(~ as.character(seq_along(.)), everything()) |> 
    dplyr::select(1, 4, 12) |> 
    rename(time = "1", vo2 = "4", po = "12") |> 
    dplyr::mutate(vo2 = as.numeric(as.character(vo2)))
  
  

  first_80_index <- which(data$po == 80)[1]
  start_index <- max(1, first_80_index - 1)  
  

  end_index <- which(data$po == 0 & seq_along(data$po) > first_80_index)[1] - 1
  if (!is.na(start_index) && !is.na(end_index) && start_index < end_index) {
    data <- data |> slice(start_index:end_index)
  }
  data$time_corrected <- seq(from = 0, by = 5, length.out = nrow(data))
 
  return(data)
}

# apply function on every file on the list
processed_CPET_CWR_list <- lapply(CPET_CWR_list, process_CPET_data)

#merge datalist
merged <- reduce(processed_CPET_CWR_list,left_join, by = "time_corrected")
CPET_CWR <- as_tibble(merged) |>  
  dplyr::select(4,2,6,9,,12,13)
  
#average
CPET_CWR <- CPET_CWR |> 
   mutate(mean = rowMeans(CPET_CWR[,2:5], na.rm = TRUE))

# plot data ---------------------------------------------------------------
ggplot(CPET_CWR, aes(x = time_corrected, y = mean )) +
  geom_point(shape = 21, size = 3, fill = "white") + 
  labs(x = "time (s)", y = "vo2 (L/min)", title = "vo2/time") +
  theme_whippr() +
  scale_y_continuous(limits = c(0,4)) +
  scale_y_continuous(breaks = seq(0, 3, by = 0.5) )

CPET_CWR <- CPET_CWR |> 
  rename("CPET Results" = time_corrected) |> 
    rename("vo2" = mean) 

  

# Regression analysis ------------------------------------------------------
Results <- vo2_kinetics(
  .data = CPET_CWR,
  intensity_domain = "moderate",
  vo2_column = vo2,
  protocol_n_transitions = 1,
  protocol_baseline_length = 5,
  protocol_transition_length = 350,
  cleaning_level = 0.95,
  cleaning_baseline_fit = "linear",
  fit_level = 0.95,
  fit_bin_average = 5,
  fit_phase_1_length = 20,
  fit_baseline_length = 5,
  fit_transition_length = 350,
  verbose = TRUE
)


Results$plot_model[[1]]
Results$model_summary

Results$plot_residuals[[1]]
Results$plot_outliers[[1]]

