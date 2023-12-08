#Simulation of Brain Growth Charts for First Year of Life, Stratified by Prematurity Status
#Assuming WM growth follows an approx linear trajectory
#Code generated with assistance from ChatGPT


# Install and load required packages
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}
library(tidyverse)

# Set seed for reproducibility
set.seed(123)

# Function to generate base dataset for wm_volume
generate_base_data <- function() {
  data.frame(
    age_weeks = seq(1, 53, by = 1),
    wm_volume = seq(27, 40, length.out = 53) + rnorm(53, mean = 1, sd = 0)
  )
}

# Function to generate datasets with different age offsets and adjusted growth rates
generate_offset_datasets <- function(base_data, offsets) {
  datasets <- list()
  for (offset in offsets) {
    datasets[[as.character(offset)]] <- base_data %>%
      mutate(dataset = as.character(offset),
             wm_volume = offset * wm_volume - (offset*10-5))  # Adjust growth rate and starting value
  }
  datasets
}

# # Function to generate datasets with different age offsets
# generate_offset_datasets <- function(base_data, offsets) {
#   datasets <- list()
#   for (offset in offsets) {
#     datasets[[as.character(offset)]] <- base_data %>%
#       mutate(dataset = as.character(offset),
#              wm_volume = wm_volume - offset)
#   }
#   datasets
# }

# Generate base dataset for wm_volume
base_data <- generate_base_data()

# Generate datasets with different age offsets
offsets <- c(0.71, 0.76, 0.82, 0.9, 1)
offset_datasets <- generate_offset_datasets(base_data, offsets)

# Combine datasets into one
all_data <- bind_rows(offset_datasets)

# Calculate percent max value of wm_volume across all datasets
all_data <- all_data %>%
  group_by(dataset)%>%
  mutate(wm_percent_max = wm_volume / max(wm_volume) * 100)%>% ungroup()

# # Calculate percent max value of wm_volume within each dataset
# all_data <- all_data %>%
#   group_by(dataset) %>%
#   mutate(wm_percent_max = wm_volume / max(wm_volume) * 100) %>%
#   ungroup()

# Plot wm_volume against percent max value with color mapping
ggplot(all_data, aes(x = age_weeks, y = wm_volume, color = factor(dataset, levels = c("0.71", "0.76", "0.82", "0.9", "1")))) +
  geom_line(size = 1) +
  scale_color_manual(values = c("0.71" = "#ADEAEA", "0.76" = "#64C8D8", "0.82" = "#29A7CA", "0.9" = "#0073D9", "1" = "#004080"),
                     labels=c('Extremely Preterm', 'Very Preterm', 'Moderately Preterm','Late Preterm', 'At Term')) +
  labs(title = "Growth Trajectory of WM Volume for Level of Prematurity",
       x = "Age (weeks)", y = "Volume (% of maximum)") +
  theme_minimal() +
  guides(color = guide_legend(title = "Prematurity",
                              reverse = TRUE,
                              labels = c("0.71" = "Extremely Preterm",
                                         "0.76" = "Very Preterm",
                                         "0.82" = "Moderately Preterm",
                                         "0.9" = "Late Preterm",
                                         "1" = "At Term")))
