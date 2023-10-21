library(ismev)
#install.packages("VGAM")
library(VGAM)
data("venice90")
venice <- venice90
library(tidyverse)
max_values <- data.frame(year = integer(0), max_value = double(0))
unique_years <- data.frame(unique(venice$year))


  

# Create an empty data frame to store the results
max_sea_levels <- data.frame(year = numeric(0), max_sea_level = numeric(0))

# Get unique year values from your "venice" data frame
unique_years <- unique(venice$year)

# Loop through each unique year
for (year in unique_years) {
  # Subset the data for the current year
  subset_data <- venice[venice$year == year, ]
  
  # Calculate the maximum sea level for the current year
  max_level <- max(subset_data$sea_level)
  
  # Create a data frame with the year and maximum sea level
  result_row <- data.frame(year = year, max_sea_level = max_level)
  
  # Append the result to the max_sea_levels data frame
  max_sea_levels <- rbind(max_sea_levels, result_row)
}

# Print the resulting data frame
print(max_sea_levels)

