#############################################
## The following loads the needed packages ##
#############################################

# load the required packages
packages <- c(
  "dplyr", # practical 1 & 2 & 3
  "QRM",  # practical 1
  "tseries",  # practical 1
  "nortest",  # practical 1
  "forecast",  # practical 1
  "plotly",  # practical 1 & 2
  "MASS",# practical 1 
  "TTR", # practical 1
  "car", # practical 1
  "fGarch", # practical 1
  "ismev", # practical 2 
  "VGAM", # practical 2 
  "tidyverse", # practical 2 
  "plotrix", # practical 2 
  "scales", # practical 2 
  "extRemes",# practical 2 
  "knitr", # practical 2 
  "xts",# practical 2 
  "kableExtra", # practical 2 
  "dygraphs",# practical 2 
  "gridExtra",# practical 2 
  "readr",# practical 2 
  "lubridate",# practical 3
  "evd",# practical 2 
  "xts",# practical 2 
  "data.table",# practical 3
  "copula",# practical 3
  "VineCopula",
  "psych")# practical 3

# Check if each package is available, and install it if not
for (package in packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  }
}

# Now load each package
lapply(packages, function(pkg) {
  library(pkg, character.only = TRUE)
})
######################################################
## The following sets a few option for nice reports ##
######################################################

# general options
options(
  digits = 9,
  str = strOptions(strict.width = "cut"),
  width = 69,
  tibble.width = 69,
  cli.unicode = FALSE
)

# ggplot options
theme_set(theme_light())


t <- list(
  family = "Times New Roman",
  size = 14)


