#### Process August nutrient data and save as csv
#### Created by Danielle Barnas
#### Created on 8/25/2022

## Load Libraries
library(tidyverse)
library(lubridate)
library(here)
library(curl) # pull data from url


## Read in data
AugChemData<-read_csv(curl('https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/August2021/Allbiogeochemdata_QC.csv')) %>% mutate(Season = "Dry")
turb <- read_csv(here("Data","Biogeochem", "Turb_NC.csv"))


## Join august and march data


## There seems to be a contaminated nutrient sample for V2 Low tide on the 8/8/2021.
# Remove outliers / irrelevant data points
removeSite1 <- AugChemData %>%
  filter(Season == "Dry",
         CowTagID == "V2",
         Tide == " Low",
         Day_Night == "Day",
         Date == ymd("2021-08-08"))

removeSite2 <- AugChemData %>%
  filter(Season == "Dry",
         CowTagID == "C4",
         Tide =="Low",
         Day_Night == "Night",
         Date == ymd("2021-08-09"))

## Filter out redundant low tide day sample, first 'low tide' was super high
removeSite3 <- AugChemData %>%
  filter(Season == "Dry",
         Tide == "Low",
         Day_Night == "Day",
         Date == ymd("2021-08-06"))

removeSite5 <- AugChemData %>%
  filter(CowTagID %in% c("VSPRING", "Varari_Well", "CSPRING"))



## Pull out Location, lat, lon of CowTagIDs
gps <- AugChemData %>% select(Location, CowTagID, lat, lon) %>% distinct()

## Remove unnecessary/redundant data
AugChem <- AugChemData %>%
  select(-c(Date,
            Time,
            DateTime,
            Plate_Seep,
            Top_Plate_ID,
            Bottom_Plate_ID,
            Jamie_Plate_ID,
            #Temperature,
            # Craig recommends to remove "because we don't hypothesize them to be
            # orthogonal to any of the other fDOM we're using"
            MarineHumic_Like, Lignin_Like)) %>%
  anti_join(removeSite1) %>% # remove outlier/irrelevant data
  anti_join(removeSite2) %>%
  anti_join(removeSite3) %>%
  anti_join(removeSite5) %>%
  #left_join(turb, by = c('CowTagID','Season')) %>% # join with T. ornata nutrient loading data
  left_join(gps, by = c('Location','CowTagID','lat','lon'))

ReducedChemData <- AugChem %>%
  relocate(Season, .after = Day_Night) %>%
  # remove clearly incorrect value of Ammonia_umolL at V17 (false analytical value)
  filter(Ammonia_umolL < 16)


##### Summarise: Max and Min of parameters across seasons ####

max_data <- ReducedChemData %>%
  group_by(Location, CowTagID, Season) %>%
  summarise_at(vars(Salinity:Tyrosine_Like), .funs = max, na.rm = T) %>%  # select for max values
  pivot_longer(cols = Salinity:Tyrosine_Like, names_to = "Parameters", values_to = "Maximum")

min_data <- ReducedChemData %>%
  group_by(Location, CowTagID, Season) %>%
  summarise_at(vars(Salinity:Tyrosine_Like), .funs = min, na.rm = T) %>%  # select for max values
  pivot_longer(cols = Salinity:Tyrosine_Like, names_to = "Parameters", values_to = "Minimum")

# join max and min values and other data sets
full_data <- full_join(max_data, min_data) %>%
  mutate(Range = Maximum - Minimum)


## Join with GPS
full_data <- gps %>%
  right_join(full_data)



##### Summarise: Mean of parameters ####

mean_data <- ReducedChemData %>%
  group_by(Location, CowTagID, Season) %>%
  # Calculate mean values and pivot longer to join for CV calculation
  summarise_at(vars(Salinity:Tyrosine_Like), .funs = mean, na.rm = T) %>%
  ungroup() %>%
  pivot_longer(cols = c(Salinity:Tyrosine_Like), names_to = "Parameters", values_to = "Mean")


##### Summarise: Coefficient of Variation of parameters ####
cv_data <- ReducedChemData %>%
  group_by(Location, CowTagID, Season) %>%
  # Calculate mean values and pivot longer to join for CV calculation
  summarise_at(vars(Salinity:Tyrosine_Like), .funs = sd, na.rm = T) %>%
  ungroup() %>%
  pivot_longer(cols = c(Salinity:Tyrosine_Like), names_to = "Parameters", values_to = "sd")

full_data <- mean_data %>%
  full_join(cv_data) %>%
  mutate(CV = sd / Mean) %>%
  select(-c(sd))
View(full_data)


## Write csv ####
write_csv(full_data, here("Data","Biogeochem","Nutrients_Processed_All.csv"))


