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
#turb1<-read_csv(here("Data","Biogeochem","August2021","Turb_NC.csv"))
MarchChemData <- read_csv(curl('https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data//March2022/CarbonateChemistry/pHProbe_Data_calculated_POcorrect.csv')) %>% mutate(Season = "Wet")
#turb2<-read_csv(here("Data","Biogeochem","March2022","Turb_NC.csv"))
#depth <- read_csv(here("Data","Adj_Sandwich_Depth.csv"))
turb <- read_csv(here("Data","Biogeochem", "July2022", "Turb_NC.csv"))


## Clean Data
## remove unnecessary/redundant turbinaria data to not include in cluster analysis
# turb1 <- turb1 %>%
#   select(CowTagID, del15N, C_N, N_percent) %>%
#   mutate(Season = "Dry")
# turb2 <- turb2 %>%
#   select(CowTagID, del15N, C_N, N_percent) %>%
#   mutate(Season = "Wet")
# turb <- full_join(turb1, turb2) %>% distinct()

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

## Remove outliers (samples run on different day from rest) from March
removeSite4 <- MarchChemData %>%
  filter(Season == "Wet",
         CowTagID %in% c("VSPRING",  # remove unnecessary sites from my analysis
                         #"V17",
                         "VRC",
                         "CRC",
                         "CPIT",
                         "CPIT_Bottom",
                         "CSPRING_ROAD",
                         "CSPRING_BEACH",
                         "CSPRING_BEACH2"))

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

MarchChem <- MarchChemData %>%
  rename(Temperature = TempInSitu) %>%
  select(-c(Date,
            SeepCode,
            SamplingTime,
            #Temperature,
            Notes)) %>%
  anti_join(removeSite4) %>%
  #left_join(turb, by = c('CowTagID','Season')) %>%
  left_join(gps, by = c('CowTagID'))


ReducedChemData <- full_join(AugChem, MarchChem) %>%
  relocate(Season, .after = Day_Night) %>%
  # remove clearly incorrect value of Ammonia_umolL at V17
  filter(Ammonia_umolL < 16)


##### Summarise: Max and Min of parameters across seasons ####
# range_data <- ReducedChemData %>%
#   group_by(Location, CowTagID, Season) %>%
#   # get parameter max and min across sampling periods by site and plate
#   #summarise_at(vars(Salinity:Tyrosine_Like), .funs = range, na.rm=T) %>% # returns two rows containing max and min value per CowTagID
#   summarise_at(vars(Salinity:Tyrosine_Like), .funs = diff, na.rm=T) %>%  # returns difference between rows per CowTagID (yields range at each location)
#   ungroup()

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


## Remove absent data (creating -Inf values)
anti <- full_data %>%
  filter(Season == "Wet",
         Parameters == "M_C" |
           Parameters == "HIX" |
           Parameters == "VisibleHumidic_Like" |
           Parameters == "Tryptophan_Like" |
           Parameters == "Tyrosine_Like")
full_data <- full_data %>%
  anti_join(anti)



##### Summarise: Mean of parameters - group by Seasons ####
mean_dataS <- ReducedChemData %>%
  group_by(Location, CowTagID, Season) %>%
  # Calculate mean values and pivot longer to join for CV calculation
  summarise_at(vars(Salinity:Tyrosine_Like), .funs = mean, na.rm = T) %>%
  ungroup() %>%
  pivot_longer(cols = c(Salinity:Tyrosine_Like), names_to = "Parameters", values_to = "MeanSeasonal")
mean_dataA <- ReducedChemData %>%
  group_by(Location, CowTagID) %>%
  # Calculate mean values and pivot longer to join for CV calculation
  summarise_at(vars(Salinity:Tyrosine_Like), .funs = mean, na.rm = T) %>%
  ungroup() %>%
  pivot_longer(cols = c(Salinity:Tyrosine_Like), names_to = "Parameters", values_to = "MeanAll")

## Join mean_data to larger data set
full_data <- full_data %>%
  left_join(mean_dataS) %>%
  left_join(mean_dataA)

##### Summarise: Coefficient of Variation of parameters - group by Seasons ####
cv_dataS <- ReducedChemData %>%
  group_by(Location, CowTagID, Season) %>%
  # Calculate mean values and pivot longer to join for CV calculation
  summarise_at(vars(Salinity:Tyrosine_Like), .funs = sd, na.rm = T) %>%
  ungroup() %>%
  pivot_longer(cols = c(Salinity:Tyrosine_Like), names_to = "Parameters", values_to = "sdS")
cv_dataA <- ReducedChemData %>%
  group_by(Location, CowTagID) %>%
  # Calculate mean values and pivot longer to join for CV calculation
  summarise_at(vars(Salinity:Tyrosine_Like), .funs = sd, na.rm = T) %>%
  ungroup() %>%
  pivot_longer(cols = c(Salinity:Tyrosine_Like), names_to = "Parameters", values_to = "sdA")

full_data <- full_data %>%
  left_join(cv_dataS) %>%
  left_join(cv_dataA) %>%
  mutate(CVSeasonal = sdS / MeanSeasonal,
         CVAll = sdA / MeanAll) %>%
  select(-c(sdS, sdA))
View(full_data)


## Write csv ####
write_csv(full_data, here("Data","Biogeochem","Nutrients_Processed_All.csv"))


