### July 2022 Site Characteristics metadata compilation
### Created by Danielle Barnas
### Created on 10/26/2022

#### LOAD LIBRARIES ####
library(tidyverse)
library(here)


#### BRING IN DATA ####
turb <- read_csv(here("Data", "Biogeochem", "July2022", "Turb_NC.csv"))
dist <- read_csv(here("Data", "Plate_Distance_to_Seep.csv"))
depth <- read_csv(here("Data", "Adj_Sandwich_Depth.csv"))
sub <- read_csv(here("Data", "Surveys", "Substrate_2022.csv"))
rug <- read_csv(here("Data","Surveys","Survey_Metadata.csv"))


#### PROCESS DATA ####
sub <- sub %>%
  group_by(Location, CowTagID) %>%
  mutate(TotalSub = sum(LiveCoral, DeadCoral, Rubble, Sand)) %>%
  mutate_at(.vars = vars(LiveCoral:Sand), .funs = ~(. / TotalSub *100)) %>% # ~ and . allow the function to run across all variables
  select(-c(Date, TotalSub))

rug <- rug %>%
  select(Location, CowTagID, Chain1, Chain2, Chain3) %>%
  group_by(Location, CowTagID) %>%
  summarise(meanChainLength = mean(c(Chain1, Chain2, Chain3))) %>%
  mutate(meanRugosity = meanChainLength / 2.03) %>%  # chain length of 2.03m
  select(-meanChainLength) %>%
  mutate(Location = if_else(Location == "Varari_Maya", "Varari", Location))

#### BRING TOGETHER ####
metadata <- turb %>%
  full_join(sub) %>%
  full_join(dist) %>%
  full_join(depth) %>%
  full_join(rug) %>%
  select(-c(Dist_CT_cm, Depth_cm, Tidal_diff))

write_csv(metadata, here("Data","Full_Metadata.csv"))

#### CREATE DATA DICTIONARY ####

Attribute <- c("Location",
               "CowTagID",
               "del15N",
               "C_N",
               "N_percent",
               "LiveCoral",
               "DeadCoral",
               "Rubble",
               "Sand",
               "lat",
               "lon",
               "dist_to_seep_m",
               "adj_CT_depth_cm",
               "meanRugosity"
               )
Units <- c(NA,
           NA,
           "unitless",
           "unitless",
           "%",
           "%",
           "%",
           "%",
           "%",
           NA,
           NA,
           "meters (m)",
           "centimeters (cm)",
           "meters (m)"
           )
Description <- c("One of two coral reef survey sites (Cabral or Varari)",
                 "Individual survey location identifier (1-20 and a Seep location)",
                 "Ratio of Nitrogen-15 / Nitrogen-14 from tissue samples of Turbinaria ornata at each survey location collected in July 2022",
                 "Ratio of Carbon / Nitrogen from tissue samples of Turbinaria ornata at each survey location collected in July 2022",
                 "Percent of Nitrogen in tissue samples of Turbinaria ornata at each survey location collected in July 2022",
                 "Percent cover of live coral substrate",
                 "Percent cover of dead coral substrate",
                 "Percent cover of rubble substrate",
                 "Percent cover of sand substrate",
                 "Latatitude of survey locations recorded on a Garmin GPS",
                 "Longitude of survey locations recorded on a Garmin GPS",
                 "Linear distance from each survey location to the location labeled as SEEP",
                 "Ajusted distance (or depth) measured from water surface to HOBO Conductiity-Temperature logger sensor at each survey location using a transect tape. Corrected for tidal changes",
                 "Average rugosity measured using a 2.03m link chain across 3 randomly chosen lines within the survey box of each survey location. Values range from 0-1, with 0 indicating higher substrate complexity and 1 indicating a flat surface"
                 )

metadata_dict <- as_tibble(cbind(Attribute, Units, Description))

write_csv(metadata_dict, here("Data", "Site_Metadata_Data_Dictionary.csv"))

