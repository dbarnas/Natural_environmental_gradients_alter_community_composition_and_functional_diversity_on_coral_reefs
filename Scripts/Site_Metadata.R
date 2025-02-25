### July 2022 Site Characteristics metadata compilation
### Created by Danielle Barnas
### Created on 10/26/2022

#### LOAD LIBRARIES ####
library(tidyverse)
library(here)


#### BRING IN DATA ####
# turb <- read_csv(here("Data", "Biogeochem", "Turb_NC.csv"))
# dist <- read_csv(here("Data", "Plate_Distance_to_Seep.csv"))
# depth <- read_csv(here("Data", "Adj_Sandwich_Depth.csv"))
sub <- read_csv(here("Data", "Surveys", "Substrate_2022.csv"))
rug <- read_csv(here("Data","Surveys","Survey_Metadata.csv"))
# alpha <- read_csv(here("Data", "CowTag_to_AlphaTag.csv"))


#### PROCESS DATA ####
# substrate
sub <- sub %>%
  group_by(Location, CowTagID) %>%
  mutate(TotalSub = sum(LiveCoral, DeadCoral, Rubble, Sand)) %>%
  mutate_at(.vars = vars(LiveCoral:Sand), .funs = ~(. / TotalSub *100)) %>% # ~ and . allow the function to run across all variables
  select(-c(Date, TotalSub))

# rugosity; complexity
rug <- rug %>%
  select(Location, CowTagID, AlphaTag, Chain1, Chain2, Chain3, lat, lon) %>%
  group_by(Location, CowTagID, AlphaTag, lat, lon) %>%
  summarise(meanChainLength = mean(c(Chain1, Chain2, Chain3))) %>%
  mutate(meanRugosity = meanChainLength / 2.03) %>%  # chain length of 2.03m
  select(-meanChainLength) %>%
  mutate(meanRugosity = if_else(CowTagID == "VSEEP", 0.97, meanRugosity)) %>% # transcribed mean rugosity at Seep
  mutate(complexity = 1-meanRugosity)

#### BRING TOGETHER ####
metadata <- sub %>%
  full_join(rug) %>%
  filter(Location == "Varari") %>%
  drop_na(AlphaTag) %>%  # only study survey sites
  relocate(AlphaTag, .after = CowTagID)

write_csv(metadata, here("Data","Full_Metadata.csv"))

#### CREATE DATA DICTIONARY ####
#
# Attribute <- c("Location",
#                "CowTagID",
#                "AlphaTag",
#                "LiveCoral",
#                "DeadCoral",
#                "Rubble",
#                "Sand",
#                "lat",
#                "lon",
#                "meanRugosity",
#                "complexity"
#                )
# Units <- c(NA,
#            NA,
#            NA,
#            "%",
#            "%",
#            "%",
#            "%",
#            NA,
#            NA,
#            "unitless",
#            "unitless"
#            )
# Description <- c("One of two coral reef survey sites (Cabral or Varari)",
#                  "Individual survey location identifier (1-20 and a Seep location)",
#                  "Individual survey location identifier, alphabetized by linear distance from seepage point (Seep location, A, and B-T)",
#                  "Percent cover of live coral substrate",
#                  "Percent cover of dead coral substrate",
#                  "Percent cover of rubble substrate",
#                  "Percent cover of sand substrate",
#                  "Latitude of survey locations recorded on a Garmin GPS",
#                  "Longitude of survey locations recorded on a Garmin GPS",
#                  "Average rugosity measured using a 2.03m link chain across 3 randomly chosen lines within the survey box of each survey location. Values range from 0-1, with 0 indicating higher substrate complexity and 1 indicating a flat surface",
#                  "Average rugosity subtracted from 1 as an indication of structural complexity of the reef surface"
#                  )
#
# metadata_dict <- as_tibble(cbind(Attribute, Units, Description))
#
# write_csv(metadata_dict, here("Data", "Site_Metadata_Data_Dictionary.csv"))

