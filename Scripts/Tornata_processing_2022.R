### July 2022 T. ornata nitrogen laoding cleaning script
### Created by Danielle Barnas
### Created on 10/13/2022

#### LOAD LIBRARIES ####
library(tidyverse)
library(here)


#### BRING IN DATA ####
raw <- read_csv(here("Data","Biogeochem","RawBarnasTornata.csv"))
rawMaya <- read_csv(here("Data","Biogeochem","RawZeffTornata.csv"))


turb <- raw %>%
  mutate(C_N = C_ug / N_ug,
         N_percent = (N_ug/1000) / Weight_mg * 100) %>%
  select(-c(Weight_mg, C_ug, N_ug)) %>%
  mutate(Location = "Varari")

turb2 <- rawMaya %>%
  mutate(C_N = C_ug / N_ug,
         N_percent = (N_ug/1000) / Weight_mg * 100) %>%
  select(-c(Weight_mg, C_ug, N_ug)) %>%
  mutate(CowTagID = if_else(CowTagID == "SAND SGD 3T", "VSEEP", CowTagID),
         CowTagID = if_else(CowTagID == "SAND SGD 1T", "Sand_SGD_1", CowTagID),
         CowTagID = if_else(CowTagID == "SAND SGD 2T", "Sand_SGD_2", CowTagID),
         CowTagID = if_else(CowTagID == "SAND AMB 1T", "Sand_Ambient_1", CowTagID),
         CowTagID = if_else(CowTagID == "SAND AMB 2T", "Sand_Ambient_2", CowTagID),
         CowTagID = if_else(CowTagID == "SAND AMB 3T", "Sand_Ambient_3", CowTagID),
         CowTagID = if_else(CowTagID == "REEF SGD 1T", "Reef_SGD_1", CowTagID),
         CowTagID = if_else(CowTagID == "REEF SGD 2T", "Reef_SGD_2", CowTagID),
         CowTagID = if_else(CowTagID == "REEF SGD 3T", "Reef_SGD_3", CowTagID),
         CowTagID = if_else(CowTagID == "REEF AMB 1T", "Reef_Ambient_1", CowTagID),
         CowTagID = if_else(CowTagID == "REEF AMB 2T", "Reef_Ambient_2", CowTagID),
         CowTagID = if_else(CowTagID == "REEF AMB 3T", "Reef_Ambient_3", CowTagID)) %>%
  mutate(Location = "Varari")

## remove Padina samples
AO <- str_detect(string = turb2$CowTagID, pattern = "AO")
SO <- str_detect(string = turb2$CowTagID, pattern = "SO")
turb2 <- turb2 %>%
  cbind(AO, SO) %>%
  filter(AO == FALSE, SO == FALSE) %>%
  select(-c(AO, SO))


## bind dataframes
Full_turb <- full_join(turb, turb2) %>%
  relocate(Location, .before = CowTagID)

write_csv(Full_turb, here("Data","Biogeochem","July2022","Turb_NC.csv"))


