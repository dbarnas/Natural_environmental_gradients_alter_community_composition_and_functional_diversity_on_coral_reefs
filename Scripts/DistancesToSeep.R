

#### LIBRARIES ####
library(tidyverse)
library(here)
library(geosphere)


#### READ IN DATA ####
locations <- read_csv(here("Data", "Sandwich_Locations_Final.csv"))


### CALCULATE DISTANCES TO SGD TREE

# isolate seep lat and lon at Varari
seepData <- locations %>%
  filter(Plate_Seep == 'Seep',
         Location == "Varari") %>%
  select(Location, lat, lon) %>%
  rename(lat_seep = lat,
         lon_seep = lon)


# select distinct points for each plate location to calculate distances
distData <- locations %>%
  filter(Plate_Seep == 'Plate') %>%
  filter(Location == "Varari") %>%
  select(Location, CowTagID, lat, lon) %>%
  distinct() %>%
  full_join(seepData)

# find Haversine distances for each site
Vdist <- distData %>%
  mutate(dist_to_seep_m = distHaversine(cbind(lon_seep, lat_seep), cbind(lon, lat))) %>%
  # group by sample Site Number
  group_by(CowTagID) %>%
  # choose only minimum distances
  slice(which.min(dist_to_seep_m)) %>%
  select(-c(lat_seep, lon_seep))

# isolate V13, which is upcurrent of seep point
V13dist <- Vdist %>%
  filter(CowTagID == 'V13') %>%
  mutate(dist_to_seep_m = -dist_to_seep_m) # get negative value because of opposite direction from other locations

# remove V13 from distData then rejoin with new value from above
Vdist <- Vdist %>%
  filter(CowTagID != 'V13') %>%
  rbind(V13dist)

# match seep point df with larger site df's
seepData <- seepData %>%
  mutate(CowTagID = "VSEEP",
         dist_to_seep_m = 0) %>%
  rename(lat = lat_seep,
         lon = lon_seep)

# bind Vseep to distance data
distData <- Vdist %>%
  bind_rows(seepData)

# write csv
write_csv(distData, here("Data", "Plate_Distance_to_Seep.csv"))


