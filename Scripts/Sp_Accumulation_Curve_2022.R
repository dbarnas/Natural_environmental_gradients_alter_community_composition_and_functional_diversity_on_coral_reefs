### Species Accumulation Curve from June 2022 ###
### Created by Danielle Barnas
### Created on November 13, 2021

library(tidyverse)
library(here)


### READ IN DATA
survey <- read_csv(here("Data","Surveys","Species_Composition_2022.csv"))

### ANALYSIS

## Total quads by species/taxa counts
## species accumulation curve shows how many new species you get by doing more surveys

# Remove abiotic from "Taxa"
abiotic <- survey %>%
  filter(Taxa %in% c("Sand", "Rubble", "Bare Rock", "Bare Rock - exposed"))

Vquad <- survey %>%
  anti_join(abiotic) %>%
  filter(Location == "Varari",
         CowTagID != "V13") %>%
  select(CowTagID, Taxa)


Vquad_dup <- duplicated(Vquad$Taxa) # returns logical vector showing which taxa are duplicated
Vquad <- Vquad %>%
  cbind(Vquad_dup) %>% # bind T/F vector relating which taxa are duplicated
  arrange(CowTagID)
Vquad_CT <- as_tibble(unique(Vquad$CowTagID)) %>% rename(CowTagID = value) # vector of vquad cowtagid's
Vquad <- Vquad %>%
  filter(Vquad_dup == FALSE) %>%  # only keep taxa that are NOT duplicated
  full_join(Vquad_CT) %>%
  mutate(Vquad_dup = if_else(is.na(Vquad_dup), TRUE, FALSE))



# create column "n" containing a value of 1 next to each species and 0 next to NAs
# n = 0 for survey locations with no new taxa
Vquad_0 <- Vquad %>%
  filter(Vquad_dup == TRUE) %>%
  mutate(n = 0) %>%
  select(-Vquad_dup)
# n = 1 for all new taxa at each survey location, then bind Vquad_0
Vquad <- Vquad %>%
  filter(Vquad_dup == FALSE) %>%
  group_by(CowTagID) %>%
  dplyr::count(Taxa) %>%
  ungroup() %>%
  rbind(Vquad_0)


# add total species per quadrat survey
Vquad <- Vquad %>%
  group_by(CowTagID) %>%
  summarise(sp.sum = sum(n)) %>%
  ungroup() %>%
  arrange(desc(sp.sum))


# add column to enumerate each quadrat survey
Vquad$row.num=seq.int(nrow(Vquad)) # add column of row numbers aka quad number

# sequentially add species richness as number of surveys increases
Vquad<-Vquad %>%
  mutate(sp_accumulation = cumsum(sp.sum))

plot1 <- Vquad %>%
  ggplot(aes(y = sp_accumulation, x = row.num)) +
  geom_point() +
  labs(x = "Total Quadrat Surveys",
       y = "Species Richness",
       title = "Varari Species Accumulation Curve",
       subtitle = "June 2022") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.background = element_rect(fill = "#e8b4b8"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  geom_smooth(color = "black", fill = "white")
plot1

# ggsave(here("Output", "Species_Accumulation_Varari.png"), plot1, height = 5, width = 5, device = "png")
 ggsave(here("Output", "PaperFigures", "Species_Accumulation_Varari.png"), plot1, height = 5, width = 5, device = "png")


#############################################################################################
#############################################################################################
