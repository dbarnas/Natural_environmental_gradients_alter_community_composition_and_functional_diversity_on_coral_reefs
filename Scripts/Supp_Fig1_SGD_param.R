### Supplemental Figure 1: SGD along the reef


############################
### LOAD LIBRARIES
############################
library(tidyverse)
library(here)
library(kableExtra)
library(curl)
library(lubridate)
library(plotrix)
library(PNWColors)
library(patchwork)




############################
### READ IN DATA
############################
chem <- read_csv(here("Data", "Biogeochem", "Nutrients_Processed_All.csv"))
alpha <- read_csv(here("Data", "CowTag_to_AlphaTag.csv"))
AugChemData<-read_csv(curl('https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/August2021/Allbiogeochemdata_QC.csv'))


############################
### CLEAN RAW AUGUST DATA FOR MAX MIN
############################
AugChemData <- AugChemData %>%
  filter(Location == "Varari") # focal site

# Remove outliers / irrelevant data points
removeSite1 <- AugChemData %>%
  filter(CowTagID == "V2",
         Tide == " Low",
         Day_Night == "Day",
         Date == ymd("2021-08-08"))

## Filter out redundant low tide day sample, first 'low tide' was super high
removeSite2 <- AugChemData %>%
  filter(Tide == "Low",
         Day_Night == "Day",
         Date == ymd("2021-08-06"))

removeSite3 <- AugChemData %>%
  filter(CowTagID %in% c("VSPRING", "Varari_Well"))


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
  filter(Ammonia_umolL < 16)


max_data <- AugChem %>%
  group_by(CowTagID) %>%
  summarise_at(vars(Salinity:Tyrosine_Like), .funs = max, na.rm = T) %>%  # select for max values
  pivot_longer(cols = Salinity:Tyrosine_Like, names_to = "Parameters", values_to = "Maximum")

min_data <- AugChem %>%
  group_by(CowTagID) %>%
  summarise_at(vars(Salinity:Tyrosine_Like), .funs = min, na.rm = T) %>%  # select for max values
  pivot_longer(cols = Salinity:Tyrosine_Like, names_to = "Parameters", values_to = "Minimum")

mean_data <- AugChem %>%
  group_by(CowTagID) %>%
  summarise_at(vars(Salinity:Tyrosine_Like), .funs = mean, na.rm = T) %>%  # select for max values
  pivot_longer(cols = Salinity:Tyrosine_Like, names_to = "Parameters", values_to = "Mean")


# join max and min values and other data sets
mm_data <- full_join(max_data, min_data) %>%
  full_join(mean_data)

onlyParam <- AugChem %>%
  select(Salinity:Ammonia_umolL)
onlyParam <- colnames(onlyParam)


############################
### PREP DATA FOR FIGURE
############################

dataChem <- chem %>%
  pivot_longer(cols = Salinity:Tyrosine_Like, names_to = "Parameters", values_to = "CV") %>%
  right_join(alpha) %>%
  left_join(mm_data) %>%
  relocate(AlphaTag, .before = Parameters) %>%
  select(-c(CowTagID,Location)) %>%
  filter(Parameters %in% onlyParam) %>%
  filter(Parameters != "TA" &
           Parameters != "Ammonia_umolL")

tabData <- dataChem %>%
  mutate(Location = if_else(AlphaTag == "A", "Seep", "Reef")) %>%
  pivot_longer(cols = c(Minimum, Maximum, Mean, CV), names_to = "StatParam", values_to = "StatVal") %>%
  group_by(Location, Parameters, StatParam) %>%
  mutate(AllMean = mean(StatVal),
         SE = std.error(StatVal)) %>%
  mutate(StatParam = if_else(StatParam == "Maximum", "Max",
                     if_else(StatParam == "Minimum", "Min", StatParam))) %>%
  mutate(Parameters = if_else(Parameters == "NN_umolL", "Nitrate+Nitrite",
                      if_else(Parameters == "Phosphate_umolL", "Phosphate",
                      if_else(Parameters == "Silicate_umolL", "Silicate", Parameters))))


#####################################################
#####################################################

mypal <- pnw_palette("Starfish", n=19)
mypal2 <- pnw_palette("Starfish", n=7)[5]

myPhos <- dataChem %>%
  select(AlphaTag, Parameters, CV) %>%
  filter(Parameters == "Phosphate_umolL") %>%
  pivot_wider(names_from = Parameters, values_from = CV)

cvreef <- tabData %>%
  left_join(myPhos) %>% #will color by phosphate
  filter(Location == "Reef") %>%
  mutate(Parameters = factor(Parameters,
                             levels = c("Nitrate+Nitrite","Silicate","Phosphate","Salinity", "pH", "Temperature"))) %>%
  filter(StatParam == "CV") %>%
  ggplot(aes(x = StatParam, y = AllMean)) +
  # raw points
  geom_point(aes(x = StatParam, y = StatVal, color = Phosphate_umolL),
             alpha = 0.6, position = position_jitter(width = 0.2)) +
  # primary points
  geom_point(aes(x = StatParam, y = AllMean), color = "black", size = 3) +
  geom_errorbar(aes(ymin = AllMean-SE, ymax = AllMean+SE), width = 0.1) +
  labs(y = "Coefficient of Variation of \nphysicochemical variables",
       color = expression("CV Phosphate ("*mu*"mol/L)")) +
  theme_classic() +
  theme(strip.background = element_rect(fill = "white"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.ticks.x = element_blank()) +
  scale_color_gradientn(colors = mypal) +
  facet_wrap(~Parameters, scales = "free_y")
cvreef

ggsave(here("Output", "PaperFigures", "Supp_Fig1_Reef_CV_Biogeochem.png"), cvreef, device = "png", height = 5, width = 7)
