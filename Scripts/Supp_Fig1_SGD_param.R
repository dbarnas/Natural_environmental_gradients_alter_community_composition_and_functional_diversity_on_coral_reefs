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
# Remove outliers / irrelevant data points
removeSite1 <- AugChemData %>%
  filter(CowTagID == "V2",
         Tide == " Low",
         Day_Night == "Day",
         Date == ymd("2021-08-08"))

removeSite2 <- AugChemData %>%
  filter(CowTagID == "C4",
         Tide =="Low",
         Day_Night == "Night",
         Date == ymd("2021-08-09"))

## Filter out redundant low tide day sample, first 'low tide' was super high
removeSite3 <- AugChemData %>%
  filter(Tide == "Low",
         Day_Night == "Day",
         Date == ymd("2021-08-06"))

removeSite5 <- AugChemData %>%
  filter(CowTagID %in% c("VSPRING", "Varari_Well", "CSPRING"))


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
  filter(Ammonia_umolL < 16)


max_data <- AugChem %>%
  group_by(Location, CowTagID) %>%
  summarise_at(vars(Salinity:Tyrosine_Like), .funs = max, na.rm = T) %>%  # select for max values
  pivot_longer(cols = Salinity:Tyrosine_Like, names_to = "Parameters", values_to = "Maximum")

min_data <- AugChem %>%
  group_by(Location, CowTagID) %>%
  summarise_at(vars(Salinity:Tyrosine_Like), .funs = min, na.rm = T) %>%  # select for max values
  pivot_longer(cols = Salinity:Tyrosine_Like, names_to = "Parameters", values_to = "Minimum")

# join max and min values and other data sets
mm_data <- full_join(max_data, min_data)

onlyParam <- AugChem %>%
  select(Salinity:Ammonia_umolL)
onlyParam <- colnames(onlyParam)


############################
### PREP DATA FOR FIGURE
############################

dataChem <- chem %>%
  filter(Season == "Dry") %>%
  right_join(alpha) %>%
  select(-c(Maximum, Minimum)) %>%
  left_join(mm_data) %>%
  select(-c(Location, lat, lon, Season, Range, MeanAll, CVAll)) %>%
  relocate(AlphaTag, .before = Parameters) %>%
  select(-CowTagID) %>%
  filter(Parameters %in% onlyParam) %>%
  filter(Parameters != "TA" &
           Parameters != "Ammonia_umolL")

tabData <- dataChem %>%
  mutate(Location = if_else(AlphaTag == "A", "Seep", "Reef")) %>%
  pivot_longer(cols = c(Minimum, Maximum, MeanSeasonal, CVSeasonal), names_to = "MeanParam", values_to = "MeanVal") %>%
  group_by(Parameters, Location, MeanParam) %>%
  mutate(AllMean = mean(MeanVal),
         SE = std.error(MeanVal)) %>%
  mutate(MeanParam = if_else(MeanParam == "MeanSeasonal", "Mean",
                             if_else(MeanParam == "Maximum", "Max",
                                     if_else(MeanParam == "Minimum", "Min",
                                             if_else(MeanParam == "CVSeasonal", "CV", MeanParam))))) %>%
  mutate(Parameters = if_else(Parameters == "NN_umolL", "Nitrate+Nitrite",
                              if_else(Parameters == "Phosphate_umolL", "Phosphate",
                                      if_else(Parameters == "Silicate_umolL", "Silicate", Parameters))))

#{r, echo = FALSE, eval = FALSE}
# mypal <- pnw_palette("Starfish", n=7)[c(4, 7)]
#
# AugChem %>%
#   mutate(Location = if_else(CowTagID == "VSEEP", "Seep", "Reef")) %>%
#   select(CowTagID:Day_Night,Location,
#          Silicate_umolL, Phosphate_umolL, NN_umolL,
#          Salinity, pH, Temperature) %>%
#   rename(Silicate = Silicate_umolL, `Nitrate+Nitrite` = NN_umolL, Phosphate = Phosphate_umolL) %>%
#   pivot_longer(cols = Silicate:Temperature, names_to = "Param", values_to = "values") %>%
#   group_by(Param, Location) %>%
#   summarise(mean = mean(values, na.rm = T),
#             min = min(values, na.rm = T),
#             max = max(values, na.rm = T),
#             se = std.error(values, na.rm = T)) %>%
#   rename(Min = min, Max = max, Mean = mean) %>%
#   pivot_longer(cols = c(Mean:Max), names_to = "metric", values_to = "values") %>%
#   mutate(se = if_else(metric == "Mean", se, 0)) %>%
#
#
#   ggplot(aes(x = metric, y = values, color = Location)) +
#   geom_point(size = 3) +
#   geom_errorbar(aes(ymin = values-se, ymax = values+se), width = 0.1) +
#   theme_classic() +
#   theme(strip.background = element_rect(fill = "white"),
#         axis.title = element_blank(),
#         axis.text.x = element_text(size = 10),
#         axis.text.y = element_text(size = 10),
#         axis.ticks.x = element_blank()) +
#   scale_color_manual(values = mypal) +
#   facet_wrap(~Param, scales = "free_y")

#####################################################
#####################################################

mypal <- pnw_palette("Starfish", n=19)
mypal2 <- pnw_palette("Starfish", n=7)[5]

myPhos <- dataChem %>%
  select(AlphaTag, Parameters,CVSeasonal) %>%
  filter(Parameters == "Phosphate_umolL") %>%
  pivot_wider(names_from = Parameters, values_from = CVSeasonal)

cvreef <- tabData %>%
  left_join(myPhos) %>% #will color by phosphate
  filter(Location == "Reef") %>%
  mutate(Parameters = factor(Parameters, levels = c("Nitrate+Nitrite","Silicate","Phosphate","Salinity", "pH", "Temperature"))) %>%
  filter(MeanParam == "CV") %>%
  ggplot(aes(x = MeanParam, y = AllMean)) +
  # raw points
  geom_point(aes(x = MeanParam, y = MeanVal, color = Phosphate_umolL),
             alpha = 0.6, position = position_jitter(width = 0.2)) +
  # primary points
  geom_point(aes(x = MeanParam, y = AllMean), color = "black", size = 3) +
  geom_errorbar(aes(ymin = AllMean-SE, ymax = AllMean+SE), width = 0.1) +
  labs(y = "Coefficient of Variation of \nphysicochemical variables", color = expression("CV Phosphate ("*mu*"mol/L)")) +
  theme_classic() +
  theme(strip.background = element_rect(fill = "white"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.ticks.x = element_blank()) +
  scale_color_gradientn(colors = mypal) +
  #scale_color_manual(values = mypal2) +
  facet_wrap(~Parameters, scales = "free_y")
#plot_annotation(tag_levels=list(c('C')))
cvreef

#ggsave(here("Output", "PaperFigures", "Fig2_Reef_CV_Biogeochem.png"), cvreef, device = "png", height = 5, width = 7)
