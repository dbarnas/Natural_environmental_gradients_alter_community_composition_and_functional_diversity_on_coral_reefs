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




onlyParam <- chem %>%
  select(Salinity:Ammonia_umolL)
onlyParam <- colnames(onlyParam)


############################
### PREP DATA FOR FIGURE
############################

dataChem <- chem %>%
  pivot_longer(cols = Salinity:Tyrosine_Like, names_to = "Parameters", values_to = "CV") %>%
  right_join(alpha) %>%
  relocate(AlphaTag, .before = Parameters) %>%
  select(-c(CowTagID,Location)) %>%
  filter(Parameters %in% onlyParam) %>%
  filter(Parameters != "TA" &
           Parameters != "Ammonia_umolL")

tabData <- dataChem %>%
  mutate(Location = if_else(AlphaTag == "A", "Seep", "Reef")) %>%
  group_by(Location, Parameters) %>%
  mutate(CVMean = mean(CV),
         SE = std.error(CV)) %>%
  mutate(Parameters = if_else(Parameters == "NN_umolL", "Nitrate+Nitrite",
                      if_else(Parameters == "Phosphate_umolL", "Phosphate",
                      if_else(Parameters == "Silicate_umolL", "Silicate", Parameters))))


#####################################################
#####################################################

mypal <- pnw_palette("Starfish", n=19)
mypal2 <- pnw_palette("Starfish", n=7)[5]

myPhos <- dataChem %>%
  filter(Parameters == "Phosphate_umolL") %>%
  pivot_wider(names_from = Parameters, values_from = CV)

cvreef <- tabData %>%
  left_join(myPhos) %>% #will color by phosphate
  filter(Location == "Reef") %>%
  mutate(Parameters = factor(Parameters,
                             levels = c("Nitrate+Nitrite","Silicate","Phosphate","Salinity", "pH", "Temperature"))) %>%
  mutate(StatParam = "CV") %>%
  ggplot(aes(x = StatParam, y = CVMean)) +
  # raw points
  geom_point(aes(x = StatParam, y = CV, color = Phosphate_umolL),
             alpha = 0.6, position = position_jitter(width = 0.2)) +
  # primary points
  geom_point(aes(x = StatParam, y = CVMean), color = "black", size = 3) +
  geom_errorbar(aes(ymin = CVMean-SE, ymax = CVMean+SE), width = 0.1) +
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
