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
  filter(Parameters != "Phosphate") %>%
  ggplot(aes(x = Phosphate_umolL, y = CV)) +
  geom_point() +
  labs(y = "CV of physicochemical variables (%)",
       x = "CV Phosphate (%)") +
  theme_classic() +
  theme(strip.background = element_rect(fill = "white"),
        #axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.ticks.x = element_blank()) +
  facet_wrap(~Parameters, scales = "free_y")
cvreef

# add trend line for significant correlations
lmcvreef <- cvreef + geom_smooth(method = "lm", color = "black")

ggsave(here("Output", "PaperFigures", "Supp_Fig1_Reef_CV_Biogeochem.png"), cvreef, device = "png", height = 5, width = 7)
#ggsave(here("Output", "PaperFigures", "Supp_Fig1_Reef_CV_Biogeochem_lm.png"), lmcvreef, device = "png", height = 5, width = 7)


# quick stats
lmdat <- dataChem %>% left_join(myPhos)
anova(lm(data = lmdat %>% filter(Parameters == "Salinity"), CV~Phosphate_umolL))
anova(lm(data = lmdat %>% filter(Parameters == "pH"), CV~Phosphate_umolL))
anova(lm(data = lmdat %>% filter(Parameters == "Silicate_umolL"), CV~Phosphate_umolL))
anova(lm(data = lmdat %>% filter(Parameters == "NN_umolL"), CV~Phosphate_umolL))
