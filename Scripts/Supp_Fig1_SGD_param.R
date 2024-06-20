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

plotData <- tabData %>%
  ungroup() %>%
  filter(Location == "Reef") %>%
  select(AlphaTag:CV) %>%
  pivot_wider(names_from = Parameters, values_from = CV)

plot_fun <- function(param = "Salinity"){
  mydata <- plotData %>%
    select(AlphaTag, Phosphate, myParam = param)

  ggplot(data = mydata,
         aes(x = Phosphate, y = myParam)) +
    geom_point() +
    labs(y = paste("CV", param, "(%)")) +
    theme_classic()
}


a <- plot_fun(param = "Silicate") + geom_smooth(method = "lm", color = "black") +theme(axis.title.x = element_blank())
b <- plot_fun(param = "Salinity") + geom_smooth(method = "lm", color = "black") +theme(axis.title.x = element_blank())
c <- plot_fun(param = "Nitrate+Nitrite") + geom_smooth(method = "lm", color = "black") +theme(axis.title.x = element_blank())
d <- plot_fun(param = "pH") + geom_smooth(method = "lm", color = "black") +theme(axis.title.x = element_blank())
e <- plot_fun(param = "Temperature") + labs(x = "CV Phosphate (%)")

fullPlot <- (a + b + c) / (d + e + plot_spacer())
fullPlot

ggsave(here("Output", "PaperFigures", "Supp_Fig1_Reef_CV_Biogeochem.png"), fullPlot, device = "png", height = 5, width = 7)


# quick stats
lmdat <- dataChem %>% left_join(myPhos)
anova(lm(data = lmdat %>% filter(Parameters == "Salinity"), CV~Phosphate_umolL))
anova(lm(data = lmdat %>% filter(Parameters == "pH"), CV~Phosphate_umolL))
anova(lm(data = lmdat %>% filter(Parameters == "Silicate_umolL"), CV~Phosphate_umolL))
anova(lm(data = lmdat %>% filter(Parameters == "NN_umolL"), CV~Phosphate_umolL))












plotData <- tabData %>%
  ungroup() %>%
  filter(Location == "Reef") %>%
  select(AlphaTag:CV) %>%
  pivot_wider(names_from = Parameters, values_from = CV)

plot_fun <- function(param = "Salinity"){
  mydata <- plotData %>%
    select(AlphaTag, Phosphate, myParam = param)

  ggplot(data = mydata,
         aes(x = Phosphate, y = myParam)) +
  geom_point() +
  labs(y = paste("CV", param, "(%)"),
       x = "CV Phosphate (%)") +
  theme_classic()
}
