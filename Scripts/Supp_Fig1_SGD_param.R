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
alpha <- read_csv(here("Data","Full_Metadata.csv")) %>% select(CowTagID, AlphaTag)




onlyParam <- chem %>%
  select(Salinity:NN_umolL)
onlyParam <- colnames(onlyParam)


############################
### PREP DATA FOR FIGURE
############################

dataChem <- chem %>%
  pivot_longer(cols = Salinity:Tyrosine_Like, names_to = "Parameters", values_to = "CV") %>%
  right_join(alpha) %>%
  relocate(AlphaTag, .before = Parameters) %>%
  select(-c(CowTagID,Location)) %>%
  filter(Parameters %in% onlyParam)

tabData <- dataChem %>%
  mutate(Location = if_else(AlphaTag == "A", "Seep", "Reef")) %>%
  group_by(Location, Parameters) %>%
  mutate(CVMean = mean(CV),
         SE = std.error(CV)) %>%
  ungroup() %>%
  mutate(Parameters = if_else(Parameters == "NN_umolL", "Nitrate+Nitrite",
                       if_else(Parameters == "Phosphate_umolL", "Phosphate",
                       if_else(Parameters == "Silicate_umolL", "Silicate", Parameters))))

#####################################################
#####################################################

plotData <- tabData %>%
  filter(Location == "Reef") %>%
  select(AlphaTag:CV) %>%
  pivot_wider(names_from = Parameters, values_from = CV)

plot_fun <- function(param){
  mydata <- plotData %>%
    select(AlphaTag, Phosphate, `Nitrate+Nitrite`, myParam = param)

  ggplot(data = mydata,
         aes(x = `Nitrate+Nitrite`, y = myParam)) + # Phosphate # `Nitrate+Nitrite`
    geom_point() +
    labs(y = paste("CV", param, "(%)")) +
    theme_classic()
}


a <- plot_fun(param = "Silicate") + geom_smooth(method = "lm", color = "black") +theme(axis.title.x = element_blank())
b <- plot_fun(param = "Salinity") + geom_smooth(method = "lm", color = "black") +theme(axis.title.x = element_blank())
c <- plot_fun(param = "Phosphate") + geom_smooth(method = "lm", color = "black") +theme(axis.title.x = element_blank())
d <- plot_fun(param = "pH") + theme(axis.title.x = element_blank())
e <- plot_fun(param = "TA") + labs(x = "CV Nitrate+Nitrite (%)")
f <- plot_fun(param = "Temperature") +theme(axis.title.x = element_blank())  # Phosphate # Nitrate+Nitrite

layout = '
ABC
DEF
'

fullPlot <- (a + b + c + d + e + f) +
  plot_layout(design = layout)
fullPlot

# ggsave(here("Output", "Supp_Fig1_Reef_CV_Biogeochem.png"), fullPlot, device = "png", height = 5, width = 7)


# quick stats
# lmdat <- dataChem %>% left_join(myPhos)
lmdat <- dataChem %>%
  filter(Parameters == "Phosphate_umolL") %>%
  rename(Phosphate_umolL=CV) %>%
  select(-Parameters) %>%
  full_join(dataChem) %>%
  filter(AlphaTag != "A")
summary(lm(data = lmdat %>% filter(Parameters == "NN_umolL"), CV~Phosphate_umolL))

lmdat <- dataChem %>%
  filter(Parameters == "Salinity") %>%
  rename(Salinity=CV) %>%
  select(-Parameters) %>%
  full_join(dataChem) %>%
  filter(AlphaTag != "A")
summary(lm(data = lmdat %>% filter(Parameters == "Phosphate_umolL"), CV~Salinity))
summary(lm(data = lmdat %>% filter(Parameters == "NN_umolL"), CV~Salinity))

lmdat <- dataChem %>%
  filter(Parameters == "Silicate_umolL") %>%
  rename(Silicate=CV) %>%
  select(-Parameters) %>%
  full_join(dataChem) %>%
  filter(AlphaTag != "A")
summary(lm(data = lmdat %>% filter(Parameters == "Phosphate_umolL"), CV~Silicate))
summary(lm(data = lmdat %>% filter(Parameters == "NN_umolL"), CV~Silicate))

lmdat <- dataChem %>%
  filter(Parameters == "pH") %>%
  rename(pH=CV) %>%
  select(-Parameters) %>%
  full_join(dataChem) %>%
  filter(AlphaTag != "A")
summary(lm(data = lmdat %>% filter(Parameters == "Phosphate_umolL"), CV~pH))
summary(lm(data = lmdat %>% filter(Parameters == "NN_umolL"), CV~pH))

lmdat <- dataChem %>%
  filter(Parameters == "Temperature") %>%
  rename(Temperature=CV) %>%
  select(-Parameters) %>%
  full_join(dataChem) %>%
  filter(AlphaTag != "A")
summary(lm(data = lmdat %>% filter(Parameters == "NN_umolL"), CV~Temperature))


