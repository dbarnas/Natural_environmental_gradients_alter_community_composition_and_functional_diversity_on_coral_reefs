#### Supplemental Figure 3: Scatter plots displaying linear regressions of the
#### coefficient of variation of silicate and salinity on CV nitrate+nitrite and phosphate.

### Created by Danielle Barnas
### Created on May 31, 2023

###############################
# LOAD LIBRARIES
###############################
library(tidyverse)
library(here)
library(ggrepel)
library(patchwork)
library(MuMIn)
library(tidytext)



###############################
# READ IN DATA
###############################
Fric <- read_csv(here("Data", "Sp_FE_Vol.csv"))
resFric <- read_csv(here("Data", "Sp_FE_Vol_res.csv"))
meta <- read_csv(here("Data", "Full_Metadata.csv"))
chem <- read_csv(here("Data","Biogeochem", "Nutrients_Processed_All.csv")) %>%
  filter(CowTagID != "V13")
alphatag <- read_csv(here("Data", "CowTag_to_AlphaTag.csv"))


### Join Sp and FE and Vol4D with metadata
resFric <- resFric %>%
  as_tibble() %>%
  left_join(meta) %>%
  filter(CowTagID != "VSEEP" &
           CowTagID != "V13") %>%
  mutate(meanRugosity = 1-meanRugosity) %>%  # instead of lower values indicating higher structure, now high value indicate high structure
  left_join(chem)


###############################
# VISUALIZATION
###############################

parameter_list <- list("NN_umolL" = expression(paste("Nitrate+Nitrite (%)")),
                       "Phosphate_umolL" = expression(paste("Phosphate (%)")))
parameter_labeller <- function(variable,value){
  return(parameter_list[value])
}

supp3a <- resFric %>%
  select(CowTagID, Salinity, Phosphate_umolL, Silicate_umolL, NN_umolL) %>%
  ggplot(aes(x = Silicate_umolL, y = NN_umolL)) +
  geom_point(color = "black") +
  geom_smooth(method = "lm", color = "black") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  labs(x = "CV Silicate (%)",
       y = "CV Nitrate+Nitrite (%)")

supp3b <- resFric %>%
  select(CowTagID, Salinity, Phosphate_umolL, Silicate_umolL, NN_umolL) %>%
  ggplot(aes(x = Silicate_umolL, y = Phosphate_umolL)) +
  geom_point(color = "black") +
  geom_smooth(method = "lm", color = "black") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  labs(x = "CV Silicate (%)",
       y = "CV Phosphate (%)")

supp3c <- resFric %>%
  select(CowTagID, Salinity, Phosphate_umolL, Silicate_umolL, NN_umolL) %>%
  ggplot(aes(x = Salinity, y = NN_umolL)) +
  geom_point(color = "black") +
  geom_smooth(method = "lm", color = "black") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  labs(x = "CV Salinity (%)",
       y = "CV Nitrate+Nitrite (%)")

supp3d <- resFric %>%
  select(CowTagID, Salinity, Phosphate_umolL, Silicate_umolL, NN_umolL) %>%
  ggplot(aes(x = Salinity, y = Phosphate_umolL)) +
  geom_point(color = "black") +
  geom_smooth(method = "lm", color = "black") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  labs(x = "CV Salinity (%)",
       y = "CV Phosphate (%)")

plot_supp3 <- (supp3a + supp3b) /
              (supp3c + supp3d) +
  plot_annotation(tag_levels = "A")
plot_supp3

ggsave(here("Output", "PaperFigures", "Supp_Fig3_Nutrients_SGD_Regressions.png"), plot_supp3, width = 6, height = 6, device = "png")

