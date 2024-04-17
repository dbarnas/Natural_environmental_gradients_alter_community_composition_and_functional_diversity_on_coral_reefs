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
  filter(Season == "Dry") %>%
  filter(Location == "Varari", CowTagID != "V13") %>%
  #filter(CowTagID != "VSEEP") %>%
  select(CowTagID, Parameters, CVSeasonal) %>%
  pivot_wider(names_from = Parameters, values_from = CVSeasonal)
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

parameter_list <- list("NN_umolL" = expression(paste("Nitrate+Nitrite ("*mu*"mol/L)")),
                       "Phosphate_umolL" = expression(paste("Phosphate ("*mu*"mol/L)")))
parameter_labeller <- function(variable,value){
  return(parameter_list[value])
}

supp3a <- resFric %>%
  select(CowTagID, Salinity, Phosphate_umolL, Silicate_umolL, NN_umolL) %>%
  pivot_longer(cols = c(NN_umolL, Phosphate_umolL), names_to = "param", values_to = "values") %>%
  ggplot(aes(x = Silicate_umolL, y = values)) +
  geom_point(color = "black") +
  geom_smooth(method = "lm", color = "black") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  labs(x = expression("Silicate ("*mu*"mol/L)"),
       y = "CV parameter values") +
  facet_wrap(~param, scales = "free_y", labeller = parameter_labeller)

supp3b <- resFric %>%
  select(CowTagID, Salinity, Phosphate_umolL, Silicate_umolL, NN_umolL) %>%
  pivot_longer(cols = c(NN_umolL, Phosphate_umolL), names_to = "param", values_to = "values") %>%
  ggplot(aes(x = Salinity, y = values)) +
  geom_point(color = "black") +
  geom_smooth(method = "lm", color = "black") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  labs(x = "Salinity (psu)",
       y = "CV parameter values") +
  facet_wrap(~param, scales = "free_y", labeller = parameter_labeller)

plot_supp3 <- supp3a / supp3b #+
  #plot_annotation(tag_levels = "A")

ggsave(here("Output", "PaperFigures", "Supp_Fig3_Nutrients_SGD_Regressions.png"), plot_supp3, width = 6, height = 6, device = "png")

