#### Supplemental Figure 2: Scatter plots displaying polynomial regressions of the
#### coefficient of variation of SGD parameters used for model testing against rugosity (A)
#### and polynomial regressions of relationship between diversity metrics and rugosity (B).

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
  filter(CowTagID != "V13",
         CowTagID != "VSEEP")


### Join Sp and FE and Vol4D with metadata
reg.Fric <- resFric %>%
  as_tibble() %>%
  left_join(meta) %>%
  filter(CowTagID != "VSEEP" &
           CowTagID != "V13") %>%
  mutate(meanRugosity = 1-meanRugosity)




###############################
# RESIDUAL MODELS (~ RUGOSITY)
###############################
resFric <- reg.Fric %>%
  left_join(chem)





### Supplemental Figure 2a:

parameter_list <- list("NN_umolL" = "Nitrate+Nitrite",
                       "Phosphate_umolL" = "Phosphate",
                       "Silicate_umolL" = "Silicate",
                       "Salinity" = "Salinity",
                       "pH" = "pH",
                       "Temperature" = "Temperature")
parameter_labeller <- function(variable,value){
  return(parameter_list[value])
}
supp2A <- resFric %>%
  select(-TA) %>%
  pivot_longer(cols = c(Salinity:NN_umolL), names_to = "Parameters", values_to = "Values") %>%
  mutate(Parameters = factor(Parameters, levels = c('NN_umolL', 'Phosphate_umolL', 'Silicate_umolL', 'Salinity', 'pH', 'Temperature'))) %>%
  ggplot(aes(x = Values, y = meanRugosity)) +
  geom_point() +
  geom_smooth(method = "lm", formula = "y~x", color = "black") +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "red") +
  facet_wrap(~Parameters, scales = "free_x", labeller = parameter_labeller) +
  theme_bw() +
  labs(x = "CV SGD Parameter Values", y = "% Mean rugosity") +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 10),
        axis.text.x = element_text(angle = 30, hjust=1))
supp2A




### Supplemental Figure 2b:

## plot richness and volume to rugosity
supp2B <- resFric %>%
  rename('% Taxon Richness' = NbSpP, '% FE Richness' = NbFEsP, '% FE Volume' = Vol8D) %>%
  pivot_longer(cols = c('% Taxon Richness', '% FE Richness', '% FE Volume'), names_to = "Parameters", values_to = "Values") %>%
  mutate(Parameters = factor(Parameters, levels = c('% Taxon Richness', '% FE Richness', '% FE Volume'))) %>%
  ggplot(aes(x = meanRugosity, y = Values)) +#, color = NN_umolL)) +
  geom_point() +
  geom_smooth(method = "lm", color = "black", formula = "y~log(x)") +
  facet_wrap(~Parameters, scales = "free") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 10))+
  labs(x = "Mean Rugosity",
       y = "Relative Diversity")
supp2B

#### save plots
ggsave(here("Output", "PaperFigures", "Supp_Fig2A_CV_SGD_Param_Rugosity.png"), supp2A, width = 6, height = 6, device = "png")
ggsave(here("Output", "PaperFigures", "Supp_Fig2B_Diversity_Rugosity.png"), supp2B, width = 6, height = 6, device = "png")


