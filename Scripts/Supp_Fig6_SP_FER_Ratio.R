#### Supplemental Figure 5: Relationship of functional entity richness to species richness
#### along the SGD gradient, with and without the seepage point.

### Created by Danielle Barnas
### Created on May 31, 2023

###############################
# LOAD LIBRARIES
###############################
library(tidyverse)
library(here)




###############################
# READ IN DATA
###############################
meta <- read_csv(here("Data","Full_metadata.csv"))
chem <- read_csv(here("Data","Biogeochem", "Nutrients_Processed_All.csv")) %>%
  filter(Season == "Dry") %>%
  filter(Location == "Varari",
         CowTagID != "V13") %>%
  select(CowTagID, Parameters, CVSeasonal) %>%
  pivot_wider(names_from = Parameters, values_from = CVSeasonal)
resFric <- read_csv(here("Data", "Sp_FE_Vol_res.csv")) %>%
  as_tibble() %>%
  left_join(meta) %>%
  left_join(chem) %>%
  filter(#CowTagID != "VSEEP",
    CowTagID != "V13")


###############################
# PROCESS DATA
###############################

### Calculate Ratios
resFric %>% mutate(FE_SP = NbFEs / NbSp) %>% select(CowTagID,NbFEs, NbSp, FE_SP)



###############################
# VISUALIZE
###############################

### Without seepage point
pRat <- resFric %>%
  filter(CowTagID != "VSEEP") %>%
  ggplot(aes(x = Phosphate_umolL,
             y = NbFEs / NbSp)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", formula = "y~x", color = "black") +
  theme_bw() +
  theme(axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title.x = element_blank()) +
  theme(panel.grid = element_blank()) +
  labs(y = "FER / SR")

### With seepage point
pRatSeep <- resFric %>%
  ggplot(aes(x = Phosphate_umolL,
             y = NbFEs / NbSp)) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_smooth(method = "lm", formula = "y~x", color = "black") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  theme(panel.grid = element_blank()) +
  labs(y = "FER / SR",
       x = expression("CV Phosphate ("*mu*"mol/L)"))

### Patch
SpFERatio <- (pRat / pRatSeep) +
  plot_annotation(tag_levels = 'A')

#### save plots
ggsave(here("Output", "PaperFigures", "Supp_Fig5_SP_FER_Ratio.png"), SpFERatio, width = 6, height = 6, device = "png")

