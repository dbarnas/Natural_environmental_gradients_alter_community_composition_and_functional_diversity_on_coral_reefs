#### Supplemental Figure 6: Stacked bar plot of proportional substrate at each
#### survey site, ordered by increasing CV phosphate.

### Created by Danielle Barnas
### Created on May 31, 2023

############################
### LOAD LIBRARIES
############################
library(tidyverse)
library(here)
library(PNWColors)



############################
### READ IN DATA
############################
sub <- read_csv(here("Data", "Surveys", "Substrate_2022.csv"))
alpha <- read_csv(here("Data", "CowTag_to_AlphaTag.csv"))
chem <- read_csv(here("Data", "Biogeochem", "Nutrients_Processed_All.csv"))


############################
### PROCESS DATA FOR FIGURE
############################
sub <- sub %>%
  right_join(alpha) %>%
  select(AlphaTag, LiveCoral:Sand) %>%
  pivot_longer(cols = c(LiveCoral:Sand), names_to = "substrate", values_to = "values") %>%
  group_by(AlphaTag) %>%
  mutate(total = sum(values),
         pcover = values/total*100) %>%
  select(AlphaTag, substrate, pcover)

chem <- chem %>%
  right_join(alpha) %>%
  filter(Season == "Dry") %>%
  select(AlphaTag, Parameters, CVSeasonal) %>%
  filter(Parameters == "Phosphate_umolL") %>%
  arrange(CVSeasonal)
myorder <- chem$AlphaTag

sub <- sub %>%
  mutate(AlphaTag = factor(AlphaTag, levels = myorder),
         substrate = if_else(substrate == "LiveCoral", "Live Coral",
                             if_else(substrate == "DeadCoral", "Dead Coral", substrate)))



############################
### VISUALIZATION
############################
mypal <- pnw_palette("Starfish", n=7)[c(3,6,4,7)]

subPlot <- sub %>%
  ggplot(aes(x = AlphaTag, y = pcover)) +
  geom_col(aes(fill = substrate), position = "stack") +
  scale_fill_manual(values = mypal) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 13)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Survey Locations", y = "% Cover", fill = "Substrate")

subPlot


ggsave(here("Output", "PaperFigures", "Supp_Fig5_Substrate.png"), subPlot, device = "png", width = 6, height = 6)
