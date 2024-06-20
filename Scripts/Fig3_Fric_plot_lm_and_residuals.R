#### Linear Regressions of Fric and Fric residuals
### species richness, functional entity richness, and functional volume

### Created by Danielle Barnas
### Created on February 26, 2023

###############################
# LOAD LIBRARIES
###############################
library(tidyverse)
library(here)
library(ggrepel)
library(patchwork)
library(MuMIn)
library(tidytext)
library(car)



###############################
# READ IN DATA
###############################
Fric <- read_csv(here("Data", "Sp_FE_Vol.csv"))
resFric <- read_csv(here("Data", "Sp_FE_Vol_res.csv"))
meta <- read_csv(here("Data", "Full_Metadata.csv"))
chem <- read_csv(here("Data","Biogeochem", "Nutrients_Processed_All.csv")) %>%
  filter(CowTagID != "V13")


### Join Sp and FE and Vol4D with metadata
resFric <- resFric %>%
  as_tibble() %>%
  left_join(meta) %>%
  filter(CowTagID != "VSEEP" &
           CowTagID != "V13") %>%
  mutate(meanRugosity = 1-meanRugosity) %>%  # instead of lower values indicating higher structure, now high value indicate high structure
  left_join(chem)


###############################
# DIVERSITY MODELS
###############################

### Calculate regressions against Rugosity

ModSpRp <- lm(NbSpP ~ poly(meanRugosity,2), data=resFric) # relative species richness
ModFERp <- lm(NbFEsP ~ poly(meanRugosity,2), data=resFric) # relative entity richness
ModVol <- lm(Vol8D ~ poly(meanRugosity,2), data=resFric) # relative entity volume

#view model summary
summary(ModSpRp) # ** strong linear significance, * poly significance
summary(ModFERp) # *** strong linear significance
summary(ModVol) # * 0.039 linear significance



###############################
# RESIDUAL MODELS (~ SALINITY)
###############################

### Calculate regressions against Salinity

#fit model: linear relationship
salModSpRp <- lm(NbSpP ~ poly(Salinity,2), data=resFric) # relative species richness
salModFERp <- lm(NbFEsP ~ poly(Salinity,2), data=resFric) # relative entity richness
salModVol <- lm(Vol8D ~ poly(Salinity,2), data=resFric) # relative entity volume

#view model summary
summary(salModSpRp)
summary(salModFERp)
summary(salModVol)

### Calculate regressions against Phosphate

salModSpRp <- lm(NbSpP ~ poly(Phosphate_umolL,2), data=resFric) # relative species richness
salModFERp <- lm(NbFEsP ~ poly(Phosphate_umolL,2), data=resFric) # relative entity richness
salModVol <- lm(Vol8D ~ poly(Phosphate_umolL,2), data=resFric) # relative entity volume

#view model summary
summary(salModSpRp) # .0499
summary(salModFERp)# 0.061
summary(salModVol) # 0.017


###############################
# MODEL RESIDUALS ~ SGD PARAMETERS
# Supplemental Figure 2A-B
###############################

#check residuals against other parameters
funData <- resFric %>%
  rename('Nitrate+Nitrite' = NN_umolL, 'Phosphate' = Phosphate_umolL,
         'Silicate' = Silicate_umolL, 'Temperature' = Temperature) %>%
  select(-TA) %>%
  pivot_longer(cols = c(Salinity:'Nitrate+Nitrite'), names_to = "Parameters", values_to = "Values") %>%
  pivot_longer(cols = c(NbSpP, NbFEsP, Vol8D, resSpp, resFEp, resVol, meanRugosity), names_to = "Dependent", values_to = "Dep_Values") %>%
  select(CowTagID, Dependent, Dep_Values, Parameters, Values)


plotFun <- function(Y){
  data <- funData %>%
    filter(Dependent == Y) %>%
    mutate(Dependent = "DependentVar") %>%
    pivot_wider(names_from = Dependent, values_from = Dep_Values)

  plot <- data %>%
    ggplot(aes(x = Values, y = DependentVar)) +
    geom_point() +
    geom_smooth(method = "lm", formula = "y~x", color = "black") +
    geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "red") +
    facet_wrap(~Parameters, scales = "free_x") +
    theme_bw() +
    labs(x = "CV SGD Parameter Values", y = Y) +
    theme(strip.background = element_rect(fill = "white"),
          strip.text = element_text(size = 12),
          axis.text.y = element_text(size = 9),
          axis.text.x = element_text(size = 9, angle = 30, hjust = 1),
          axis.title = element_text(size = 12))

  return(plot)
}


###############################
# MODEL RESIDUALS ~ SGD PARAMETERS
# Supplemental Figure 1A-F
###############################

a <- plotFun("NbSpP") + labs(y = "% Taxon Richness") + theme(axis.title.x = element_blank())
b <- plotFun("NbFEsP") + labs(y = "% FE Richness") + theme(axis.title.x = element_blank())
c <- plotFun("Vol8D") + labs(y = "% FE Volume")

d <- plotFun("resSpp") + labs(y = "% Taxon Richness (res)") + theme(axis.title.x = element_blank())
e <- plotFun("resFEp") + labs(y = "% FE Richness (res)") + theme(axis.title.x = element_blank())
f <- plotFun("resVol") + labs(y = "% FE Volume (res)")

rawPlots <- a/b/c +
  plot_annotation(tag_levels = list(c('A', 'B', 'C')))
resPlots <- d/e/f +
  plot_annotation(tag_levels = list(c('D', 'E', 'F')))



#####################################################################
### MODELING
#####################################################################


##########################
### PHOSPHATE
##########################

plotfun <-function(data = resFric, y){

  y <- enquo(y)

  plota <- data %>%
    ggplot(aes(x = Phosphate_umolL,
               y = !!y)) +
    geom_point(size = 2.5) +
    geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
    theme_bw() +
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 14)) +
    theme(panel.grid.major = element_blank(),
          #axis.title.y = element_blank(),
          axis.text.x = element_text(hjust = 1)) +
    labs(x = "")

  return(plota)

}

## Relative non-normalized richness
rug_SpR_plot <- plotfun(y = NbSpP) +
  labs(y = "% Taxon Richness") +
  theme(axis.title.x = element_blank()) +
  ylim(min = 0, max = 100)

rug_FER_plot <- plotfun(y = NbFEsP) +
  labs(y = "% FE Richness") +
  theme(axis.title.x = element_blank()) +
  ylim(min = 0, max = 100)

rug_Vol_plot <- plotfun(y = Vol8D) +
  labs(y = "% FE Volume") +
  theme(axis.title.x = element_blank()) +
  ylim(min = 0, max = 100)

rugosityplot <- rug_SpR_plot + rug_FER_plot + rug_Vol_plot
rugosityplot


## Relative richness and volume residuals
rug_res_SpRp_plot <- plotfun(y = resSpp) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "% Taxon Richness residuals",
       x = "") +
  ylim(min = -30, max = 40)

rug_res_FERp_plot <- plotfun(y = resFEp) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "% FE Richness residuals",
       x = expression("CV Phosphate (%)")) +
  ylim(min = -30, max = 40)

rug_res_Vol_plot <- plotfun(y = resVol) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "% FE Volume residuals",
       x = "") +
  ylim(min = -30, max = 40)

rugosityresplot <- rug_res_SpRp_plot + rug_res_FERp_plot + rug_res_Vol_plot
rugosityresplot

#### SAVE FIGURE 3

divPlots <- rugosityplot / rugosityresplot +
  plot_annotation(tag_levels = 'A')
divPlots
ggsave(here("Output", "PaperFigures", "Fig3_LM_diversity_Phosphate.png"), divPlots, height = 6, width = 10)



