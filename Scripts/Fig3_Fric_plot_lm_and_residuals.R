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
# DIVERSITY MODELS
###############################

### Calculate regressions against Rugosity

ModSpRp <- lm(NbSpP ~ poly(meanRugosity,2), data=resFric) # relative species richness
ModFERp <- lm(NbFEsP ~ poly(meanRugosity,2), data=resFric) # relative entity richness
ModVol <- lm(Vol8D ~ poly(meanRugosity,2), data=resFric) # relative entity volume

#view model summary
summary(ModSpRp) # ** strong significance
summary(ModFERp) # *** strong significance
summary(ModVol) # * 0.49 weak significance



###############################
# RESIDUAL MODELS (~ SALINITY)
###############################

### Calculate regressions against Salinity

#fit model: linear relationship
# resModSpR <- lm(NbSp ~ meanRugosity, data=reg.Fric) # species richness
# resModFER <- lm(NbFEs ~ meanRugosity, data=reg.Fric) # entity richness
salModSpRp <- lm(NbSpP ~ poly(Salinity,2), data=resFric) # relative species richness
salModFERp <- lm(NbFEsP ~ poly(Salinity,2), data=resFric) # relative entity richness
salModVol <- lm(Vol8D ~ poly(Salinity,2), data=resFric) # relative entity volume

#view model summary
summary(salModSpRp)
summary(salModFERp)
summary(salModVol)


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

########################################
### Rugosity LM (Supplemental Figure 2A)
########################################

SuppFig2A <- plotFun("meanRugosity") + labs(y = "Mean rugosity") +
  plot_annotation(tag_levels = list(c('A')))
SuppFig2A
ggsave(here("Output", "PaperFigures", "SuppFig2A_regression.png"), SuppFig2A, device = "png", width = 6, height = 6)

# modeling
moddata <- funData %>%
  filter(Dependent == "meanRugosity") %>%
  pivot_wider(names_from = Dependent, values_from = Dep_Values) %>%
  pivot_wider(names_from = Parameters, values_from = Values)
mymod <- lm(data = moddata %>% rename(NN = 'Nitrate+Nitrite'),
            meanRugosity ~ Phosphate)
anova(mymod)

#check our assumptions
plot(mymod)
qqp(mymod)
resid1<-residuals(mymod)
qqp(resid1, "norm")

########################################
## plot richness and volume to rugosity
## (Supplemental Figure 2B)
########################################
SuppFig2B <- resFric %>%
  rename('%Sp Richness' = NbSpP, '%FE Richness' = NbFEsP, '%FE Volume' = Vol8D) %>%
  pivot_longer(cols = c('%Sp Richness', '%FE Richness', '%FE Volume'), names_to = "Parameters", values_to = "Values") %>%
  mutate(Parameters = factor(Parameters, levels = c('%Sp Richness', '%FE Richness', '%FE Volume'))) %>%
  ggplot(aes(x = meanRugosity, y = Values)) +#, color = NN_umolL)) +
  geom_point() +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  facet_wrap(~Parameters, scales = "fixed") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  labs(x = "Mean rugosity", y = "Relative Diversity") +
  ylim(0,80) +
  plot_annotation(tag_levels = list(c('B')))
SuppFig2B

ggsave(here("Output", "PaperFigures", "SuppFig2B_rugosity.png"), SuppFig2B, device = "png", height = 6, width = 6)


########################################
## plot richness and volume to salinity
## (Supplemental Figure 2C)
########################################
SuppFig2C <- resFric %>%
  rename('%Sp Richness' = NbSpP, '%FE Richness' = NbFEsP, '%FE Volume' = Vol8D) %>%
  pivot_longer(cols = c('%Sp Richness', '%FE Richness', '%FE Volume'), names_to = "Parameters", values_to = "Values") %>%
  mutate(Parameters = factor(Parameters, levels = c('%Sp Richness', '%FE Richness', '%FE Volume'))) %>%
  ggplot(aes(x = Salinity, y = Values)) +#, color = NN_umolL)) +
  geom_point() +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  facet_wrap(~Parameters, scales = "fixed") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45,hjust = 1),
        axis.title = element_text(size = 14)) +
  labs(x = "CV Salinity", y = "Relative Diversity") +
  ylim(0,80) +
  plot_annotation(tag_levels = list(c('C')))
SuppFig2C

#ggsave(here("Output", "PaperFigures", "SuppFig2C_salinity.png"), SuppFig2C, device = "png", height = 6, width = 6)



## WHEN SEEP IS REMOVED, DIST TO SEEP IS NO LONGER SIGNIFICANT FOR ALL PARAMETERS
## PHOSPHATE AND NN SIGNIFICANT FOR ALL (POLYNOMIAL): increase Rich and Vol with elevating NN and Phosphate, then rich and vol drop off
## SpR: increase with increasing ammonia and visible humidics
## Vol: Decrease with increasing M_C, Inc with Inc Ammonia, Inc with elevating Salinity, then vol drops again





## Other checks
###############################
# MODEL RESIDUALS ~ SGD PARAMETERS
# Supplemental Figure 1A-F
###############################

a <- plotFun("NbSpP") + labs(y = "% Sp Richness") + theme(axis.title.x = element_blank())
b <- plotFun("NbFEsP") + labs(y = "% FE Richness") + theme(axis.title.x = element_blank())
c <- plotFun("Vol8D") + labs(y = "% FE Volume")

d <- plotFun("resSpp") + labs(y = "% Sp Richness (res)") + theme(axis.title.x = element_blank())
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

plotfun <-function(data = resFric %>% left_join(alphatag), y){

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
  labs(y = "% SR") +
  theme(axis.title.x = element_blank()) +
  ylim(min = 0, max = 100)

rug_FER_plot <- plotfun(y = NbFEsP) +
  labs(y = "% FER") +
  theme(axis.title.x = element_blank()) +
  ylim(min = 0, max = 100)

rug_Vol_plot <- plotfun(y = Vol8D) +
  labs(y = "% FEV") +
  theme(axis.title.x = element_blank()) +
  ylim(min = 0, max = 100)

rugosityplot <- rug_SpR_plot + rug_FER_plot + rug_Vol_plot
rugosityplot


## Relative richness and volume residuals
rug_res_SpRp_plot <- plotfun(y = resSpp) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "% SR residuals",
       x = "")

rug_res_FERp_plot <- plotfun(y = resFEp) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "% FER residuals",
       x = expression("CV Phosphate ("*mu*"mol/L)"))

rug_res_Vol_plot <- plotfun(y = resVol) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "% FEV residuals",
       x = "")

rugosityresplot <- rug_res_SpRp_plot + rug_res_FERp_plot + rug_res_Vol_plot
rugosityresplot




divPlots <- rugosityplot / rugosityresplot +
  plot_annotation(tag_levels = 'A')
divPlots
#ggsave(here("Output", "PaperFigures", "Fig3_LM_diversity_Phosphate.png"), divPlots, height = 6, width = 10)

##########################
### NITRATES + NITRITES
##########################

## Raw richness residuals
rug_res_SpR_plot <- resFric %>%
  ggplot(aes(x = NN_umolL,
             y = resSp)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(panel.grid = element_blank()) +
  labs(y = "SR (Rugosity-normalized)",
       x = "CV of NN (umol/L)")




rug_res_FER_plot <- resFric %>%
  ggplot(aes(x = NN_umolL,
             y = resFE)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(panel.grid = element_blank()) +
  labs(y = "FER (Rugosity-normalized)",
       x = "CV of NN (umol/L)")



## Relative richness and volume residuals
rug_res_SpRp_plot <- resFric %>%
  ggplot(aes(x = NN_umolL,
             y = resSpp)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(panel.grid = element_blank()) +
  labs(y = "Relative SR (%, Rugosity-normalized)",
       x = "CV of NN (umol/L)")


rug_res_FERp_plot <- resFric %>%
  ggplot(aes(x = NN_umolL,
             y = resFEp)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(panel.grid = element_blank()) +
  labs(y = "Relative FER (%, Rugosity-normalized)",
       x = "CV of NN (umol/L)")


rug_res_Vol_plot <- resFric %>%
  ggplot(aes(x = NN_umolL,
             y = resVol)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(panel.grid = element_blank()) +
  labs(y = "Relative FE Volume (%, Rugosity-normalized)",
       x = "CV of NN (umol/L)")


(rug_res_SpR_plot + rug_res_FER_plot + rug_res_SpRp_plot)
(rug_res_FERp_plot + rug_res_Vol_plot)



# ratio of FE:Sp richness ~ distance to seep INCLUDING SEEP
# value of 1 would show that each species has its own function, and any lower values show lower functional diversity

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

resFric %>% mutate(FE_SP = NbFEs / NbSp) %>% select(CowTagID,NbFEs, NbSp, FE_SP)

pRat <- resFric %>%
  filter(CowTagID != "VSEEP") %>%
  ggplot(aes(x = Phosphate_umolL,
             y = NbFEs / NbSp)) +
  geom_point(size = 2.5) +
  #geom_smooth(method = "lm", formula = "y~x", color = "black") +
  theme_bw() +
  theme(axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title.x = element_blank()) +
  theme(panel.grid = element_blank()) +
  labs(y = "FER / SR")

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


SpFERatio <- (pRat / pRatSeep) +
   plot_annotation(tag_levels = 'A')
SpFERatio

#ggsave(here("Output", "PaperFigures", "Sp_FE_Ratio.png"), SpFERatio, device = "png", width = 6, height = 6)


summary(lm(data = resFric %>% filter(CowTagID != "VSEEP"), (NbFEs / NbSp) ~ Phosphate_umolL))
summary(lm(data = resFric, (NbFEs / NbSp) ~ Phosphate_umolL))

