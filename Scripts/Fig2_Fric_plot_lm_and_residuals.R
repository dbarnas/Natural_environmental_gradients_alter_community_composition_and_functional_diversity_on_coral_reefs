#### Linear Regressions of taxonomic and functional diversity residuals
### Taxonomic richness, functional entity richness, and functional volume
### Full community, stony coral, and macroalgae

### Created by Danielle Barnas
### Created on February 26, 2023
### Modified on February 23, 2025

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
resFric <- read_csv(here("Data", "Sp_FE_Vol_res.csv"))
resFric.coral <- read_csv(here("Data", "Sp_FE_res_coral.csv"))
resFric.ma <- read_csv(here("Data", "Sp_FE_res_macroalgae.csv"))
meta <- read_csv(here("Data", "Full_Metadata.csv"))
chem <- read_csv(here("Data","Biogeochem", "Nutrients_Processed_All.csv")) %>%
  filter(CowTagID != "V13")


### Join Sp and FE and Vol4D with metadata
resFric <- resFric %>%
  as_tibble() %>%
  left_join(meta) %>%
  filter(CowTagID != "VSEEP" &
           CowTagID != "V13") %>%
  left_join(chem)

## coral and MA

resFric.coral <- resFric.coral %>%
  as_tibble() %>%
  left_join(meta) %>%
  filter(CowTagID != "VSEEP" &
           CowTagID != "V13") %>%
  left_join(chem)

resFric.ma <- resFric.ma %>%
  as_tibble() %>%
  left_join(meta) %>%
  filter(CowTagID != "VSEEP" &
           CowTagID != "V13") %>%
  left_join(chem)

###############################
# DIVERSITY MODELS
###############################

### Calculate regressions against Rugosity

ModSpRp <- lm(NbSpP ~ poly(complexity,2), data=resFric) # relative taxonomic richness
ModFERp <- lm(NbFEsP ~ poly(complexity,2), data=resFric) # relative entity richness
ModVol <- lm(Vol8D ~ poly(complexity,2), data=resFric) # relative entity volume

#view model summary
summary(ModSpRp) # ** strong linear significance, * poly significance
summary(ModFERp) # *** strong linear significance
summary(ModVol) # * 0.039 linear significance


## coral and MA
### Calculate regressions against Rugosity

ModSpRp <- lm(NbSpP ~ poly(complexity,2), data=resFric.coral) # relative taxonomic richness
ModFERp <- lm(NbFEsP ~ poly(complexity,2), data=resFric.coral) # relative entity richness

#view model summary
summary(ModSpRp) # ns
summary(ModFERp) # * linear significance

### Calculate regressions against Rugosity

ModSpRp <- lm(NbSpP ~ poly(complexity,2), data=resFric.ma) # relative taxonomic richness
ModFERp <- lm(NbFEsP ~ poly(complexity,2), data=resFric.ma) # relative entity richness

#view model summary
summary(ModSpRp) # ** strong poly significance
summary(ModFERp) # ns

###############################
# RESIDUAL MODELS (~ NUTRIENTS)
###############################

# residual values N+N
resnnModSpRp <- lm(resSp ~ poly(NN_umolL,2), data=resFric) # relative taxonomic richness
resnnModFERp <- lm(resFEp ~ poly(NN_umolL,2), data=resFric) # relative entity richness
resnnModVol <- lm(resVol ~ poly(NN_umolL,2), data=resFric) # relative entity volume

#view model summary
summary(resnnModSpRp) # poly *
summary(resnnModFERp) # poly *
summary(resnnModVol) # poly *


###############################
# Calculate regressions against NN
## coral and MA
###############################

# Coral raw values
nnModSpRp <- lm(NbSpP ~ poly(NN_umolL,2), data=resFric.coral) # relative taxonomic richness
nnModFERp <- lm(NbFEsP ~ poly(NN_umolL,2), data=resFric.coral) # relative entity richness
#view model summary
summary(nnModSpRp) # ns
summary(nnModFERp) # ns

### Calculate residual regressions
resnnModSpRp <- lm(resSp ~ poly(NN_umolL,2), data=resFric.coral) # relative taxonomic richness
resnnModFERp <- lm(resFEp ~ poly(NN_umolL,2), data=resFric.coral) # relative entity richness
#view model summary
summary(resnnModSpRp) # * poly significance
summary(resnnModFERp) # * poly significance

# MA raw values
nnModSpRp <- lm(NbSpP ~ poly(NN_umolL,2), data=resFric.ma) # relative taxonomic richness
nnModFERp <- lm(NbFEsP ~ poly(NN_umolL,2), data=resFric.ma) # relative entity richness
#view model summary
summary(nnModSpRp) # * poly significance phosphate
summary(nnModFERp) # * poly significance NN

### Calculate residual regressions
resnnModSpRp <- lm(resSp ~ poly(NN_umolL,2), data=resFric.ma) # relative taxonomic richness
resnnModFERp <- lm(resFEp ~ poly(NN_umolL,2), data=resFric.ma) # relative entity richness
#view model summary
summary(resnnModSpRp) # ns
summary(resnnModFERp) # poly significance NN


###############################
# MODEL RESIDUALS ~ SGD PARAMETERS
###############################

#check residuals against other parameters
funData <- resFric %>%
  rename('Nitrate+Nitrite' = NN_umolL, 'Phosphate' = Phosphate_umolL,
         'Silicate' = Silicate_umolL, 'Temperature' = Temperature) %>%
  select(-TA) %>%
  pivot_longer(cols = c(Salinity:'Nitrate+Nitrite'), names_to = "Parameters", values_to = "Values") %>%
  pivot_longer(cols = c(NbSpP, NbFEsP, Vol8D, resSpp, resFEp, resVol, complexity),
               names_to = "Dependent", values_to = "Dep_Values") %>%
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




##########################
### NITRATE + NITRITE
##########################

plotfun <-function(data = resFric, y, x = NN_umolL){

  y <- enquo(y)
  x <- enquo(x)

  plota <- data %>%
    ggplot(aes(x = !!x,
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
rug_SpR_plot <- plotfun(y = NbSpP, x = NN_umolL) +
  labs(y = "% Taxon Richness") +
  theme(axis.title.x = element_blank()) +
  ylim(min = 0, max = 100)

rug_FER_plot <- plotfun(y = NbFEsP, x = NN_umolL) +
  labs(y = "% FE Richness") +
  theme(axis.title.x = element_blank()) +
  ylim(min = 0, max = 100)

rug_Vol_plot <- plotfun(y = Vol8D, x = NN_umolL) +
  labs(y = "% FE Volume") +
  theme(axis.title.x = element_blank()) +
  ylim(min = 0, max = 100)

rugosityplot <- rug_SpR_plot + rug_FER_plot + rug_Vol_plot


## Relative richness and volume residuals
rug_res_SpRp_plot <- plotfun(y = resSpp, x = NN_umolL) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "% Taxon Richness residuals",
       x = "") +
  ylim(min = -30, max = 40)

rug_res_FERp_plot <- plotfun(y = resFEp, x = NN_umolL) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "% FE Richness residuals",
       x = expression("CV N+N (%)")) +
  ylim(min = -30, max = 40)

rug_res_Vol_plot <- plotfun(y = resVol, x = NN_umolL) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "% FE Volume residuals",
       x = "") +
  ylim(min = -30, max = 40)

rugosityresplot <- rug_res_SpRp_plot + rug_res_FERp_plot + rug_res_Vol_plot


# divPlots <- rugosityplot / rugosityresplot +
#   plot_annotation(tag_levels = 'A')
# divPlots



## coral and MA (Phosphate)

plotfuncma <-function(data = resFric, y, x = NN_umolL){

  y <- enquo(y)
  x <- enquo(x)

  plota <- data %>%
    ggplot(aes(x = !!x,
               y = !!y)) +
    geom_point(size = 2.5) +
    # geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
    theme_bw() +
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 14)) +
    theme(panel.grid.major = element_blank(),
          #axis.title.y = element_blank(),
          axis.text.x = element_text(hjust = 1)) +
    labs(x = "")

  return(plota)

}

###############################
# LM Diversity ~ SGD Parameters
###############################

### N + N
test.p.1.nn <- plotfuncma(data = resFric.coral,
                           y = NbSpP,
                           x = NN_umolL) +
  labs(subtitle = "Stony Coral", x = "N + N") # NS
test.p.2.nn <- plotfuncma(data = resFric.coral,
                           y = NbFEsP,
                           x = NN_umolL) +
  labs(subtitle = "Stony Coral", x = "N + N") # NS
test.p.3.nn <- plotfuncma(data = resFric.coral,
                           y = resSpp,
                           x = NN_umolL) +
  labs(subtitle = "Stony Coral", x = "N + N") +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") # poly
test.p.4.nn <- plotfuncma(data = resFric.coral,
                           y = resFEp,
                           x = NN_umolL) +
  labs(subtitle = "Stony Coral", x = "N + N") +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") # poly

test.p.5.nn <- plotfuncma(data = resFric.ma,
                           y = NbSpP,
                           x = NN_umolL) +
  labs(subtitle = "Macroalgae", x = "N + N") # NS
test.p.6.nn <- plotfuncma(data = resFric.ma,
                           y = NbFEsP,
                           x = NN_umolL) +
  labs(subtitle = "Macroalgae", x = "N + N") +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") # poly
test.p.7.nn <- plotfuncma(data = resFric.ma,
                           y = resSpp,
                           x = NN_umolL) +
  labs(subtitle = "Macroalgae", x = "N + N") # NS
test.p.8.nn <- plotfuncma(data = resFric.ma,
                           y = resFEp,
                           x = NN_umolL) +
  labs(subtitle = "Macroalgae", x = "N + N") +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") # poly

(test.p.1.nn + test.p.2.nn) / (test.p.3.nn + test.p.4.nn)
summary(lm(data = resFric.coral, NbSpP ~ poly(NN_umolL,2))) # NS
summary(lm(data = resFric.coral, NbFEsP ~ poly(NN_umolL,2))) # NS
summary(lm(data = resFric.coral, resSpp ~ poly(NN_umolL,2))) # poly
summary(lm(data = resFric.coral, resFEp ~ poly(NN_umolL,2))) # poly

(test.p.5.nn + test.p.6.nn) / (test.p.7.nn + test.p.8.nn)
summary(lm(data = resFric.ma, NbSpP ~ poly(NN_umolL,2))) # NS
summary(lm(data = resFric.ma, NbFEsP ~ poly(NN_umolL,2))) # poly
summary(lm(data = resFric.ma, resSpp ~ poly(NN_umolL,2))) # NS
summary(lm(data = resFric.ma, resFEp ~ poly(NN_umolL,2))) # poly


### RUGOSITY
anova(lm(data = resFric.coral, NbSpP ~ complexity)) # p=0.07
anova(lm(data = resFric.coral, NbFEsP ~ complexity)) # * linear
anova(lm(data = resFric.ma, NbSpP ~ poly(complexity,2))) # * poly
anova(lm(data = resFric.ma, NbFEsP ~ poly(complexity,2))) # NS


### TEMPERATURE
summary(lm(data = resFric.coral, NbSpP ~ poly(Temperature,2))) # NS
summary(lm(data = resFric.coral, NbFEsP ~ poly(Temperature,2))) # linear
summary(lm(data = resFric.coral, resSpp ~ poly(Temperature,2))) # NS
summary(lm(data = resFric.coral, resFEp ~ poly(Temperature,2))) # linear

summary(lm(data = resFric.ma, NbSpP ~ poly(Temperature,2))) # NS
summary(lm(data = resFric.ma, NbFEsP ~ poly(Temperature,2))) # NS
summary(lm(data = resFric.ma, resSpp ~ poly(Temperature,2))) # NS
summary(lm(data = resFric.ma, resFEp ~ poly(Temperature,2))) # NS


### SILICATE
summary(lm(data = resFric.coral, NbSpP ~ poly(Silicate_umolL,2))) # NS
summary(lm(data = resFric.coral, NbFEsP ~ poly(Silicate_umolL,2))) # NS
summary(lm(data = resFric.coral, resSpp ~ poly(Silicate_umolL,2))) # NS
summary(lm(data = resFric.coral, resFEp ~ poly(Silicate_umolL,2))) # NS

summary(lm(data = resFric.ma, NbSpP ~ poly(Silicate_umolL,2))) # NS
summary(lm(data = resFric.ma, NbFEsP ~ poly(Silicate_umolL,2))) # NS
summary(lm(data = resFric.ma, resSpp ~ poly(Silicate_umolL,2))) # NS
summary(lm(data = resFric.ma, resFEp ~ poly(Silicate_umolL,2))) # NS


### PHOSPHATE
summary(lm(data = resFric.coral, NbSpP ~ poly(Phosphate_umolL,2))) # NS
summary(lm(data = resFric.coral, NbFEsP ~ poly(Phosphate_umolL,2))) # NS
summary(lm(data = resFric.coral, resSpp ~ poly(Phosphate_umolL,2))) # NS
summary(lm(data = resFric.coral, resFEp ~ poly(Phosphate_umolL,2))) # NS

summary(lm(data = resFric.ma, NbSpP ~ poly(Phosphate_umolL,2))) # poly
summary(lm(data = resFric.ma, NbFEsP ~ poly(Phosphate_umolL,2))) # NS
summary(lm(data = resFric.ma, resSpp ~ poly(Phosphate_umolL,2))) # NS
summary(lm(data = resFric.ma, resFEp ~ poly(Phosphate_umolL,2))) # NS


### pH
summary(lm(data = resFric.ma, NbSpP ~ poly(pH,2))) # NS
summary(lm(data = resFric.ma, NbFEsP ~ poly(pH,2))) # NS
summary(lm(data = resFric.ma, resSpp ~ poly(pH,2))) # NS
summary(lm(data = resFric.ma, resFEp ~ poly(pH,2))) # poly




###############################
# PLOT FIGURE 2
###############################

layout <- '
ABC
DEF
'

rug_res_SpRp_plot_titled <- rug_res_SpRp_plot + labs(subtitle = "Full Community")
rug_res_FERp_plot_titled <- rug_res_FERp_plot + labs(subtitle = "Full Community")

cma.res.plots <- (rug_res_SpRp_plot_titled + test.p.3.nn + test.p.7.nn +
                    rug_res_FERp_plot_titled + test.p.4.nn + test.p.8.nn) +
  plot_layout(design = layout)
cma.res.plots

ggsave(here("Output", "Fig2_cma_residual_plots.png"), cma.res.plots, device = "png", height = 6, width = 9)
