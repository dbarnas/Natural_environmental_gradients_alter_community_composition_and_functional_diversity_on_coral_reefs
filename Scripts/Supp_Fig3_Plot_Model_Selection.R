#### Model selection on linear and nonlinear regressions of Fric residuals and SGD parameters
#### Taxonomic Richness, FE Richness, and FE Volume vs
#### N+N, Phosphate, Silicate, Salinity, Temperature, and pH
#### Taxonomic and Functional Entity Richness for stony coral and macroalgae subgroups

### Created by Danielle Barnas
### Created on May 31, 2023
### Updated on Feb 23, 2025

###############################
# LOAD LIBRARIES
###############################
library(tidyverse)
library(here)
library(ggrepel)
library(patchwork)
library(tidytext)
library(AICcmodavg)
library(kableExtra)
library(PNWColors)


###############################
# READ IN DATA
###############################
resFric <- read_csv(here("Data", "Sp_FE_Vol_res.csv"))
resFric.coral <- read_csv(here("Data", "Sp_FE_res_coral.csv"))
resFric.ma <- read_csv(here("Data", "Sp_FE_res_macroalgae.csv"))
meta <- read_csv(here("Data", "Full_Metadata.csv"))
chem <- read_csv(here("Data","Biogeochem", "Nutrients_Processed_All.csv")) %>%
  filter(CowTagID != "V13",
         CowTagID != "VSEEP")


### Join Sp and FE and Vol4D with metadata
reg.Fric <- resFric %>%
  as_tibble() %>%
  left_join(meta) %>%
  filter(CowTagID != "VSEEP" &
           CowTagID != "V13")



# use above to justify use of resFric
resFric <- reg.Fric %>%
  left_join(chem)


#####################################################################
### MODEL SELECTION (SUPPLEMENTAL FIGURE 3)
#####################################################################

# for no interactive effect
modTable <- tibble(Y = as.character(),
                     Parameter = as.character(),
                     Reg_Type = as.character(), # linear or nonlinear
                     AICc = as.numeric(),
                     R2 = as.numeric(),
                     pVal.P = as.numeric()
                   )

myDep <- colnames(resFric %>% select(NbSpP, NbFEsP, Vol8D,
                                     resSpp, resFEp, resVol
                                     ))
mydata <- resFric %>%
  select(CowTagID,NbSpP, NbFEsP, Vol8D,
         resSpp, resFEp, resVol,
         Salinity, Temperature, pH, Phosphate_umolL,
         Silicate_umolL, NN_umolL, complexity)


for(i in myDep){
  Y <- as.character(i)
  k <- which(myDep == i) # use as multiplier for list

  for(j in 8:ncol(mydata)){
    Parameter <- colnames(mydata[j])

    # with complexity as covariate
    model1R <- lm(paste0(Y, "~", Parameter, ""), data = mydata)
    subdata1R <- as_tibble(cbind(Y,Parameter)) %>%
      mutate(Parameter = paste0(Parameter,"")) %>%
      mutate(Reg_Type = "Linear",
             AICc = AICc(model1R),
             R2 = summary(model1R)$r.squared,
             pVal.P = summary(model1R)$coefficients[8]
             )

    # with complexity as covariate
    model2R <- lm(paste0(Y, "~ poly(", Parameter, ",2)"), data = mydata)
    subdata2R <- as_tibble(cbind(Y,Parameter)) %>%
      mutate(Parameter = paste0(Parameter,"")) %>%
      mutate(Reg_Type = "Polynomial",
             AICc = AICc(model2R),
             R2 = summary(model2R)$r.squared,
             pVal.P = summary(model2R)$coefficients[12]
             )


    modTable <- modTable %>%
      rbind(subdata1R) %>%
      rbind(subdata2R)

  }
}


modelTable <- modTable %>%
  group_by(Y) %>%
  mutate(minAIC = min(AICc)) %>%
  mutate(delAICc = AICc - minAIC) %>%
  select(-c(minAIC)) %>%
  mutate(Parameter = if_else(Parameter == "NN_umolL", "Nitrate+Nitrite",
                     if_else(Parameter == "Phosphate_umolL", "Phosphate",
                     if_else(Parameter == "Silicate_umolL", "Silicate",
                     if_else(Parameter == "complexity", "Complexity", Parameter)))))


# create AIC table for paper
write_csv(modelTable, here("Output", "Model_Selection_Table.csv"))
#write_csv(modelTable_rug, here("Output", "PaperFigures","Model_Selection_Table_rugosity.csv"))

modelTable <- read_csv(here("Output", "Model_Selection_Table.csv"))
#modelTable_rug <- read_csv(here("Output", "PaperFigures","Model_Selection_Table_rugosity.csv"))

mypal <- pnw_palette("Starfish", n = 2)

modelTableData <- modelTable %>%
  group_by(Y) %>%
  mutate(Y = if_else(Y == "resSpp", "% Taxon Richness (res)",
             if_else(Y == "resFEp", "% FE Richness (res)",
             if_else(Y == "resVol", "% FE Volume (res)",
             if_else(Y == "NbSpP", "% Taxon Richness",
             if_else(Y == "NbFEsP", "% FE Richness",
             if_else(Y == "Vol8D", "% FE Volume", Y)))))),
         Y = factor(Y, levels = c('% Taxon Richness', '% FE Richness', '% FE Volume',
                                  '% Taxon Richness (res)', '% FE Richness (res)', '% FE Volume (res)')),
         Reg_Type = factor(Reg_Type, levels = c("Polynomial", "Linear")),
         Parameter = factor(Parameter,
                            levels = c("Complexity", "Phosphate", "Nitrate+Nitrite",
                                       "pH", "Salinity", "Silicate", "Temperature"))) %>%
  mutate(minAICc = min(AICc),
         deltaAICc = AICc - minAICc)

AICplot <- modelTableData %>%
  filter(Y== "% Taxon Richness" | Y == "% FE Richness") %>% #| Y== "% FE Volume") %>%
  ggplot(aes(x = deltaAICc, y = fct_reorder(Parameter, desc(Parameter)), fill = Reg_Type)) +
  geom_col(position = "dodge", color = "black") +
  facet_wrap(~Y) +
  labs(x = expression(Delta*"AICc"),
       y= "Parameters",
       fill = "Regression") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 10),
        legend.position = "top",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  scale_fill_manual(values = mypal) +
  geom_vline(xintercept = 2, linetype = "dashed", size = 0.7)


AICplotres <- modelTableData %>%
  filter(Y== "% Taxon Richness (res)" | Y == "% FE Richness (res)") %>% #| Y== "% FE Volume (res)") %>%
  ggplot(aes(x = deltaAICc, y = fct_reorder(Parameter, desc(Parameter)), fill = Reg_Type)) +
  geom_col(position = "dodge", color = "black") +
  facet_wrap(~Y) +
  labs(x = expression(Delta*"AICc"),
       y= "Parameters",
       fill = "Regression") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 10),
        legend.position = "top",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14)) +
  scale_fill_manual(values = mypal) +
  geom_vline(xintercept = 2, linetype = "dashed", size = 0.7)


AICplotAll <- AICplot / AICplotres +
  plot_annotation(title = "Full community") +
  plot_layout(guides = 'collect') & theme(legend.position = 'top')
AICplotAll



ggsave(here("Output", "SuppFig3ad_AIC_model_community.png"), AICplotAll, width = 4.5, height = 6, device = "png") # 6.5



####################################################################################
####################################################################################
## Separate Stony Corals and Macroalgae

### Join Sp and FE and Vol4D with metadata
reg.Fric.coral <- resFric.coral %>%
  as_tibble() %>%
  left_join(meta) %>%
  filter(CowTagID != "VSEEP" &
           CowTagID != "V13")
reg.Fric.ma <- resFric.ma %>%
  as_tibble() %>%
  left_join(meta) %>%
  filter(CowTagID != "VSEEP" &
           CowTagID != "V13")



# use above to justify use of resFric
resFric.coral <- reg.Fric.coral %>%
  left_join(chem)
resFric.ma <- reg.Fric.ma %>%
  left_join(chem)



# for no interactive effect
modTable <- tibble(Y = as.character(),
                   Parameter = as.character(),
                   Reg_Type = as.character(), # linear or nonlinear
                   AICc = as.numeric(),
                   R2 = as.numeric(),
                   pVal.P = as.numeric(),
                   #pVal.Rug = as.numeric()
)
modTable.coral <- modTable
modTable.ma <- modTable

myDep.cma <- colnames(resFric.coral %>% select(NbSpP, NbFEsP,
                                     resSpp, resFEp
                                     ))

mydata.coral <- resFric.coral %>%
  select(CowTagID,NbSpP, NbFEsP,
         resSpp, resFEp,
         Salinity, Temperature, pH, Phosphate_umolL,
         Silicate_umolL, NN_umolL, complexity)
mydata.ma <- resFric.ma %>%
  select(CowTagID,NbSpP, NbFEsP,
         resSpp, resFEp,
         Salinity, Temperature, pH, Phosphate_umolL,
         Silicate_umolL, NN_umolL, complexity)


for(i in myDep.cma){
  Y <- as.character(i)

  for(j in 6:ncol(mydata.coral)){ # envi parameters
    Parameter <- colnames(mydata.coral[j])

    # with complexity as covariate
    model1R <- lm(paste0(Y, "~", Parameter, ""), data = mydata.coral)
    subdata1R <- as_tibble(cbind(Y,Parameter)) %>%
      mutate(Parameter = paste0(Parameter,"")) %>%
      mutate(Reg_Type = "Linear",
             AICc = AICc(model1R),
             R2 = summary(model1R)$r.squared,
             pVal.P = summary(model1R)$coefficients[8]
      )

    # with complexity as covariate
    model2R <- lm(paste0(Y, "~ poly(", Parameter, ",2)"), data = mydata.coral)
    subdata2R <- as_tibble(cbind(Y,Parameter)) %>%
      mutate(Parameter = paste0(Parameter,"")) %>%
      mutate(Reg_Type = "Polynomial",
             AICc = AICc(model2R),
             R2 = summary(model2R)$r.squared,
             pVal.P = summary(model2R)$coefficients[12]
      )


    modTable.coral <- modTable.coral %>%
      rbind(subdata1R) %>%
      rbind(subdata2R)

  }
}

for(i in myDep.cma){
  Y <- as.character(i)

  for(j in 6:ncol(mydata.ma)){ # envi parameters
    Parameter <- colnames(mydata.ma[j])

    # with complexity as covariate
    model1R <- lm(paste0(Y, "~", Parameter, ""), data = mydata.ma)
    subdata1R <- as_tibble(cbind(Y,Parameter)) %>%
      mutate(Parameter = paste0(Parameter,"")) %>%
      mutate(Reg_Type = "Linear",
             AICc = AICc(model1R),
             R2 = summary(model1R)$r.squared,
             pVal.P = summary(model1R)$coefficients[8]
      )

    # with complexity as covariate
    model2R <- lm(paste0(Y, "~ poly(", Parameter, ",2)"), data = mydata.ma)
    subdata2R <- as_tibble(cbind(Y,Parameter)) %>%
      mutate(Parameter = paste0(Parameter,"")) %>%
      mutate(Reg_Type = "Polynomial",
             AICc = AICc(model2R),
             R2 = summary(model2R)$r.squared,
             pVal.P = summary(model2R)$coefficients[12]
      )


    modTable.ma <- modTable.ma %>%
      rbind(subdata1R) %>%
      rbind(subdata2R)

  }
}


# raw values
modelTable.coral.raw <- modTable.coral %>%
  filter(Y == "NbSpP" | Y == "NbFEsP") %>%
  group_by(Y) %>%
  mutate(minAIC = min(AICc)) %>%
  mutate(delAICc = AICc - minAIC) %>%
  select(-c(minAIC)) %>%
  mutate(Parameter = if_else(Parameter == "NN_umolL", "Nitrate+Nitrite",
                     if_else(Parameter == "Phosphate_umolL", "Phosphate",
                     if_else(Parameter == "Silicate_umolL", "Silicate",
                     if_else(Parameter == "complexity", "Complexity", Parameter)))))
modelTable.ma.raw <- modTable.ma %>%
  filter(Y == "NbSpP" | Y == "NbFEsP") %>%
  group_by(Y) %>%
  mutate(minAIC = min(AICc)) %>%
  mutate(delAICc = AICc - minAIC) %>%
  select(-c(minAIC)) %>%
  mutate(Parameter = if_else(Parameter == "NN_umolL", "Nitrate+Nitrite",
                     if_else(Parameter == "Phosphate_umolL", "Phosphate",
                     if_else(Parameter == "Silicate_umolL", "Silicate",
                     if_else(Parameter == "complexity", "Complexity", Parameter)))))
# residual values
modelTable.coral.res <- modTable.coral %>%
  filter(Y == "resSpp" | Y == "resFEp") %>%
  filter(Parameter != "complexity") %>%
  group_by(Y) %>%
  mutate(minAIC = min(AICc)) %>%
  mutate(delAICc = AICc - minAIC) %>%
  select(-c(minAIC)) %>%
  mutate(Parameter = if_else(Parameter == "NN_umolL", "Nitrate+Nitrite",
                     if_else(Parameter == "Phosphate_umolL", "Phosphate",
                     if_else(Parameter == "Silicate_umolL", "Silicate",
                     if_else(Parameter == "complexity", "Complexity", Parameter)))))
modelTable.ma.res <- modTable.ma %>%
  filter(Y == "resSpp" | Y == "resFEp") %>%
  filter(Parameter != "complexity") %>%
  group_by(Y) %>%
  mutate(minAIC = min(AICc)) %>%
  mutate(delAICc = AICc - minAIC) %>%
  select(-c(minAIC)) %>%
  mutate(Parameter = if_else(Parameter == "NN_umolL", "Nitrate+Nitrite",
                     if_else(Parameter == "Phosphate_umolL", "Phosphate",
                     if_else(Parameter == "Silicate_umolL", "Silicate",
                     if_else(Parameter == "complexity", "Complexity", Parameter)))))

modelTable.coral <- rbind(modelTable.coral.raw, modelTable.coral.res)
modelTable.ma <- rbind(modelTable.ma.raw, modelTable.ma.res)

# create AIC table for paper
write_csv(modelTable.coral, here("Output", "Model_Selection_Table_Coral.csv"))
write_csv(modelTable.ma, here("Output", "Model_Selection_Table_MA.csv"))


# quick check of similar models (AIC <=2)
modelTable.coral %>% filter(delAICc <=2) %>% arrange(Y, delAICc)
modelTable.ma %>% filter(delAICc <=2) %>% arrange(Y, delAICc)


### CORAL ###
modelTableData <- modelTable.coral %>%
  group_by(Y) %>%
  mutate(Y = if_else(Y == "resSpp", "% Taxon Richness (res)",
             if_else(Y == "resFEp", "% FE Richness (res)",
             if_else(Y == "NbSpP", "% Taxon Richness",
             if_else(Y == "NbFEsP", "% FE Richness", Y
             )))),
         Y = factor(Y, levels = c('% Taxon Richness', '% FE Richness',
                                  '% Taxon Richness (res)', '% FE Richness (res)'
                                  )),
         Parameter = if_else(Parameter == "complexity", "Complexity", Parameter),
         Reg_Type = factor(Reg_Type, levels = c("Polynomial", "Linear")),
         Parameter = factor(Parameter,
                            levels = c("Complexity", "Phosphate", "Nitrate+Nitrite",
                                       "pH", "Salinity", "Silicate", "Temperature"
                                       ))) %>%
  mutate(minAICc = min(AICc),
         deltaAICc = AICc - minAICc)

AICplot.coral <- modelTableData %>%
  filter(Y== "% Taxon Richness" | Y == "% FE Richness" #| Y== "% FE Volume"
         ) %>%
  ggplot(aes(x = deltaAICc, y = fct_reorder(Parameter, desc(Parameter)), fill = Reg_Type)) +
  geom_col(position = "dodge", color = "black") +
  facet_wrap(~Y) +
  labs(x = expression(Delta*"AICc"),
       y= "Parameters",
       fill = "Regression") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 10),
        legend.position = "top",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  scale_fill_manual(values = mypal) +
  geom_vline(xintercept = 2, linetype = "dashed", size = 0.7)


AICplotres.coral <- modelTableData %>%
  filter(Y== "% Taxon Richness (res)" | Y == "% FE Richness (res)" #| Y== "% FE Volume (res)"
         ) %>%
  ggplot(aes(x = deltaAICc, y = fct_reorder(Parameter, desc(Parameter)), fill = Reg_Type)) +
  geom_col(position = "dodge", color = "black") +
  facet_wrap(~Y) +
  labs(x = expression(Delta*"AICc"),
       y= "Parameters",
       fill = "Regression") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 10),
        legend.position = "top",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14)) +
  scale_fill_manual(values = mypal) +
  geom_vline(xintercept = 2, linetype = "dashed", size = 0.7)


AICplotAll.coral <- AICplot.coral / AICplotres.coral +
  plot_annotation(title = "Stony Coral") +
  plot_layout(guides = 'collect') & theme(legend.position = 'top')
AICplotAll.coral



### MACROALGAE ###
modelTableData <- modelTable.ma %>%
  group_by(Y) %>%
  mutate(Y = if_else(Y == "resSpp", "% Taxon Richness (res)",
             if_else(Y == "resFEp", "% FE Richness (res)",
             if_else(Y == "NbSpP", "% Taxon Richness",
             if_else(Y == "NbFEsP", "% FE Richness", Y
             )))),
         Y = factor(Y, levels = c('% Taxon Richness', '% FE Richness',
                                  '% Taxon Richness (res)', '% FE Richness (res)'
         )),
         #Parameter = if_else(Parameter == "complexity", "Complexity", Parameter),
         Reg_Type = factor(Reg_Type, levels = c("Polynomial", "Linear")),
         Parameter = factor(Parameter,
                            levels = c("Complexity", "Phosphate", "Nitrate+Nitrite",
                                       "pH", "Salinity", "Silicate", "Temperature"
                                       ))) %>%
  mutate(minAICc = min(AICc),
         deltaAICc = AICc - minAICc)

AICplot.ma <- modelTableData %>%
  filter(Y== "% Taxon Richness" | Y == "% FE Richness"
  ) %>%
  ggplot(aes(x = deltaAICc, y = fct_reorder(Parameter, desc(Parameter)), fill = Reg_Type)) +
  geom_col(position = "dodge", color = "black") +
  facet_wrap(~Y) +
  labs(x = expression(Delta*"AICc"),
       y= "Parameters",
       fill = "Regression") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 10),
        legend.position = "top",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  scale_fill_manual(values = mypal) +
  geom_vline(xintercept = 2, linetype = "dashed", size = 0.7)


AICplotres.ma <- modelTableData %>%
  filter(Y== "% Taxon Richness (res)" | Y == "% FE Richness (res)"
  ) %>%
  ggplot(aes(x = deltaAICc, y = fct_reorder(Parameter, desc(Parameter)), fill = Reg_Type)) +
  geom_col(position = "dodge", color = "black") +
  facet_wrap(~Y) +
  labs(x = expression(Delta*"AICc"),
       y= "Parameters",
       fill = "Regression") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 10),
        legend.position = "top",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14)) +
  scale_fill_manual(values = mypal) +
  geom_vline(xintercept = 2, linetype = "dashed", size = 0.7)


AICplotAll.ma <- AICplot.ma / AICplotres.ma +
  plot_annotation(title = "Macroalgae") +
  plot_layout(guides = 'collect') & theme(legend.position = 'top')
AICplotAll.ma


ggsave(here("Output", "SuppFig3be_AIC_model_coral.png"), AICplotAll.coral, width = 4.5, height = 6, device = "png")
ggsave(here("Output", "SuppFig3cf_AIC_model_ma.png"), AICplotAll.ma, width = 4.5, height = 6, device = "png")
