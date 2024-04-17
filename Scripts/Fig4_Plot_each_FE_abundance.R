####  Plot relative abundance of individual functional groups across sites

### LOAD LIBRARIES ###

library(tidyverse)
library(here)
library(RColorBrewer)
library(patchwork)
library(PNWColors)


### READ IN DATA ###

ab.sgd <- read_csv(here("Data", "Species_Abundances_wide.csv"))
fes_traits.sgd <- read_csv(here("Data", "Distinct_FE.csv"))
spe_fes.sgd <- as.data.frame(read_csv(here("Data", "Species_FE.csv")))
alphatag <- read_csv(here("Data", "CowTag_to_AlphaTag.csv"))
meta <- read_csv(here("Data","Full_metadata.csv")) %>%
  filter(Location == "Varari",
         CowTagID != "V13") %>%
  select(CowTagID, meanRugosity)
fullchem <- read_csv(here("Data", "Biogeochem", "Nutrients_Processed_All.csv")) %>%
  filter(Location == "Varari",
         Season == "Dry")
chem <- fullchem %>%
  filter(Parameters == "Phosphate_umolL" | Parameters == "NN_umolL") %>%
  select(CowTagID, Parameters, CVSeasonal) %>%
  pivot_wider(names_from = Parameters, values_from = CVSeasonal)
fullchem <- fullchem %>%
  select(CowTagID, Parameters, CVSeasonal) %>%
  pivot_wider(names_from = Parameters, values_from = CVSeasonal)

### CREATE PALETTES FOR FIGURES ###

taxonpalette <- c("#0f4a6f", "#6aa4b0", "#D4F1F4",
                  "#98d7c2",
                  "#738fa7", "#2e8bc0", "#7391c8", "#52688f")
morphpalette <- rev(pnw_palette(name = "Moth", n = 11))
calcpalette <- c("#f9eac2","#b2d2a4", "#96ad90", "#1a4314")
#symbpalette <- c("#98d7c2", "#167d7f", "#29a0b1", "#05445e")
erpalette <- c("#fbc490", "#fbaa60", "#f67b50")
#fmpalette <- c( "#deb3ad","#de847b", "#c85250")


### JOIN DATA AND CREATE PLOT FUNCTION ###

Full_data <- ab.sgd %>%
  pivot_longer(cols = 2:ncol(ab.sgd), names_to = "Taxa", values_to = "pCover") %>%
  filter(pCover > 0) %>%
  left_join(spe_fes.sgd) %>%
  left_join(fes_traits.sgd) %>%
  left_join(alphatag) %>%
  left_join(chem) %>%
  left_join(meta)

Full_data$Morph2 <- factor(Full_data$Morph2,
                           levels = c('Br', 'Dig', 'Fol', 'Fil', 'Stol',
                                      'Mush', 'Poly', 'Cushion', 'Mas', 'Enc', 'Sph'))
# arrange alphatag factor for nutrient levels
AlphaOrder <- Full_data %>%
  left_join(chem) %>%
  arrange(Phosphate_umolL) %>%
  distinct(AlphaTag)
AlphaOrder <- AlphaOrder$AlphaTag
Full_data$AlphaTag <- factor(Full_data$AlphaTag, levels = AlphaOrder)



##############################
#  CREATE FUNCTIONS
##############################

myplot <- function(param, pal){

  my_data <- Full_data %>%
    group_by(AlphaTag,{{param}}) %>%
    summarise(pCover = sum(pCover))

  myfill = enquo(param)

  plota <- my_data %>%
    ggplot(aes(x = AlphaTag,
               y = pCover,
               fill = !!myfill)) +
    geom_col(position = "stack",
             color = "white") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = "top") +
    scale_fill_manual(values = pal) +
    labs(x = "", y = "% Cover")

  return(plota)

}




ptplot <- function(entity, param){

  my_data <- Full_data %>%
    filter(Morph2 != "Mush" & Taxon_Group != "Cyanobacteria") %>%  # remove single points
    filter(CowTagID != "VSEEP") %>%
    group_by(AlphaTag,{{entity}}, {{param}}) %>%
    summarise(pCover = sum(pCover)) %>%
    mutate(indep = round(as.numeric({{param}}),2))

  myfacet <- enquo(entity)
  x.var <- colnames(my_data[,3])

  plota <- my_data %>%
    ggplot(aes(x = indep, #!!independent,
               y = pCover,
               color = !!myfacet
               )) +
    geom_point() +
    #geom_smooth(method = "lm", formula = "y~x", color = "black") +
    #geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "red") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = "none") +
    facet_wrap(myfacet, scales = "free") +
    labs(x = paste(x.var), y = "% Cover")

  return(plota)

}

# Param_data <- Full_data %>%
#   filter(CowTagID != "VSEEP") %>%
#   mutate(entity = as.character(Morph2)) %>%
#   group_by(AlphaTag,entity,Phosphate_umolL) %>%
#   summarise(pCover = sum(pCover)) %>%
#   filter(entity != "Cyanobacteria" & entity != "Mush")


pval <- function(data = Full_data, entity, param, form = "poly"){ # poly or lm


  Param_data <- data %>%
    filter(CowTagID != "VSEEP") %>%
    mutate(entity = as.character({{entity}})) %>%
    group_by(AlphaTag,entity, {{param}}) %>%
    summarise(pCover = sum(pCover)) %>%
    filter(entity != "Cyanobacteria" & # single point in Taxon_Group
           entity != "Mush") %>% # single point in Morph2
    mutate(entity = if_else(entity == "Non-AC", "NAC", entity))

  distinctFT <- unique(Param_data$entity)

  p_df <- tibble(FTrait = as.character(),
                 Parameter = as.character(),
                 pvalue1 = as.numeric(),
                 pvalue2 = as.numeric(),
                 r_squared = as.numeric(),
                 adj_r_squared = as.numeric())

  for(i in 1:length(distinctFT)){

    mydata <- Param_data %>%
      filter(entity %in% distinctFT[i]) %>%
      pivot_wider(names_from = entity, values_from = pCover)

      Parameter = colnames(mydata)[2]
      FTrait = distinctFT[i]
      if(form == "poly"){
      mod <- lm(paste(FTrait, "~ poly(", Parameter, ",2)"), data = mydata)
      pvalue1 <- summary(mod)[4]$coefficients[11]
      pvalue2 <- summary(mod)[4]$coefficients[12]
      r_squared <- summary(mod)[8]$r.squared
      adj_r_squared <- summary(mod)[9]$adj.r.squared
      } else if(form == "lm"){
        mod <- lm(paste(FTrait, "~", Parameter), data = mydata)
        pvalue1 <- summary(mod)$coefficients[8]
        pvalue2 <- NA
        r_squared <- summary(mod)$r.squared
        adj_r_squared <- summary(mod)$adj.r.squared
      }

      temp <- as_tibble(cbind(FTrait, Parameter,
                              pvalue1, pvalue2, r_squared, adj_r_squared)) %>%
        mutate(pvalue1 = as.numeric(pvalue1),
               pvalue2 = as.numeric(pvalue2),
               r_squared = as.numeric(r_squared),
               adj_r_squared = as.numeric(adj_r_squared)) %>%
        mutate(FTrait = if_else(FTrait == "NAC", "Non-AC", FTrait))

      p_df <- p_df %>%
        rbind(temp)
      }
  return(p_df)
}

### PLOT FUNCTIONS ACROSS SEEP ###

### CALCULATE RESIDUALS OF PCOVER ~ RUGOSITY
Full_data_fe <- Full_data %>%
  group_by(CowTagID, AlphaTag,
           #FE, Morph2, Calc, ER,
           Morph2,
           Phosphate_umolL, meanRugosity) %>%
  summarise(pCover = sum(pCover))
Full_data_fe %>%
  ggplot(aes(x = meanRugosity, y = pCover)) +
  geom_point() +
  facet_wrap(~Morph2, scales = "free_y")
mymod <- lm(data = Full_data_fe %>%  filter(CowTagID != "VSEEP"),
   pCover ~ meanRugosity)
summary(mymod)

### FUNCTIONAL ENTITIES
anova(lm(pCover ~ poly(Phosphate_umolL,2)*FE, data = Full_data %>%
             filter(CowTagID != "VSEEP"))) # using all abundance data at each location
### SUMMARY
# Does an individual entity change its relative abundance along the phosphate gradient? Yes!
## Chlorophyta, Fil, NC, Auto p=0.005 (no interaction...not sure what this means then)
## Cnidaria, Mas, Herm, Mix p=0.003 (no interaction...not sure what this means then)
# Does an individual entity change along the NN gradient? Yes!
## Cnidaria, Mas, Herm, Mix p=5.4e-5 (interactio bw NN and FE)
## Phaeophyta, Fol, Non-AC, Auto p=0.002 (no interaction...not sure what this means then)
## Chlorophyta, Fil, NC, Auto p=0.01 (interaction bw NN and FE)
### ANOVA
# NN:FE interaction p<0.0008, F=2.08, DF=34
# P:FE interaction p<0.002, F=1.99, DF=34

anova(lm(pCover ~ poly(Phosphate_umolL,2)*FE, data = Full_data)) # including the seep
### ANOVA
# NN:FE interaction p=0.001, F=2.03, DF=34
# P:FE no interaction


### TAXA
# stacked bar
pt <- myplot(Taxon_Group, taxonpalette) + labs(fill = "Phyla")
pt2 <- myplot(Taxon_Group, taxonpalette) + theme(legend.position = "none")
pt
# regression
ppt <- ptplot(Taxon_Group, Phosphate_umolL) +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  theme(strip.background = element_rect(fill = "white")) +
  labs(x = expression("CV Phosphate ("*mu*"mol/L)"), title = "Phyla")
npt <- ptplot(Taxon_Group, NN_umolL) +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  theme(strip.background = element_rect(fill = "white")) +
  labs(x = expression("CV Nitrate+Nitrite ("*mu*"mol/L)"), title = "Phyla")
npt
ppt
ggsave(here("Output", "PaperFigures", "taxa_nn.png"), npt, device = "png", width = 6, height = 6)
# pvalues
pv1 <- pval(entity = Taxon_Group, param = Phosphate_umolL)
pv2 <- pval(entity = Taxon_Group, param = NN_umolL)
# anova
anova(lm(pCover~poly(NN_umolL,2)*Taxon_Group, data = Full_data %>%
           filter(CowTagID != "VSEEP") %>%
           group_by(CowTagID, Taxon_Group, NN_umolL, Phosphate_umolL) %>%
           summarise(pCover=sum(pCover))))
anova(lm(pCover~poly(Phosphate_umolL,2)*Taxon_Group, data = Full_data %>%
           filter(CowTagID != "VSEEP") %>%
           group_by(CowTagID, Taxon_Group, NN_umolL, Phosphate_umolL) %>%
           summarise(pCover=sum(pCover))))
## no interaction of NN:FE or P:FE
# summary
summary(lm(pCover~poly(NN_umolL,2)*Taxon_Group, data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             group_by(CowTagID, Taxon_Group, NN_umolL, Phosphate_umolL) %>%
             summarise(pCover=sum(pCover))))
summary(lm(pCover~poly(Phosphate_umolL,2)*Taxon_Group, data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             group_by(CowTagID, Taxon_Group, NN_umolL, Phosphate_umolL) %>%
             summarise(pCover=sum(pCover))))
summary(lm(pCover~poly(NN_umolL,2), data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             group_by(CowTagID, Taxon_Group, NN_umolL, Phosphate_umolL) %>%
             summarise(pCover=sum(pCover)) %>%
             filter(Taxon_Group == "Turf")))
Full_data %>%
  filter(CowTagID!= "VSEEP") %>%
  group_by(CowTagID, Taxon_Group, NN_umolL, Phosphate_umolL) %>%
  summarise(pCover=sum(pCover)) %>%
  ggplot(aes(x = NN_umolL, y = pCover)) +geom_point()+geom_smooth(method = "lm", formula="y~poly(x,2)") + facet_wrap(~Taxon_Group)


### MORPHOLOGY
# stacked bar
pm <- myplot(Morph2, morphpalette) + labs(fill = "Morphology")
pm2 <- myplot(Morph2, morphpalette) + theme(legend.position = "right") + labs(fill = "Morphology")
pm3 <- myplot(Morph2, morphpalette) + theme(legend.position = "none")
pm
# regression
ppm <- ptplot(Morph2, Phosphate_umolL) +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  theme(strip.background = element_rect(fill = "white")) +
  labs(x = expression("CV Phosphate ("*mu*"mol/L)"), title = "Morphology")
npm <- ptplot(Morph2, NN_umolL) +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  theme(strip.background = element_rect(fill = "white")) +
  labs(x = expression("CV Nitrate+Nitrite ("*mu*"mol/L)"), title = "Morphology")
npm
ppm
ggsave(here("Output", "PaperFigures", "morphology_phos.png"), ppm, device = "png", width = 6, height = 6)
ggsave(here("Output", "PaperFigures", "morphology_nn.png"), npm, device = "png", width = 6, height = 6)
# pvalues
pv3 <- pval(entity = Morph2, param = Phosphate_umolL)
pv4 <- pval(entity = Morph2, param = NN_umolL)
# anova
anova(lm(pCover~poly(NN_umolL,2)*Morph2, data = Full_data %>%
           filter(CowTagID!= "VSEEP") %>%
           group_by(CowTagID, Morph2, NN_umolL, Phosphate_umolL) %>%
           summarise(pCover=sum(pCover))))
anova(lm(pCover~poly(Phosphate_umolL,2)*Morph2, data = Full_data %>%
           filter(CowTagID!= "VSEEP") %>%
           group_by(CowTagID, Morph2, NN_umolL, Phosphate_umolL) %>%
           summarise(pCover=sum(pCover))))
## NN:FE p<0.04, F=1.83, DF=18
## P:FE p<0.03, F=1.87, DF=18
# summary
summary(lm(pCover~poly(NN_umolL,2)*Morph2, data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             group_by(CowTagID, Morph2, NN_umolL, Phosphate_umolL) %>%
             summarise(pCover=sum(pCover))))
summary(lm(pCover~poly(Phosphate_umolL,2)*Morph2, data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             group_by(CowTagID, Morph2, NN_umolL, Phosphate_umolL) %>%
             summarise(pCover=sum(pCover))))
# by morphology trait
summary(lm(pCover~poly(Phosphate_umolL,2), data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             filter(Morph2 == "Br") %>%
             group_by(CowTagID, Phosphate_umolL) %>%
             summarise(pCover=sum(pCover))))

summary(lm(pCover~poly(Phosphate_umolL,2), data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             filter(Morph2 == "Mas") %>%
             group_by(CowTagID, Phosphate_umolL) %>%
             summarise(pCover=sum(pCover))))

summary(lm(pCover~poly(Phosphate_umolL,2), data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             filter(Morph2 == "Fil") %>%
             group_by(CowTagID, Phosphate_umolL) %>%
             summarise(pCover=sum(pCover))))

summary(lm(pCover~poly(Phosphate_umolL,2), data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             filter(Morph2 == "Fol") %>%
             group_by(CowTagID, Phosphate_umolL) %>%
             summarise(pCover=sum(pCover))))

summary(lm(pCover~poly(Phosphate_umolL,2), data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             filter(Morph2 == "Cushion") %>%
             group_by(CowTagID, Phosphate_umolL) %>%
             summarise(pCover=sum(pCover))))

summary(lm(pCover~poly(NN_umolL,2), data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             filter(Morph2 == "Poly") %>%
             group_by(CowTagID, NN_umolL) %>%
             summarise(pCover=sum(pCover))))

summary(lm(pCover~poly(NN_umolL,2), data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             filter(Morph2 == "Mas") %>%
             group_by(CowTagID, NN_umolL) %>%
             summarise(pCover=sum(pCover))))

summary(lm(pCover~poly(NN_umolL,2), data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             filter(Morph2 == "Cushion") %>%
             group_by(CowTagID, NN_umolL) %>%
             summarise(pCover=sum(pCover))))


# CALCIFICATION
# stacked bar
pc <- myplot(Calc, calcpalette) + labs(fill = "Calcification")
pc
# regression
ppc <- ptplot(Calc, Phosphate_umolL) +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  theme(strip.background = element_rect(fill = "white")) +
  labs(x = expression("CV Phosphate ("*mu*"mol/L)"), title = "Calcification")
npc <- ptplot(Calc, NN_umolL) +
  geom_smooth(method = "lm", formula = "y~x", color = "black") +
  theme(strip.background = element_rect(fill = "white")) +
  labs(x = expression("CV Nitrate+Nitrite ("*mu*"mol/L)"), title = "Calcification")
#npc2 <- ptplot(Calc, NN_umolL) + geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black")
#pper + nper
npc
ppc
ggsave(here("Output", "PaperFigures", "calc_nn.png"), npc, device = "png", width = 6, height = 6)

# pvalues
pv5 <- pval(entity = Calc, param = Phosphate_umolL)
pv6 <- pval(entity = Calc, param = NN_umolL)
# anova
anova(lm(pCover~NN_umolL*Calc, data = Full_data %>%
           filter(CowTagID!= "VSEEP") %>%
           group_by(CowTagID, Calc, NN_umolL, Phosphate_umolL) %>%
           summarise(pCover=sum(pCover))))
anova(lm(pCover~poly(Phosphate_umolL,2)*Calc, data = Full_data %>%
           filter(CowTagID!= "VSEEP") %>%
           group_by(CowTagID, Calc, NN_umolL, Phosphate_umolL) %>%
           summarise(pCover=sum(pCover))))
## NN:FE p<0.01, F=3.14, DF=6
## P:FE no interaction
# summary
summary(lm(pCover~poly(NN_umolL,2)*Calc, data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             group_by(CowTagID, Calc, NN_umolL, Phosphate_umolL) %>%
             summarise(pCover=sum(pCover))))
summary(lm(pCover~NN_umolL*Calc, data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             group_by(CowTagID, Calc, NN_umolL, Phosphate_umolL) %>%
             summarise(pCover=sum(pCover))))
summary(lm(pCover~poly(Phosphate_umolL,2)*Calc, data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             group_by(CowTagID, Calc, NN_umolL, Phosphate_umolL) %>%
             summarise(pCover=sum(pCover))))
summary(lm(pCover~NN_umolL, data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             group_by(CowTagID, Calc, NN_umolL, Phosphate_umolL) %>%
             summarise(pCover=sum(pCover)) %>%
             filter(Calc == "Herm")))
summary(lm(pCover~Phosphate_umolL, data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             group_by(CowTagID, Calc, NN_umolL, Phosphate_umolL) %>%
             summarise(pCover=sum(pCover)) %>%
             filter(Calc == "NC")))
summary(lm(pCover~poly(Phosphate_umolL,2), data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             group_by(CowTagID, Calc, NN_umolL, Phosphate_umolL) %>%
             summarise(pCover=sum(pCover)) %>%
             filter(Calc == "Non-AC")))


# ENERGETIC RESOURCE
# stacked bar
per <- myplot(ER, erpalette) + labs(fill = "Energetic \nResource")
per
# regression
pper <- ptplot(ER, Phosphate_umolL) +
  #geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  #geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "red") +
  theme(strip.background = element_rect(fill = "white")) +
  labs(x = expression("CV Phosphate ("*mu*"mol/L)"), title = "Energetic Resource")
nper <- ptplot(ER, NN_umolL) +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  theme(strip.background = element_rect(fill = "white")) +
  labs(x = expression("CV Nitrate+Nitrite ("*mu*"mol/L)"), title = "Energetic Resource")
nper
pper
ggsave(here("Output", "PaperFigures", "ER_nn.png"), nper, device = "png", width = 6, height = 6)
# pvalues
pv7 <- pval(entity = ER, param = Phosphate_umolL)
pv8 <- pval(entity = ER, param = NN_umolL)
# anova
anova(lm(pCover~poly(NN_umolL,2)*ER, data = Full_data %>%
           filter(CowTagID!= "VSEEP") %>%
           group_by(CowTagID, ER, NN_umolL, Phosphate_umolL) %>%
           summarise(pCover=sum(pCover))))
anova(lm(pCover~poly(Phosphate_umolL,2)*ER, data = Full_data %>%
           filter(CowTagID!= "VSEEP") %>%
           group_by(CowTagID, ER, NN_umolL, Phosphate_umolL) %>%
           summarise(pCover=sum(pCover))))
## NN:FE p<2e-6, F=3.67, DF=18
## P:FE p<5e-7, F=3.87, DF=18
# summary
summary(lm(pCover~NN_umolL*ER, data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             group_by(CowTagID, ER, NN_umolL, Phosphate_umolL) %>%
             summarise(pCover=sum(pCover))))
summary(lm(pCover~poly(Phosphate_umolL,2)*ER, data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             group_by(CowTagID, ER, NN_umolL, Phosphate_umolL) %>%
             summarise(pCover=sum(pCover))))
summary(lm(pCover~poly(NN_umolL,2), data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             group_by(CowTagID, ER, NN_umolL, Phosphate_umolL) %>%
             summarise(pCover=sum(pCover)) %>%
             filter(ER == "Auto")))



# bind pvalue df
mypval <- rbind(pv1,pv2,pv3,pv4,pv5,pv6,pv7,pv8)
TraitSignif <- mypval %>%
  filter(pvalue1 < 0.07 | pvalue2 < 0.07)
write_csv(TraitSignif, here("Output", "PaperFigures", "Trait_pVal.csv"))

#######################
# FIGURE 4 BARPLOT
#######################

Figure4 <- (pm) /
  (per + pc) /
  (pt) +
  plot_annotation(tag_levels = 'A') +
  theme(plot.tag = element_text(size = 10))
Figure4


 ggsave(here("Output", "PaperFigures", "Fig5_Plot_FEgroups.png"), Figure4, width = 8, height = 10)
 ggsave(here("Output", "PaperFigures", "Plot_FEgroups_3.png"), Figure4, width = 6, height = 6)

 ggsave(here("Output", "PaperFigures", "Plot_Taxon_dist.png"), pt, width = 6, height = 3.5)
 ggsave(here("Output", "PaperFigures", "Plot_Taxon_dist.png"), pt2, width = 6, height = 3)
 ggsave(here("Output", "PaperFigures", "Plot_Morph_dist.png"), pm2, width = 6, height = 3.5)
 ggsave(here("Output", "PaperFigures", "Plot_Morph_dist.png"), pm3, width = 6, height = 3)
 ggsave(here("Output", "PaperFigures", "Plot_Calc_dist.png"), pc, width = 6, height = 3.5)
 ggsave(here("Output", "PaperFigures", "Plot_ER_dist.png"), per, width = 6, height = 3.5)

#######################
# SUPPLEMENTAL FIGURE 4
#######################

SuppFig4 <- (ppm) / (pper + ppc) / (ppt)
SuppFig4
ggsave(here("Output", "PaperFigures", "SuppFig4_Trait_LM.png"),SuppFig4, width = 8, height = 12)


####################################################
####################################################










 ### DISTANCE ###


 FE_data <- Full_data %>%
   group_by(CowTagID, AlphaTag, FE, Taxon_Group, Morph2, Calc, ER) %>%
   summarise(pCover = sum(pCover)) %>%
   ungroup() %>%
   left_join(dist) %>%
   left_join(chem)

 morphd <- FE_data %>%
   select(dist_to_seep_m, NN_umolL, Phosphate_umolL, Morph2, pCover) %>%
   group_by(dist_to_seep_m, Morph2) %>%
   mutate(pCover = sum(pCover)) %>%
   distinct()
 anova(lm(pCover ~ dist_to_seep_m*Morph2, data = morphd))

 taxond <- FE_data %>%
   select(dist_to_seep_m, NN_umolL, Phosphate_umolL, Taxon_Group, pCover) %>%
   group_by(dist_to_seep_m, Taxon_Group) %>%
   mutate(pCover = sum(pCover)) %>%
   distinct()
 anova(lm(pCover ~ dist_to_seep_m*Taxon_Group, data = taxond))

 calcd <- FE_data %>%
   select(dist_to_seep_m, NN_umolL, Phosphate_umolL, Calc, pCover) %>%
   group_by(dist_to_seep_m, Calc) %>%
   mutate(pCover = sum(pCover)) %>%
   distinct()
 anova(lm(pCover ~ dist_to_seep_m*Calc, data = calcd))


 erd <- FE_data %>%
   select(dist_to_seep_m, NN_umolL, Phosphate_umolL, ER, pCover) %>%
   group_by(dist_to_seep_m, ER) %>%
   mutate(pCover = sum(pCover)) %>%
   distinct()
 anova(lm(pCover ~ dist_to_seep_m*ER, data = erd))


 anova(lm(pCover ~ dist_to_seep_m*FE, data = Full_data %>% left_join(dist)))

 Full_data %>%
   left_join(dist) %>%
   left_join(chem) %>%
   group_by(dist_to_seep_m, FE) %>%
   summarise(pCover = sum(pCover)) %>%
   ggplot(aes(x = dist_to_seep_m, y = pCover, color = FE)) +
   geom_point()+
   geom_smooth(method = "lm", formula = "y~x", color = "black") +
   geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "red") +
   theme_bw() +
   facet_wrap(~FE, scales = "free_y") +
   theme(legend.position = "none")


 # nutrients relative to distance
 summary(lm(data = Full_data %>% filter(CowTagID != "VSEEP") %>% left_join(dist), NN_umolL ~ poly(dist_to_seep_m,2)))
 Full_data %>%
   filter(CowTagID != "VSEEP") %>%
   left_join(dist) %>%
   select(AlphaTag, Phosphate_umolL, NN_umolL, dist_to_seep_m) %>%
   rename('Nitrate+Nitrite' = NN_umolL,
          'Phosphate' = Phosphate_umolL) %>%
   pivot_longer(cols = c('Phosphate', 'Nitrate+Nitrite'), names_to = "Parameters", values_to = "Values") %>%
   ggplot(aes(x = dist_to_seep_m, y = Values)) +
   geom_point(aes(color = Parameters), show.legend = FALSE) +
   geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
   facet_wrap(~Parameters, scales = "free") +
   theme_bw() +
   labs(x = "Distance to seep (m)", y = "CV of Nutrient Values (umol/L)")


 ### PHOSPHATE ###



 FE_data <- Full_data %>%
   group_by(CowTagID, AlphaTag, FE, Taxon_Group, Morph2, Calc, ER) %>%
   summarise(pCover = sum(pCover)) %>%
   ungroup() %>%
   left_join(chem)

 morphd <- FE_data %>%
   select(Phosphate_umolL, Morph2, pCover) %>%
   group_by(Phosphate_umolL, Morph2) %>%
   mutate(pCover = sum(pCover)) %>%
   distinct()
 anova(lm(pCover ~ Phosphate_umolL*Morph2, data = morphd))

 taxond <- FE_data %>%
   select(Phosphate_umolL, Taxon_Group, pCover) %>%
   group_by(Phosphate_umolL, Taxon_Group) %>%
   mutate(pCover = sum(pCover)) %>%
   distinct()
 anova(lm(pCover ~ Phosphate_umolL*Taxon_Group, data = taxond)) # taxa NS ~ phosphate

 calcd <- FE_data %>%
   select(Phosphate_umolL, Calc, pCover) %>%
   group_by(Phosphate_umolL, Calc) %>%
   mutate(pCover = sum(pCover)) %>%
   distinct()
 anova(lm(pCover ~ Phosphate_umolL*Calc, data = calcd))


 erd <- FE_data %>%
   select(Phosphate_umolL, ER, pCover) %>%
   group_by(Phosphate_umolL, ER) %>%
   mutate(pCover = sum(pCover)) %>%
   distinct()
 anova(lm(pCover ~ Phosphate_umolL*ER, data = erd)) # er NS ~ phosphate


 anova(lm(pCover ~ Phosphate_umolL*FE, data = Full_data %>% left_join(chem)))


 Full_data %>%
   left_join(chem) %>%
   group_by(Phosphate_umolL, FE) %>%
   summarise(pCover = sum(pCover)) %>%
   ggplot(aes(x = Phosphate_umolL, y = pCover, color = FE)) +
   geom_point()+
   geom_smooth(method = "lm", formula = "y~x", color = "black") +
   geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "red") +
   theme_bw() +
   facet_wrap(~FE, scales = "free") +
   theme(legend.position = "none")



# nutrients relative to silicate and salinity
summary(lm(data = Full_data %>% filter(CowTagID != "VSEEP") %>% left_join(fullchem), Phosphate_umolL ~ Silicate_umolL))
Full_data %>%
  filter(CowTagID != "VSEEP") %>%
  left_join(fullchem) %>%
  select(AlphaTag, Phosphate_umolL, NN_umolL, Silicate_umolL) %>%
  rename('Nitrate+Nitrite' = NN_umolL,
         'Phosphate' = Phosphate_umolL) %>%
  pivot_longer(cols = c('Phosphate', 'Nitrate+Nitrite'), names_to = "Parameters", values_to = "Values") %>%
  ggplot(aes(x = Silicate_umolL, y = Values)) +
  geom_point(aes(color = Parameters), show.legend = FALSE) +
  geom_smooth(method = "lm", formula = "y~x", color = "black") +
  facet_wrap(~Parameters, scales = "free") +
  theme_bw() +
  labs(x = "CV of Silicate (umol/L)", y = "CV of Nutrient Values (umol/L)")




View(Full_data %>%
  filter(CowTagID!= "VSEEP") %>%
  group_by(FE) %>%
  summarise(sumCover = sum(pCover),
            meanCover=mean(pCover),
            sdCover=sd(pCover),
            seCover = sdCover / sqrt(nrow(.))))

Full_data %>%
  filter(FE == "Chlorophyta,Fil,NC,Auto") %>%
  distinct(Taxa)
Full_data %>%
  filter(FE == "Phaeophyta,Br,NC,Auto") %>%
  group_by(Taxa) %>%
  summarise(sumCover = sum(pCover),
            meanCover=mean(pCover),
            sdCover=sd(pCover),
            seCover = sdCover / sqrt(nrow(.)))
Full_data %>%
  filter(FE == "Cnidaria,Mas,Herm,Mix") %>%
  group_by(Taxa) %>%
  summarise(sumCover = sum(pCover),
            meanCover=mean(pCover),
            sdCover=sd(pCover),
            seCover = sdCover / sqrt(nrow(.)))

Full_data %>%
  group_by(CowTagID) %>%
  distinct(FE) %>%
  ungroup() %>%
  dplyr::count(FE) %>%
  arrange(desc(n))

# View(meta %>%
#   filter(Location == "Varari") %>%
#     drop_na(lat) %>%
#     select(CowTagID, Sand) %>%
#     mutate(NotSand = (100-Sand)/100) %>%
#     right_join(Full_data) %>%
#     filter(CowTagID != "VSEEP") %>%
#     select(CowTagID:AlphaTag))

#### FUNCTIONAL REDUNDANCY ACROSS GRADIENT WITHIN EACH PLOT
View(Full_data %>%
  group_by(CowTagID, AlphaTag) %>%
  dplyr::count(FE) %>%
  mutate(totalFE = sum(n)) %>%
  filter(n > 1) %>%
  mutate(totalRedundantFE = sum(n)) %>%
  distinct(CowTagID, AlphaTag, totalRedundantFE, totalFE) %>%
  mutate(totalRedundancy = totalRedundantFE / totalFE * 100) %>%
  arrange(desc(totalRedundancy)))

View(Full_data %>%
  group_by(CowTagID, AlphaTag) %>%
  dplyr::count(FE) %>%
  mutate(FER = 1) %>% # get FE richness
  mutate(totalFE = sum(FER)) %>%
  filter(n > 1) %>% # only keep redundant FE
  mutate(totalRedundantFE = sum(FER)) %>%
  distinct(CowTagID, AlphaTag, totalRedundantFE, totalFE) %>%
  mutate(relativeRedundantFE = totalRedundantFE / totalFE * 100) %>%
  arrange(desc(relativeRedundantFE)))


#############################
# MACROALGAE ALONG GRADIENT
#############################
macroalg <- c("Turf", "Rhodophyta", "Chlorophyta", "Phaeophyta")

Full_data %>%
  filter(CowTagID != "VSEEP") %>%
  mutate(Taxon_Group = if_else(Taxon_Group %in% macroalg, "Macroalgae", Taxon_Group)) %>%
  group_by(CowTagID, AlphaTag, Phosphate_umolL, Taxon_Group) %>%
  summarise(pCover = sum(pCover)) %>%
  filter(Taxon_Group != "Cyanobacteria") %>%
  ggplot(aes(x = Phosphate_umolL, y = pCover)) +
  geom_point() +
  #geom_smooth(method = "lm", formula = "y~poly(x,2)") +
  facet_wrap(~Taxon_Group, scales = "free") +
  theme_bw()+
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white"),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 15)) +
  labs(x = "CV Phosphate", y = "% Cover")

MAabundance <- Full_data %>%
  filter(CowTagID != "VSEEP") %>%
  mutate(Taxon_Group = if_else(Taxon_Group %in% macroalg, "Macroalgae", Taxon_Group)) %>%
  group_by(CowTagID, AlphaTag, Phosphate_umolL, Taxon_Group) %>%
  summarise(pCover = sum(pCover))

summary(lm(data = MAabundance, pCover ~ poly(Phosphate_umolL,2)))

### RICHNESS
Full_data %>%
  filter(CowTagID != "VSEEP") %>%
  mutate(Taxon_Group = if_else(Taxon_Group %in% macroalg, "Macroalgae", Taxon_Group)) %>%
  group_by(CowTagID, AlphaTag, Phosphate_umolL, Taxa, Taxon_Group) %>%
  summarise(pCover = sum(pCover)) %>%
  ungroup() %>%
  count(CowTagID, AlphaTag, Phosphate_umolL, Taxa) %>%
  group_by(CowTagID, AlphaTag, Phosphate_umolL) %>%
  summarise(richness = sum(n)) %>%
  ggplot(aes(x = Phosphate_umolL, y = richness)) +
  geom_point(color = "red", size = 5) +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  theme_bw()+
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18)) +
  labs(x = "CV Phosphate", y = "Macroalgal Richness")


