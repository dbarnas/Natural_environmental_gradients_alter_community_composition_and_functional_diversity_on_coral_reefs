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
meta <- read_csv(here("Data","Full_metadata.csv")) %>%
  filter(CowTagID != "V13") %>%
  select(CowTagID, AlphaTag, meanRugosity)
fullchem <- read_csv(here("Data", "Biogeochem", "Nutrients_Processed_All.csv"))
chem <- fullchem %>%
  select(CowTagID, Phosphate_umolL, NN_umolL)


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
  left_join(chem) %>%
  left_join(meta)

# full trait names
Full_data <- Full_data %>%
  mutate(Morph2 = if_else(Morph2 == "Br", "Branching",
                  if_else(Morph2 == "Dig", "Digitate",
                  if_else(Morph2 == "Fol", "Foliose",
                  if_else(Morph2 == "Stol", "Stolonial",
                  if_else(Morph2 == "Mush", "Mushroom",
                  if_else(Morph2 == "Poly", "Polypoid",
                  if_else(Morph2 == "Mas", "Massive",
                  if_else(Morph2 == "Enc", "Encrusting",
                  if_else(Morph2 == "Sph", "Spherical",
                  if_else(Morph2 == "Fil", "Filamentous", Morph2)))))))))),
         Calc = if_else(Calc == "NC", "Non-calcifying",
                if_else(Calc == "AC", "Articulated",
                if_else(Calc == "Non-AC", "Non-articulated",
                if_else(Calc == "Herm", "Hermatypic", Calc)))),
         ER = if_else(ER == "Auto", "Autotroph",
              if_else(ER == "Het", "Heterotroph",
              if_else(ER == "Mix", "Mixotroph", ER))))


Full_data$Morph2 <- factor(Full_data$Morph2,
                           levels = c('Branching', 'Digitate', 'Foliose', 'Filamentous', 'Stolonial',
                                      'Mushroom', 'Polypoid', 'Cushion', 'Massive', 'Encrusting', 'Spherical'))


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
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "", y = "% Cover")

  return(plota)

}




ptplot <- function(entity, param){

  my_data <- Full_data %>%
    filter(Morph2 != "Mushroom" & Taxon_Group != "Cyanobacteria") %>%  # remove single points
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
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = "none") +
    facet_wrap(myfacet, scales = "free") +
    labs(x = paste(x.var), y = "% Cover")

  return(plota)

}




pval <- function(data = Full_data, entity, param, form = "poly"){ # poly or lm


  Param_data <- data %>%
    filter(CowTagID != "VSEEP") %>%
    mutate(entity = as.character({{entity}})) %>%
    group_by(AlphaTag,entity, {{param}}) %>%
    summarise(pCover = sum(pCover)) %>%
    filter(entity != "Cyanobacteria" & # single point in Taxon_Group
             entity != "Mushroom") %>% # single point in Morph2
    mutate(entity = if_else(entity == "Non-articulated", "NAC", entity),
           entity = if_else(entity == "Non-calcifying", "NC", entity))

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
        mutate(FTrait = if_else(FTrait == "NAC", "Non-articulated",
                        if_else(FTrait == "NC", "Non-calcifying", FTrait)))

      p_df <- p_df %>%
        rbind(temp)
      }
  return(p_df)
}

### PLOT FUNCTIONS ACROSS SEEP ###

### CALCULATE RESIDUALS OF PCOVER ~ RUGOSITY
Full_data_fe <- Full_data %>%
  group_by(CowTagID,
           AlphaTag,
           Morph2,
           Phosphate_umolL,
           meanRugosity) %>%
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
             filter(CowTagID != "VSEEP")))
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

# regression
ppt <- ptplot(Taxon_Group, Phosphate_umolL) +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  theme(strip.background = element_rect(fill = "white")) +
  labs(x = expression("CV Phosphate ("*mu*"mol/L)"), title = "Phyla")
npt <- ptplot(Taxon_Group, NN_umolL) +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  theme(strip.background = element_rect(fill = "white")) +
  labs(x = expression("CV Nitrate+Nitrite ("*mu*"mol/L)"), title = "Phyla")

#ggsave(here("Output", "PaperFigures", "taxa_nn.png"), npt, device = "png", width = 6, height = 6)

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


### MORPHOLOGY
# stacked bar
pm <- myplot(Morph2, morphpalette) + labs(fill = "Morphology")
pm2 <- myplot(Morph2, morphpalette) + theme(legend.position = "right") + labs(fill = "Morphology")
pm3 <- myplot(Morph2, morphpalette) + theme(legend.position = "none")

# regression
ppm <- ptplot(Morph2, Phosphate_umolL) +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  theme(strip.background = element_rect(fill = "white")) +
  labs(x = expression("CV Phosphate ("*mu*"mol/L)"), title = "Morphology")
npm <- ptplot(Morph2, NN_umolL) +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  theme(strip.background = element_rect(fill = "white")) +
  labs(x = expression("CV Nitrate+Nitrite ("*mu*"mol/L)"), title = "Morphology")

#ggsave(here("Output", "PaperFigures", "morphology_phos.png"), ppm, device = "png", width = 6, height = 6)

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
             filter(Morph2 == "Branching") %>%
             group_by(CowTagID, Phosphate_umolL) %>%
             summarise(pCover=sum(pCover))))

summary(lm(pCover~poly(Phosphate_umolL,2), data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             filter(Morph2 == "Massive") %>%
             group_by(CowTagID, Phosphate_umolL) %>%
             summarise(pCover=sum(pCover))))

summary(lm(pCover~poly(Phosphate_umolL,2), data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             filter(Morph2 == "Filamentous") %>%
             group_by(CowTagID, Phosphate_umolL) %>%
             summarise(pCover=sum(pCover))))

summary(lm(pCover~poly(Phosphate_umolL,2), data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             filter(Morph2 == "Foliose") %>%
             group_by(CowTagID, Phosphate_umolL) %>%
             summarise(pCover=sum(pCover))))

summary(lm(pCover~poly(Phosphate_umolL,2), data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             filter(Morph2 == "Cushion") %>%
             group_by(CowTagID, Phosphate_umolL) %>%
             summarise(pCover=sum(pCover))))

summary(lm(pCover~poly(NN_umolL,2), data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             filter(Morph2 == "Polypoid") %>%
             group_by(CowTagID, NN_umolL) %>%
             summarise(pCover=sum(pCover))))

summary(lm(pCover~poly(NN_umolL,2), data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             filter(Morph2 == "Massive") %>%
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

# regression
ppc <- ptplot(Calc, Phosphate_umolL) +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  theme(strip.background = element_rect(fill = "white")) +
  labs(x = expression("CV Phosphate ("*mu*"mol/L)"), title = "Calcification")
npc <- ptplot(Calc, NN_umolL) +
  geom_smooth(method = "lm", formula = "y~x", color = "black") +
  theme(strip.background = element_rect(fill = "white")) +
  labs(x = expression("CV Nitrate+Nitrite ("*mu*"mol/L)"), title = "Calcification")

#ggsave(here("Output", "PaperFigures", "calc_nn.png"), npc, device = "png", width = 6, height = 6)

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
             filter(Calc == "Hermatypic")))
summary(lm(pCover~Phosphate_umolL, data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             group_by(CowTagID, Calc, NN_umolL, Phosphate_umolL) %>%
             summarise(pCover=sum(pCover)) %>%
             filter(Calc == "Non-calcifying")))
summary(lm(pCover~poly(Phosphate_umolL,2), data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             group_by(CowTagID, Calc, NN_umolL, Phosphate_umolL) %>%
             summarise(pCover=sum(pCover)) %>%
             filter(Calc == "Non-articulated")))


# ENERGETIC RESOURCE
# stacked bar
per <- myplot(ER, erpalette) + labs(fill = "Energetic \nResource")

# regression
pper <- ptplot(ER, Phosphate_umolL) +
  theme(strip.background = element_rect(fill = "white")) +
  labs(x = expression("CV Phosphate ("*mu*"mol/L)"), title = "Energetic Resource")
nper <- ptplot(ER, NN_umolL) +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  theme(strip.background = element_rect(fill = "white")) +
  labs(x = expression("CV Nitrate+Nitrite ("*mu*"mol/L)"), title = "Energetic Resource")

#ggsave(here("Output", "PaperFigures", "ER_nn.png"), nper, device = "png", width = 6, height = 6)

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
             filter(ER == "Autotroph")))



# bind pvalue df
mypval <- rbind(pv1,pv2,pv3,pv4,pv5,pv6,pv7,pv8)
TraitSignif <- mypval %>%
  filter(pvalue1 < 0.07 | pvalue2 < 0.07)
write_csv(TraitSignif, here("Output", "PaperFigures", "Trait_pVal.csv"))

#######################
# FIGURE 4 BARPLOT
#######################

Figure4 <- (pt) /
  (pm) /
  (pc + per)  +
  #plot_annotation(tag_levels = 'A') +
  theme(plot.tag = element_text(size = 10))
Figure4

Figure4_ii <- (pt) /
  (pm) /
  (pc) /
  (per)  +
  #plot_annotation(tag_levels = 'A') +
  theme(plot.tag = element_text(size = 10))
Figure4_ii


 ggsave(here("Output", "PaperFigures", "Fig4_Plot_FEgroups.png"), Figure4, width = 7, height = 9)
 ggsave(here("Output", "PaperFigures", "Fig4_Plot_FEgroups_ii.png"), Figure4_ii, width = 7, height = 10)

 # ggsave(here("Output", "PaperFigures", "Plot_FEgroups_3.png"), Figure4, width = 6, height = 6)
 #
 # ggsave(here("Output", "PaperFigures", "Plot_Taxon_dist.png"), pt, width = 6, height = 3.5)
 # ggsave(here("Output", "PaperFigures", "Plot_Taxon_dist.png"), pt2, width = 6, height = 3)
 # ggsave(here("Output", "PaperFigures", "Plot_Morph_dist.png"), pm2, width = 6, height = 3.5)
 # ggsave(here("Output", "PaperFigures", "Plot_Morph_dist.png"), pm3, width = 6, height = 3)
 # ggsave(here("Output", "PaperFigures", "Plot_Calc_dist.png"), pc, width = 6, height = 3.5)
 # ggsave(here("Output", "PaperFigures", "Plot_ER_dist.png"), per, width = 6, height = 3.5)

#######################
# SUPPLEMENTAL FIGURE 4
#######################

SuppFig4 <- (ppm) / (pper + ppc) / (ppt)
SuppFig4
ggsave(here("Output", "PaperFigures", "Supp_Fig4_Trait_LM.png"),SuppFig4, width = 8, height = 12)



