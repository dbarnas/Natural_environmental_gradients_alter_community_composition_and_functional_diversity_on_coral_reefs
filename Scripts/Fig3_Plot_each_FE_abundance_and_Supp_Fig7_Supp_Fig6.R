####  Plot relative abundance of individual functional groups across sites


### Created by Danielle Barnas
### Created on February 26, 2023
### Modified on February 23, 2025

###############################
# LOAD LIBRARIES
###############################

library(tidyverse)
library(here)
library(RColorBrewer)
library(patchwork)
library(PNWColors)


###############################
# READ IN DATA
###############################
ab.sgd <- read_csv(here("Data", "Species_Abundances_wide.csv"))
spe_fes.sgd <- as.data.frame(read_csv(here("Data", "Distinct_Taxa_FEtraits.csv")))
meta <- read_csv(here("Data","Full_metadata.csv")) %>%
  filter(CowTagID != "V13") %>%
  select(CowTagID, AlphaTag, complexity)
fullchem <- read_csv(here("Data", "Biogeochem", "Nutrients_Processed_All.csv"))
chem <- fullchem %>%
  select(CowTagID, Phosphate_umolL, NN_umolL)


###############################
# CREATE PALETTES FOR FIGURES
###############################
taxonpalette <- c("#0f4a6f", "#6aa4b0", "#D4F1F4",
                  "#98d7c2",
                  "#738fa7", "#2e8bc0", "#7391c8", "#52688f")
morphpalette <- rev(pnw_palette(name = "Moth", n = 11))
calcpalette <- c("#f9eac2","#b2d2a4", "#96ad90", "#1a4314")
erpalette <- c("#fbc490", "#fbaa60", "#f67b50")


######################################
# JOIN DATA AND CREATE PLOT FUNCTION
######################################
Full_data <- ab.sgd %>%
  pivot_longer(cols = 2:ncol(ab.sgd), names_to = "Taxa", values_to = "pCover") %>%
  filter(pCover > 0) %>%
  left_join(spe_fes.sgd) %>%
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
  arrange(NN_umolL) %>%
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
               y = pCover)) +
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
           NN_umolL,
           complexity) %>%
  summarise(pCover = sum(pCover))

Full_data_fe %>%
  ggplot(aes(x = complexity, y = pCover)) +
  geom_point() +
  facet_wrap(~Morph2, scales = "free_y")

mymod <- lm(data = Full_data_fe %>%  filter(CowTagID != "VSEEP"),
   pCover ~ complexity*Morph2)

summary(mymod)

### FUNCTIONAL ENTITIES
anova(lm(pCover ~ poly(NN_umolL,2)*FE, data = Full_data %>%
             filter(CowTagID != "VSEEP")))


anova(lm(pCover ~ poly(NN_umolL,2)*FE, data = Full_data)) # including the seep


### TAXA
# stacked bar
pt <- myplot(Taxon_Group, taxonpalette) + labs(fill = "Phyla")
pt2 <- myplot(Taxon_Group, taxonpalette) + theme(legend.position = "none")

# regression
ppt <- ptplot(Taxon_Group, Phosphate_umolL) +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  theme(strip.background = element_rect(fill = "white")) +
  labs(x = "CV Phosphate (%)", title = "Phyla")
npt <- ptplot(Taxon_Group, NN_umolL) +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  theme(strip.background = element_rect(fill = "white")) +
  labs(x = "CV Nitrate+Nitrite (%)", title = "Phyla")
ppt_2 <- ptplot(Taxon_Group, NN_umolL) +
  theme(strip.background = element_rect(fill = "white")) +
  labs(x = "CV Nitrate+Nitrite (%)", title = "Phyla")

# ggsave(here("Output", "taxa_nn.png"), npt, device = "png", width = 6, height = 6)
# ggsave(here("Output", "taxa_nn_noLM.png"), ppt_2, device = "png", width = 6, height = 6)

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

# summary
summary(lm(pCover~poly(NN_umolL,2)*Taxon_Group, data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             group_by(CowTagID, Taxon_Group, NN_umolL, Phosphate_umolL) %>%
             summarise(pCover=sum(pCover))))
summary(lm(pCover~poly(Phosphate_umolL,2)*Taxon_Group, data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             group_by(CowTagID, Taxon_Group, NN_umolL, Phosphate_umolL) %>%
             summarise(pCover=sum(pCover))))
# Turf~(NN,2), check alone:
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
  labs(x = "CV Phosphate (%)", title = "Morphology")
npm <- ptplot(Morph2, NN_umolL) +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  theme(strip.background = element_rect(fill = "white")) +
  labs(x = "CV Nitrate+Nitrite (%)", title = "Morphology")
ppm_2 <- ptplot(Morph2, NN_umolL) +
  theme(strip.background = element_rect(fill = "white")) +
  labs(x = "CV Nitrate+Nitrite (%)", title = "Morphology")

# ggsave(here("Output", "morphology_phos.png"), npm, device = "png", width = 6, height = 6)
# ggsave(here("Output", "morphology_phos_noLM.png"), ppm_2, device = "png", width = 6, height = 6)

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
summary(lm(pCover~poly(NN_umolL,2), data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             filter(Morph2 == "Branching") %>%
             group_by(CowTagID, NN_umolL) %>%
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
summary(lm(pCover~poly(NN_umolL,2), data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             filter(Morph2 == "Filamentous") %>%
             group_by(CowTagID, NN_umolL) %>%
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
  labs(x = "CV Phosphate (%)", title = "Calcification")
npc <- ptplot(Calc, NN_umolL) +
  geom_smooth(method = "lm", formula = "y~x", color = "black") +
  theme(strip.background = element_rect(fill = "white")) +
  labs(x = "CV Nitrate+Nitrite (%)", title = "Calcification")
ppc_2 <- ptplot(Calc, NN_umolL) +
  theme(strip.background = element_rect(fill = "white")) +
  labs(x = "CV Phosphate (%)", title = "Calcification")

# ggsave(here("Output", "calc_nn.png"), npc, device = "png", width = 6, height = 6)
# ggsave(here("Output", "calc_nn_noLM.png"), ppc_2, device = "png", width = 6, height = 6)

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


# summary
summary(lm(pCover~poly(NN_umolL,2)*Calc, data = Full_data %>%
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
summary(lm(pCover~NN_umolL, data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             group_by(CowTagID, Calc, NN_umolL, Phosphate_umolL) %>%
             summarise(pCover=sum(pCover)) %>%
             filter(Calc == "Hermatypic")))
# get logged stats for a logarithmic relationship (levels off at high sgd)
summary(lm(pCover~log(NN_umolL), data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             group_by(CowTagID, Calc, NN_umolL, Phosphate_umolL) %>%
             summarise(pCover=sum(pCover)) %>%
             filter(Calc == "Hermatypic")))


summary(lm(pCover~poly(Phosphate_umolL,2), data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             group_by(CowTagID, Calc, NN_umolL, Phosphate_umolL) %>%
             summarise(pCover=sum(pCover)) %>%
             filter(Calc == "Non-calcifying")))
summary(lm(pCover~NN_umolL, data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             group_by(CowTagID, Calc, NN_umolL, Phosphate_umolL) %>%
             summarise(pCover=sum(pCover)) %>%
             filter(Calc == "Non-calcifying")))

summary(lm(pCover~poly(Phosphate_umolL,2), data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             group_by(CowTagID, Calc, NN_umolL, Phosphate_umolL) %>%
             summarise(pCover=sum(pCover)) %>%
             filter(Calc == "Non-articulated")))


# TROPHIC GROUP
# stacked bar
per <- myplot(ER, erpalette) + labs(fill = "Trophic \nGroup")

# regression
pper <- ptplot(ER, Phosphate_umolL) +
  theme(strip.background = element_rect(fill = "white")) +
  labs(x = "CV Phosphate (%)", title = "Trophic Group")
nper <- ptplot(ER, NN_umolL) +
  # geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  theme(strip.background = element_rect(fill = "white")) +
  labs(x = "CV Nitrate+Nitrite (%)", title = "Trophic Group")

# ggsave(here("Output", "ER_nn_noLM.png"), nper, device = "png", width = 6, height = 6)

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

# summary
summary(lm(pCover~NN_umolL*ER, data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             group_by(CowTagID, ER, NN_umolL, Phosphate_umolL) %>%
             summarise(pCover=sum(pCover))))
summary(lm(pCover~poly(Phosphate_umolL,2)*ER, data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             group_by(CowTagID, ER, NN_umolL, Phosphate_umolL) %>%
             summarise(pCover=sum(pCover))))
# check (NS)
summary(lm(pCover~poly(NN_umolL,2), data = Full_data %>%
             filter(CowTagID!= "VSEEP") %>%
             group_by(CowTagID, ER, NN_umolL, Phosphate_umolL) %>%
             summarise(pCover=sum(pCover)) %>%
             filter(ER == "Autotroph")))



# bind pvalue df
mypval <- rbind(pv1,pv2,pv3,pv4,pv5,pv6,pv7,pv8)
TraitSignif <- mypval %>%
  filter(pvalue1 < 0.07 | pvalue2 < 0.07)
# write_csv(TraitSignif, here("Output", "Trait_pVal.csv"))

#######################
# FIGURE 3 BARPLOT
#######################

Figure3 <- (pt) /
  (pm) /
  (pc + per)  +
  #plot_annotation(tag_levels = 'A') +
  theme(plot.tag = element_text(size = 10))
Figure3



 ggsave(here("Output", "Fig3_Plot_FEgroups_NN.png"), Figure3, width = 7, height = 9)



 # ggsave(here("Output", "PaperFigures", "Plot_FEgroups_3.png"), Figure4, width = 6, height = 6)
 #
 # ggsave(here("Output", "PaperFigures", "Plot_Taxon_dist.png"), pt, width = 6, height = 3.5)
 # ggsave(here("Output", "PaperFigures", "Plot_Taxon_dist.png"), pt2, width = 6, height = 3)
 # ggsave(here("Output", "PaperFigures", "Plot_Morph_dist.png"), pm2, width = 6, height = 3.5)
 # ggsave(here("Output", "PaperFigures", "Plot_Morph_dist.png"), pm3, width = 6, height = 3)
 # ggsave(here("Output", "PaperFigures", "Plot_Calc_dist.png"), pc, width = 6, height = 3.5)
 # ggsave(here("Output", "PaperFigures", "Plot_ER_dist.png"), per, width = 6, height = 3.5)

#######################
# SUPPLEMENTAL FIGURE 7
#######################

SuppFig7 <- (npt) / (npm) / (npc + nper)
   # plot_annotation(tag_levels = 'A')
SuppFig7
# ggsave(here("Output", "Supp_Fig7_Trait_LM_NN.png"),SuppFig7, width = 8, height = 12)

SuppFig7_noregression <- (ppt_2) / (ppm_2) / (ppc_2 + pper)
  # plot_annotation(tag_levels = 'A')
SuppFig7_noregression
ggsave(here("Output", "Supp_Fig7_Trait_noLM_NN.png"),SuppFig7_noregression, width = 8, height = 12)



#######################
# SUPPLEMENTAL FIGURE 6
#######################
### DOMINANT TAXA ALONG SGD NUTRIENT GRADIENT

taxa_data_sum <- Full_data %>%
  filter(CowTagID != "VSEEP") %>%
  group_by(CowTagID, Taxa, Phosphate_umolL, NN_umolL) %>%
  summarise(pCover = sum(pCover, na.rm = TRUE))

### get significance values across taxa
summary(lm(data = taxa_data_sum, pCover ~ poly(NN_umolL,2)*Taxa))

### get specific r2
summary(lm(data = Full_data %>% filter(Taxa == "Montipora grisea"), pCover ~ poly(NN_umolL,2)))

summary(lm(data = Full_data %>% filter(Taxa == "Turf"), pCover ~ poly(NN_umolL,2))) # 0.065

summary(lm(data = Full_data %>% filter(Taxa == "Turbinaria ornata"), pCover ~ poly(NN_umolL,2))) # poly

summary(lm(data = Full_data %>% filter(Taxa == "Padina boryana"), pCover ~ NN_umolL)) # lm

summary(lm(data = Full_data %>% filter(Taxa == "Porites rus"), pCover ~ poly(NN_umolL,2))) # 0.061


# N+N
taxa.nn.poly <- taxa_data_sum %>%
  filter(Taxa == "Padina boryana" | # Taxa == "Montipora grisea" |
           Taxa == "Porites rus" | Taxa == "Turbinaria ornata" |
           Taxa == "Turf") %>%
  ggplot(aes(x = NN_umolL, y = pCover)) +
  geom_point() +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  facet_wrap(~Taxa, scales = "free") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank()) +
  labs(x = "CV Nitrate+Nitrite (%)",
       y = "% Cover")

taxa.nn.lm <- taxa_data_sum %>%
  filter(Taxa == "Padina boryana" | # Taxa == "Montipora grisea" |
           Taxa == "Porites rus" | Taxa == "Turbinaria ornata" |
           Taxa == "Turf") %>%
  ggplot(aes(x = NN_umolL, y = pCover)) +
  geom_point() +
  geom_smooth(method = "lm", formula = "y~x", color = "black") +
  facet_wrap(~Taxa, scales = "free") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank()) +
  labs(x = "CV Nitrate+Nitrite (%)",
       y = "% Cover")

taxa.nn.noLM <- taxa_data_sum %>%
  filter(Taxa == "Padina boryana" | # Taxa == "Montipora grisea" |
           Taxa == "Porites rus" | Taxa == "Turbinaria ornata" |
           Taxa == "Turf") %>%
  ggplot(aes(x = NN_umolL, y = pCover)) +
  geom_point() +
  # geom_smooth(method = "lm", formula = "y~x", color = "black") +
  facet_wrap(~Taxa, scales = "free") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank()) +
  labs(x = "CV Nitrate+Nitrite (%)",
       y = "% Cover")

# ggsave(here("Output", "Supp_Fig6_taxa.nn.poly.png"), taxa.nn.poly, device = "png", width = 5, height = 6)
# ggsave(here("Output", "Supp_Fig6_taxa.nn.lm.png"), taxa.nn.lm, device = "png", width = 5, height = 6)
# ggsave(here("Output", "Supp_Fig6_taxa.nn.noLM.png"), taxa.nn.noLM, device = "png", width = 5, height = 6)
