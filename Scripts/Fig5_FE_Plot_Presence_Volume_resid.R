### Comparing functional presence and volume across each quadrat

### Created by Danielle Barnas
### Created on December 14, 2022
### Modified March 5, 2023

###############################
# LOAD LIBRARIES
###############################
library(tidyverse)
library(here)
library(FD)
library(tripack) # Triangulation of Irregularly Spaced Data
library(geometry) # Mesh Generation and Surface Tessellation
library(matrixStats) # Functions that Apply to Rows and Columns of Matrices (and to Vectors)
library(patchwork)
library(PNWColors)
library(ggrepel)
library(pairwiseAdonis)

###############################
# READ IN DATA
###############################
alphatag <- read_csv(here("Data","CowTag_to_AlphaTag.csv")) %>% mutate(AlphaTag = if_else(CowTagID == "VSEEP", "A-Seep",AlphaTag))
traits <- read_csv(here("Data", "Surveys","Distinct_Taxa.csv"))
comp <- read_csv(here("Data", "Surveys", "Species_Composition_2022.csv"))
meta <- read_csv(here("Data", "Full_Metadata.csv"))
chem <- read_csv(here("Data","Biogeochem", "Nutrients_Processed_All.csv"))
# richness, % richness of community pool, and % volume of community pool
Fric <- read_csv(here("Data", "Sp_FE_Vol.csv"))
resFric <- read_csv(here("Data", "Sp_FE_Vol_res.csv"))
# PCoA axes of traits
fd.coord.sgd <- read.csv(here("Data","FE_4D_coord.csv"), row.names = 1) # class data.frame
# FE with trait groups
fes_traits.sgd <- read_csv(here("Data", "Distinct_FE.csv"))
# species abundances (%) wide format
myspecies <- read_csv(here("Data", "Species_Abundances_wide.csv"))
# species and functional entities
species_entities <- read_csv(here("Data", "Species_FE.csv"))




###############################
# CLEANING AND ANALYSIS
###############################
meta <- meta %>% # fill in seep rugosity value
  mutate(meanRugosity = if_else(CowTagID == "VSEEP", 0.97, meanRugosity))

chem <- chem %>%
  filter(Location == "Varari",
         CowTagID != "V13")
redchem <- chem %>%
  select(CowTagID, Phosphate_umolL, NN_umolL)

comp <- comp %>%
  filter(Location == "Varari",
         CowTagID != "V13")

myspecies <- myspecies %>%
  arrange(CowTagID)

# df with unique functional entities for each row
entity <- fes_traits.sgd


# CowTagIDs as rownames
sgd.sp <- column_to_rownames(.data = myspecies, var = "CowTagID")
sgd.sp <- as.data.frame(sgd.sp)



## Plot convex hull (modified Teixido script)
### Bar plot

# order CowTagID's by distance from the seepage point (matching the map orientation)
tagOrder <- meta %>%
  filter(CowTagID != "V13") %>%
  arrange(dist_to_seep_m) %>% # set arrange factor
  select(CowTagID) %>%
  left_join(alphatag)
# set cowtag order as arrange factor order
tagOrder <- tagOrder$AlphaTag[1:20] # exclude maya's sites


### Functional space using PCoA (modified Teixido script)
### plot convex hull

q <- list()
All.ch.tib <- tibble(x = as.numeric(),
                     y = as.numeric(),
                     CowTagID = as.character())
All.m.sgd <- tibble(PC1 = as.numeric(),
                    PC2 = as.numeric(),
                    PC3 = as.numeric(),
                    PC4 = as.numeric(),
                    CowTagID = as.character(),
                    FE = as.character())

# needs to be df before running for loop
species_entities <- as.data.frame(column_to_rownames(species_entities, var = 'Taxa'))

cowtagOrder <- alphatag %>%
  arrange(AlphaTag) %>%
  select(CowTagID)
cowtagOrder <- cowtagOrder$CowTagID[2:20] # removes seep



for(i in cowtagOrder) {

  tag = i # use for rbinding data below

  species.sgd <- colnames(sgd.sp)[which(sgd.sp[i,] > 0)]
  # only species present in each survey location

  fes_cond.sgd <- species_entities[rownames(species_entities) %in% species.sgd, ]
  # functional entities characterizing those species

  m.sgd <- fd.coord.sgd[rownames(fd.coord.sgd) %in% fes_cond.sgd, ]
  m.sgd <- data.matrix(m.sgd) # parse from data frame to matrix array

  mid.m.sgd <- as_tibble(m.sgd) %>%
    mutate(CowTagID = tag, FE = rownames(m.sgd))
  All.m.sgd <- All.m.sgd %>% rbind(mid.m.sgd)

  tr.sgd <- tri.mesh(m.sgd[,1],m.sgd[,2], duplicate = "remove") # duplicate: default = "error", "strip" = removes all duplicate points, "remove" = leaves one point of duplicate points

  ch.sgd <- convex.hull(tr.sgd)

  ch.tib <- cbind(ch.sgd$x, ch.sgd$y, ch.sgd$i) # parse as tibble df
  colnames(ch.tib) <- c("x", "y", "i")
  ch.tib <- as_tibble(ch.tib) %>%
    select(x,y) %>%
    mutate(CowTagID = tag)

  All.ch.tib <- All.ch.tib %>% rbind(ch.tib)
}

# add VSEEP coordinates
seep.species.sgd <- colnames(sgd.sp)[which(sgd.sp['VSEEP',] > 0)] # select present species from seep
seep.fes_cond.sgd <- species_entities[rownames(species_entities) %in% seep.species.sgd, ]
seep.m.sgd <- fd.coord.sgd[rownames(fd.coord.sgd) %in% seep.fes_cond.sgd, ]
seep.mid.m.sgd <- as_tibble(seep.m.sgd) %>%
  mutate(CowTagID = "VSEEP", FE = rownames(seep.m.sgd))
sub.ch.tib <- seep.mid.m.sgd %>%
  select(CowTagID, x = PC1, y = PC2)


#### COLOR PALETTE
mypalette <- (pnw_palette(name = "Bay", n = 20))
alphapalette <- c(paste0(mypalette[1],"70"), paste0(mypalette[2],"70"), paste0(mypalette[3],"70"),
                  paste0(mypalette[4],"70"), paste0(mypalette[4],"70"), paste0(mypalette[6],"70"),
                  paste0(mypalette[7],"70"), paste0(mypalette[8],"70"), paste0(mypalette[9],"70"),
                  paste0(mypalette[10],"70"), paste0(mypalette[11],"70"), paste0(mypalette[12],"70"),
                  paste0(mypalette[13],"70"), paste0(mypalette[14],"70"), paste0(mypalette[15],"70"),
                  paste0(mypalette[16],"70"), paste0(mypalette[17],"70"), paste0(mypalette[18],"70"),
                  paste0(mypalette[19],"70"), paste0(mypalette[20], "70"))


### Arrange AlphaTag ID's by CV Phosphate
disc.redchem <- redchem %>%
  left_join(resFric) %>%
  select(-c(NbSp:NbFEsP, resSp:resFEp)) %>%
  mutate(Phosphate_umolL = as.character(round(Phosphate_umolL,3))) %>%
  mutate(NN_umolL = as.character(round(NN_umolL,3))) %>%
  arrange(resVol) %>%
  left_join(alphatag)
names(alphapalette) <- (disc.redchem  %>% arrange(Phosphate_umolL) %>% mutate(nAlphaTag = AlphaTag))$nAlphaTag
names(mypalette) <- (disc.redchem %>% arrange(Phosphate_umolL) %>% mutate(nAlphaTag = AlphaTag))$nAlphaTag
myOrder <- (disc.redchem %>% arrange(Phosphate_umolL))$AlphaTag
#volumeOrder <- (disc.redchem %>% arrange(resVol))$AlphaTag
volumeOrder <- (disc.redchem %>% arrange(Vol8D))$AlphaTag


# add VSEEP point
All.ch.tib <- All.ch.tib %>%
  rbind(sub.ch.tib) %>%
  left_join(alphatag) %>%
  select(-CowTagID) %>%
  group_by(AlphaTag) %>%
  mutate(Area = polyarea(x,y)) %>%
  mutate(nAlphaTag = factor(AlphaTag, levels = names(alphapalette))) %>%
  mutate(AlphaTag = factor(AlphaTag, levels = myOrder))
#AreaOrder <- unique(All.ch.tib$AlphaTag)
#All.ch.tib$AlphaTag <- factor(All.ch.tib$AlphaTag, levels = AreaOrder) # arrange facet by polygon area

#PolyArea <- All.ch.tib %>%
#  mutate(AlphaTag, Area)

All.m.sgd <- All.m.sgd %>%
  rbind(seep.mid.m.sgd) %>%
  left_join(alphatag) %>%
  select(-CowTagID) %>%
#  left_join(PolyArea) %>%
  mutate(nAlphaTag = factor(AlphaTag, levels = names(alphapalette)),
         AlphaTag = factor(AlphaTag, levels = myOrder))
#         AlphaTag = factor(AlphaTag, levels = AreaOrder))



# graph faceted polygons showing functional volume
qAll <- All.ch.tib %>%
  ggplot(aes(x = x, y = y)) +
  geom_polygon(aes(fill = nAlphaTag, color = nAlphaTag), alpha = 0.5) + # create polygon using product of convex.hull(tri.mesh)
  labs(x = "PCoA 1", y = "PCoA 2") +
  geom_point(data = as_tibble(All.m.sgd), aes(x = PC1, y = PC2, color = nAlphaTag)) +
  theme_bw() +
  theme(#legend.position = "none",
    panel.grid = element_blank(),
    panel.spacing = unit(1, "lines")) + # increase facet wrap panel spacing
  scale_fill_manual(values = alphapalette) +
  scale_color_manual(values = mypalette) +
  facet_wrap(~AlphaTag) +
  xlim(min(fd.coord.sgd[,1]), max(fd.coord.sgd[,1])) + ylim(min(fd.coord.sgd[,2]), max(fd.coord.sgd[,2])) +
  theme(strip.background = element_rect(fill = "white"))
qAll


### COMBINE VOLUME FIGURE FROM Presence SCRIPT WITH ABOVE TO INCLUDE POLYGON


#load data
ab.sgd <- read_csv(here("Data", "Species_Abundances_wide.csv"))
ab.sgd <- as.data.frame(column_to_rownames(ab.sgd, var = 'CowTagID')) # move tag names to rownames and make data.frame class
spe_fes.sgd <- as.data.frame(read_csv(here("Data", "Species_FE.csv")))
#alphatag <- read_csv(here("Data","CowTag_to_AlphaTag.csv"))

################################## Data manipulation and arrangements

ab.conditions.sgd <- ab.sgd

################################# compute abundance of FEs

fes.sgd <- levels(as_factor(spe_fes.sgd$FE))

ab.conditions.sgd <- rownames_to_column(ab.conditions.sgd, var = "CowTagID")
ab.conditions.sgd2 <- ab.conditions.sgd %>%
  pivot_longer(names_to = "Taxa", values_to = "pCover", cols = 2:ncol(ab.conditions.sgd)) %>%
  left_join(spe_fes.sgd) %>%
  group_by(CowTagID, FE) %>%
  summarise(pCover = sum(pCover)) %>%
  ungroup()
ab.conditions.sgd <- ab.conditions.sgd2 %>%
  pivot_wider(names_from = FE, values_from = pCover)




###########################################################################
## FE VOLUME AND NUTRIENTS
###########################################################################

## FACET ORDERED BY RESIDUAL VOLUME
## COLORED BY PHOSPHATE

# Figure 5A. Overall distribution of FE abundance across the functional space colored by Phosphate_umolL

# order alpha tag by phosphate levels and volume residuals
n.alphatag <- alphatag %>% mutate(nAlphaTag = factor(alphatag$AlphaTag, levels = names(alphapalette))) %>% select(-AlphaTag)
alphatag$AlphaTag <- factor(alphatag$AlphaTag, levels = myOrder)

## relative abundance in ggplot
fig5.fd.sgd <- rownames_to_column(as.data.frame(fd.coord.sgd), var = "FE") %>%
  full_join(ab.conditions.sgd2) %>%
  left_join(alphatag) %>%
  left_join(n.alphatag) %>%
  arrange(AlphaTag)

fig5adist <- fig5.fd.sgd %>%
  filter(pCover > 0) %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(size = pCover,
                 color = nAlphaTag,
                 fill = nAlphaTag),
             shape = 21,
             show.legend = FALSE) + # shape of a filable circle. lets us fill with alpha values
  geom_polygon(data = All.ch.tib %>% left_join(n.alphatag),
               aes(x = x, y = y,
                   color = nAlphaTag),
               alpha = 0.5,
               fill = NA) + # no fill on the polygon
  labs(x = "PCoA1", y = "PCoA2", color = "Site") +
  facet_wrap(~AlphaTag) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 9)) +
  scale_fill_manual(values = alphapalette) +
  scale_color_manual(values = mypalette)

fig5adist

#ggsave(here("Output", "PaperFigures", "Fig5a_Vol_Abund_PCoA_PO4.png"), fig5adist, height = 6, width = 7)



#####################################
# FIGURE 5B FE POINTS PCOA
#####################################

fd.coord.sgd.fe <- read.csv(here("Data","FE_4D_coord.csv"), row.names = 1) # class data.frame

### View location of each functional entity:
fd.coord.sgd.tibble <- as_tibble(rownames_to_column(as.data.frame(fd.coord.sgd.fe))) %>%
  rename(FE = "rowname")


## View functional trait faceted figures
fe_group_pcoa <- fd.coord.sgd.tibble %>%
  separate(FE, into = c('Phyla','Morphology','Calcification','Trophic Group'),
           sep = ",", remove = F) %>%
  pivot_longer(cols = 'Phyla':'Trophic Group', names_to = "Group", values_to = "Trait")
fe_group_pcoa$Group <- factor(fe_group_pcoa$Group,
                              levels = c("Phyla", "Morphology", "Calcification", "Trophic Group"))
plot_fe_group_pcoa <- fe_group_pcoa %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(shape = 21, fill = "darkgrey") +
  geom_text_repel(aes(label = Trait),
                  size = 3,
                  max.overlaps = 18) +
  labs(x = "PCoA1", y = "PCoA2") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 9)) +
  facet_wrap(~Group)
plot_fe_group_pcoa

ggsave(here("Output", "PaperFigures", "Fig5b_FE_grouped_pcoa.png"), plot_fe_group_pcoa, width = 6, height = 6)


#### PATCH PLOTS TOGETHER FOR FIGURE 5

Figure5 <- fig5adist / plot_fe_group_pcoa +
  plot_annotation(tag_levels = "A")
Figure5

ggsave(here("Output", "PaperFigures", "Fig5_FEV_pcoa.png"),Figure5, device = "png", width = 6, height = 10)


# just to get the Phosphate scale bar
apal<-(pnw_palette("Bay"))
po4_scale <- redchem %>%
  ggplot(aes(x = NN_umolL, y = Phosphate_umolL, color = Phosphate_umolL)) +
  geom_point() +
  scale_colour_gradientn(colours = apal) +
  labs(color = expression("CV Phosphate (%)")) +
  theme(legend.key.size = unit(2, 'cm'),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20))
ggsave(here("Output", "PaperFigures", "Fig5_Scale_bar.png"),po4_scale, device = "png", width = 6, height = 6)


#########################################################################################################
#########################################################################################################


#####################
# Blank figure with points and outline
#####################
blankPoly <- fig5.fd.sgd %>%
  filter(pCover > 0) %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(shape = 21, size = 4, fill = "grey") +
  geom_polygon(data = All.ch.tib %>% filter(AlphaTag == "K"), # K, N, and S all have the same largest area
               aes(x = x, y = y, color = AlphaTag),
               alpha = 0.5,
               fill = NA) + # no fill on the polygon
  labs(x = "PCoA1", y = "PCoA2") +
  theme_bw() +
  xlim(-0.35, 0.65) +
  ylim(-0.5, 0.5) +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  scale_color_manual(values = "black")


blankPoly_small <- fig5.fd.sgd %>%
  filter(pCover > 0) %>%
  filter(AlphaTag == "D",
         AlphaTag != "A-Seep") %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(shape = 21, size = 4, fill = "grey") +
  geom_polygon(data = All.ch.tib %>% filter(AlphaTag == "D"), # D has the smallest area after Seep
               aes(x = x, y = y, color = AlphaTag),
               alpha = 0.5,
               fill = NA) + # no fill on the polygon
  labs(x = "PCoA1", y = "PCoA2") +
  xlim(-0.3, 0.6) +
  ylim(-0.5, 0.5) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  scale_color_manual(values = "black")


#####################
# VARIANCE IN HIGH, MODERATE, LOW
#####################

# to look at variation of FE dispersion, look at a distance matrix
Alpha_FE <- ab.conditions.sgd2 %>%
  left_join(alphatag) %>%
  filter(pCover > 0) %>%
  select(AlphaTag, FE) %>%
  distinct()
Alpha_chem <- redchem %>%
  left_join(alphatag) %>%
  select(AlphaTag, NN_umolL,Phosphate_umolL)

#dist(fd.coord.sgd %>% select(PC1, PC2), method = "euclidean")
dist_FE <- as.data.frame(as.matrix(dist(fd.coord.sgd, method = "euclidean"))) # distance across 4D
dist_chem <- rownames_to_column(dist_FE, var = "FE") %>%
  right_join(Alpha_FE) %>%
  filter(AlphaTag != "A-Seep") %>%
  pivot_longer(cols = 2:(ncol(.)-1), names_to = "crossFE", values_to = "distances") %>%
  group_by(AlphaTag) %>%
  summarise(Avg_dist = mean(distances),
            sd_dist = sd(distances),
            se_dist = plotrix::std.error(distances)) %>%
  arrange(desc(Avg_dist)) %>%
  left_join(Alpha_chem) %>%
  ungroup() %>%
  mutate(comm_mean = mean(Avg_dist),
         var_from_mean = Avg_dist - comm_mean)

summary(lm(data = dist_chem, var_from_mean ~ poly(NN_umolL,2)))

ggplot(data = dist_chem,
       aes(x = NN_umolL, y = Avg_dist)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_text_repel(aes(label = AlphaTag))



###########################################################################
## RESIDUALS AND NUTRIENTS
###########################################################################


# prepare chem data for joining other df for graphing
redchem <- redchem %>%
  left_join(alphatag) %>%
  arrange(alphatag)


summary(lm(data = resFric %>% left_join(redchem) %>% filter(CowTagID != "VSEEP"),
           resSpp ~ poly(Phosphate_umolL, 2)))
summary(lm(data = resFric %>% left_join(redchem) %>% filter(CowTagID != "VSEEP"),
           resFEp ~ poly(Phosphate_umolL, 2)))
summary(lm(data = resFric %>% left_join(redchem) %>% filter(CowTagID != "VSEEP"),
           resVol ~ poly(Phosphate_umolL, 2)))

# color palette assignment
cols <- pnw_palette("Bay",20,type="continuous")
#cols <- rev(cols) # reverse color pattern so high sgd gets red
tagorder <- n.alphatag %>% arrange(nAlphaTag)
names(cols) <- tagorder$nAlphaTag # name colors by Phosphate gradient

# Supplemental Figure. Species and functional diversity changes along SGD gradient
# All volumes in distinct plots
# RESIDUALS
pphos <- resFric %>%
  left_join(redchem) %>%
  left_join(n.alphatag) %>%
  select('% SR' = resSpp, '% FER' = resFEp, '% FEV' = resVol, AlphaTag, nAlphaTag) %>%
  pivot_longer(cols = 1:3, names_to = "Parameters", values_to = "Values") %>%
  mutate(Parameters = factor(Parameters, levels = c("% SR", "% FER", "% FEV"))) %>%
  # plot facet_wrapped
  ggplot(aes(x = Parameters, y = Values, fill = nAlphaTag)) +
  geom_col(color = "black") +
  facet_wrap(~nAlphaTag) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90),
        legend.position = "none") +
  #ylim(0,100) +
  scale_fill_manual(values = cols) +
  geom_hline(yintercept = 0) +
  labs(fill = "Survey Location", x = "", y = "Residual proportions (%)") +
  theme(strip.background = element_rect(fill = "white"))
pphos

# RAW
pphosRaw <- resFric %>%
  left_join(redchem) %>%
  left_join(n.alphatag) %>%
  select('% SR' = NbSpP, '% FER' = NbFEsP, '% FEV' = Vol8D, AlphaTag, nAlphaTag) %>%
  pivot_longer(cols = 1:3, names_to = "Parameters", values_to = "Values") %>%
  mutate(Parameters = factor(Parameters, levels = c("% SR", "% FER", "% FEV"))) %>%
  # plot facet_wrapped
  ggplot(aes(x = Parameters, y = Values, fill = nAlphaTag)) +
  geom_col(color = "black") +
  facet_wrap(~nAlphaTag) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90),
        legend.position = "none") +
  ylim(0,100) +
  scale_fill_manual(values = cols) +
  labs(fill = "Survey Location", x = "", y = "Community Proportions (%)") +
  geom_text(aes(x = Parameters, label = round(Values,0)),
           size = 3, vjust = -0.4) +
  theme(strip.background = element_rect(fill = "white"))
pphosRaw

pphosPlots <- pphosRaw / pphos


