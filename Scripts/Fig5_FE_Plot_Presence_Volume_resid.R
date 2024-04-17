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
fd.coord.sgd <- read.csv(here("Data","FE_4D_coord_dmb.csv"), row.names = 1) # class data.frame
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
  filter(Season == "Dry") %>%
  filter(Location == "Varari",
         CowTagID != "V13") %>%
  select(CowTagID, Parameters, CVSeasonal) %>%
  pivot_wider(names_from = Parameters, values_from = CVSeasonal)
redchem <- chem %>% select(CowTagID, Phosphate_umolL, NN_umolL)

comp <- comp %>%
  filter(Location == "Varari") %>%  # only analyze varari for now
  filter(CowTagID != "V13")

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
  filter(Location == "Varari",
         CowTagID != "V13") %>%
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
        strip.background = element_rect(fill = "white")) +
  scale_fill_manual(values = alphapalette) +
  scale_color_manual(values = mypalette)

fig5adist

ggsave(here("Output", "PaperFigures", "Fig5a_Vol_Abund_PCoA_PO4.png"), fig2bdist, height = 6, width = 7)



#####################################
# FIGURE 5B FE POINTS PCOA
#####################################

fd.coord.sgd.fe <- read.csv(here("Data","FE_4D_coord_dmb.csv"), row.names = 1) # class data.frame

### View location of each functional entity:
fd.coord.sgd.tibble <- as_tibble(rownames_to_column(as.data.frame(fd.coord.sgd.fe))) %>%
  rename(FE = "rowname")


## View functional trait faceted figures
fe_group_pcoa <- fd.coord.sgd.tibble %>%
  separate(FE, into = c('Phyla','Morphology','Calcification','Energetic Resource'),
           sep = ",", remove = F) %>%
  pivot_longer(cols = 'Phyla':'Energetic Resource', names_to = "Group", values_to = "Trait")
fe_group_pcoa$Group <- factor(fe_group_pcoa$Group,
                              levels = c("Phyla", "Morphology", "Calcification", "Energetic Resource"))
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
        strip.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  facet_wrap(~Group)
plot_fe_group_pcoa

ggsave(here("Output", "PaperFigures", "Fig5b_FE_grouped_pcoa.png"), plot_fe_group_pcoa, width = 6, height = 6)


#### PATCH PLOTS TOGETHER FOR FIGURE 5

Figure5 <- fig2bdist / plot_fe_group_pcoa +
  plot_annotation(tag_levels = "A")

ggsave(here("Output", "PaperFigures", "Fig5_FEV_pcoa.png"),Figure5, device = "png", width = 6, height = 10)


# just to get the Phosphate scale bar
apal<-(pnw_palette("Bay"))
redchem %>%
  ggplot(aes(x = NN_umolL, y = Phosphate_umolL, color = Phosphate_umolL)) +
  geom_point() +
  scale_colour_gradientn(colours = apal) +
  labs(color = "Phosphate (umol/L)") +
  theme(legend.key.size = unit(2, 'cm'),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20))


#########################################################################################################
#########################################################################################################


#####################
# Blank figure with points and outline
#####################
blankPoly <- fig2.fd.sgd %>%
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

ggsave(here("Output","Blank_polygon.png"), blankPoly, height = 6, width = 7)



blankPoly_small <- fig2.fd.sgd %>%
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

ggsave(here("Output", "Blank_polygon_small.png"), blankPoly_small, height = 6, width = 7)


#####################
# VARIANCE IN HIGH, MODERATE, LOW
#####################
library(vegan) # nMDS and permanova
library(pairwiseAdonis)

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

## check this using ANOSIM?
## look up other literature on this


#### PREP DATA FOR MODEL
VarResFric <- resFric %>%
  left_join(alphatag) %>%
  select(AlphaTag, Vol8D)
low <- VarResFric %>%
  filter(AlphaTag == "B" | AlphaTag == "F" | AlphaTag == "E" | AlphaTag == "Q") %>%
  mutate(Varcat = "Low")
high <- VarResFric %>%
  filter(AlphaTag == "D" | AlphaTag == "C" | AlphaTag == "N" | AlphaTag == "H" | AlphaTag == "A-Seep") %>%
  mutate(Varcat = "High")
LH <- rbind(low,high)

### try with abundances
LH <- LH %>% select(-Vol8D)
permData <- rownames_to_column(as.data.frame(fd.coord.sgd), var = "FE") %>%
  full_join(ab.conditions.sgd2) %>%
  select(FE, CowTagID, pCover) %>%
  left_join(alphatag) %>%
  filter(AlphaTag != "A-Seep") %>%
  select(-CowTagID) %>% # using alpha tag as the identifier
  pivot_wider(names_from = FE, values_from = pCover) %>%
  left_join(LH) %>%
  mutate(Varcat = if_else(is.na(Varcat), "Moderate", Varcat)) %>%
  relocate(Varcat, .after = AlphaTag) %>%
  select(-AlphaTag)


permanovamodel<-adonis2(permData[,-1]~Varcat, permData, permutations = 999,method="bray")
permanovamodel

#assume that the dispersion among data is the same in each group. We can test with assumption with a PermDisp test:
disper<-vegdist(permData[,-1])
betadisper(disper, permData$Varcat)
#Look at the Average distance to median...these numbers should be reasonably similar
#A rule of thumb is that one number should not be twice as high as any other

#An option for doing post-hoc pairwise comparisons in R
library(pairwiseAdonis)
pairwise.adonis(permData[-1], permData$Varcat, perm=999)






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
# geom_text(aes(x = Parameters, label = round(Values,0)),
#           size = 3, vjust = -0.4)
pphos

ggsave(here("Output", "xSupp_barplot_res_phosphate.png"), pphos, width = 6, height = 5)

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

ggsave(here("Output", "xSupp_barplot_raw_phosphate.png"), pphosRaw, width = 6, height = 5)

pphosPlots <- pphosRaw / pphos +
  plot_annotation(tag_levels = "A")
ggsave(here("Output", "xSuppFig2_barplot_phosphate.png"), pphosPlots, width = 5, height = 8)











###########################################################################
## SPECIES VOLUME AND NUTRIENTS
###########################################################################



### Functional space using PCoA (modified Teixido script)

q <- list()
All.ch.tib <- tibble(x = as.numeric(),
                     y = as.numeric(),
                     CowTagID = as.character())
All.m.sgd <- tibble(PC1 = as.numeric(),
                    PC2 = as.numeric(),
                    PC3 = as.numeric(),
                    PC4 = as.numeric(),
                    CowTagID = as.character(),
                    Sp = as.character())


cowtagOrder <- alphatag %>%
  arrange(AlphaTag) %>%
  select(CowTagID) %>%
  filter(CowTagID != "VSEEP", CowTagID != "V14")
cowtagOrder <- cowtagOrder$CowTagID[1:18] # removes seep

fd.coord.sgd <- read_csv(here("Data", "FE_4D_coord_species.csv"))
fd.coord.sgd <- column_to_rownames(fd.coord.sgd, var = '...1')

for(i in cowtagOrder) {

  tag = i # use for rbinding data below

  species.sgd <- colnames(sgd.sp)[which(sgd.sp[i,] > 0)]
  # only species present in each treatment

  fes_cond.sgd <- species.sgd #species_entities[rownames(species_entities) %in% species.sgd, ]

  m.sgd <- fd.coord.sgd[rownames(fd.coord.sgd) %in% fes_cond.sgd, ]
  m.sgd <- data.matrix(m.sgd) # parse from data frame to matrix array

  mid.m.sgd <- as_tibble(m.sgd) %>%
    mutate(CowTagID = tag, Species = rownames(m.sgd))
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
seep.fes_cond.sgd <- seep.species.sgd
seep.m.sgd <- fd.coord.sgd[rownames(fd.coord.sgd) %in% seep.fes_cond.sgd, ]
seep.mid.m.sgd <- as_tibble(seep.m.sgd) %>%
  mutate(CowTagID = "VSEEP", Species = rownames(seep.m.sgd))
seep.sub.ch.tib <- seep.mid.m.sgd %>%
  select(CowTagID, x = PC1, y = PC2)

# add V14 coordinates
v14.species.sgd <- colnames(sgd.sp)[which(sgd.sp['V14',] > 0)] # select present species from seep
v14.fes_cond.sgd <- v14.species.sgd
v14.m.sgd <- fd.coord.sgd[rownames(fd.coord.sgd) %in% v14.fes_cond.sgd, ]
v14.mid.m.sgd <- as_tibble(v14.m.sgd) %>%
  mutate(CowTagID = "V14", Species = rownames(v14.m.sgd))
v14.sub.ch.tib <- v14.mid.m.sgd %>%
  select(CowTagID, x = PC1, y = PC2)
nopoly<- v14.sub.ch.tib[c(3, 5, 7,8),] # removes extra points in the polygon muddling the volume graph
v14.sub.ch.tib <- v14.sub.ch.tib %>%
  anti_join(nopoly)

### Arrange AlphaTag ID's by CV Phosphate
names(alphapalette) <- (disc.redchem %>% left_join(alphatag) %>% arrange(Phosphate_umolL) %>% mutate(nAlphaTag = AlphaTag))$nAlphaTag
names(mypalette) <- (disc.redchem %>% left_join(alphatag) %>% arrange(Phosphate_umolL) %>% mutate(nAlphaTag = AlphaTag))$nAlphaTag
myOrder <- (disc.redchem %>% left_join(alphatag) %>% arrange(Phosphate_umolL))$AlphaTag

# add VSEEP and V14 points

All.ch.tib <- All.ch.tib %>%
  rbind(seep.sub.ch.tib) %>%
  rbind(v14.sub.ch.tib) %>%
  left_join(alphatag) %>%
  left_join(redchem) %>%
  select(-CowTagID) %>%
  mutate(nAlphaTag = factor(AlphaTag, levels = names(alphapalette)),
         AlphaTag = factor(AlphaTag, levels = myOrder))
All.m.sgd <- All.m.sgd %>%
  rbind(seep.mid.m.sgd) %>%
  rbind(v14.mid.m.sgd) %>%
  left_join(alphatag) %>%
  left_join(redchem) %>%
  select(-CowTagID) %>%
  mutate(nAlphaTag = factor(AlphaTag, levels = names(alphapalette)),
         AlphaTag = factor(AlphaTag, levels = myOrder))

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
  facet_wrap(~nAlphaTag) +
  xlim(min(fd.coord.sgd[,1]), max(fd.coord.sgd[,1])) + ylim(min(fd.coord.sgd[,2]), max(fd.coord.sgd[,2])) +
  theme(strip.background = element_rect(fill = "white"))
qAll


#ggsave(here("Output", "PaperFigures", "Teixido_Figure1volume_dmb_CowTags.png"), qAll, width = 8, height = 5)



### COMBINE VOLUME FIGURE FROM Presence SCRIPT WITH ABOVE TO INCLUDE POLYGON


#load data
ab.sgd <- read_csv(here("Data", "Species_Abundances_wide.csv"))
ab.sgd <- as.data.frame(column_to_rownames(ab.sgd, var = 'CowTagID')) # move tag names to rownames and make data.frame class


################################## Data manipulation and arrangements

ab.conditions.sgd <- ab.sgd

################################# compute abundance of FEs for the three conditions

fes.sgd <- levels(as_factor(colnames(sgd.sp)))

ab.conditions.sgd <- rownames_to_column(ab.conditions.sgd, var = "CowTagID")
ab.conditions.sgd2 <- ab.conditions.sgd %>%
  pivot_longer(names_to = "Taxa", values_to = "pCover", cols = 2:ncol(ab.conditions.sgd)) %>%
  #left_join(spe_fes.sgd) %>%
  group_by(CowTagID, Taxa) %>%
  summarise(pCover = sum(pCover)) %>%
  ungroup()
ab.conditions.sgd <- ab.conditions.sgd2 %>%
  pivot_wider(names_from = Taxa, values_from = pCover)




######################

# Figure 2. Overall distribution of Species abundance across the species composition space


## relative abundance in ggplot
fig2.fd.sgd <- rownames_to_column(as.data.frame(fd.coord.sgd), var = "Taxa") %>%
  full_join(ab.conditions.sgd2) %>%
  left_join(alphatag)

fig4B_species <- fig2.fd.sgd %>%
  filter(pCover > 0) %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(size = pCover,
                 color = AlphaTag,
                 fill = AlphaTag),
             shape = 21) + # shape of a fillable circle. lets us fill with alpha values
  geom_polygon(data = All.ch.tib,
               aes(x = x, y = y,
                   color = AlphaTag),
               alpha = 0.5,
               fill = NA) + # no fill on the polygon
  labs(x = "PCoA1", y = "PCoA2") +
  facet_wrap(~AlphaTag) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        strip.background = element_rect(fill = "white")) +
  scale_fill_manual(values = alphapalette) +
  scale_color_manual(values = mypalette)

fig4B_species

ggsave(here("Output", "xFig_species_Vol_Abund_PCoA_PO4.png"), fig4B_species, height = 5, width = 7)

##########################################################################################
##########################################################################################




###########################################################################################################################
###########################################################################################################################
### Biplot showing one point per cowtag id colored by nutrient values
###########################################################################################################################

ct.coord <- read_csv(here("Data","FE_4D_coord_cowtags.csv"))


## relative abundance in ggplot
fig2.fd.sgd <- ct.coord %>%
  rename(CowTagID = '...1') %>%
  filter(CowTagID != "VSEEP") %>%
  #full_join(ab.conditions.sgd2) %>%
  left_join(alphatag) %>%
  arrange(AlphaTag)

fig2b_cowtag_biplot <- fig2.fd.sgd %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = AlphaTag),
             shape = 21,
             color = "black",
             size = 4) + # shape of a fillable circle. lets us fill with alpha values
  labs(x = "PCoA1", y = "PCoA2",
       fill = "Phosphate (umol/L") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_fill_manual(values = alphapalette)
#scale_color_manual(values = mypalette)

fig2b_cowtag_biplot

ggsave(here("Output", "CowTag_biplot_pcoa.png"), fig2b_cowtag_biplot, width = 6, height = 6)

###########################################################################################################################
### Biplot showing one point per species id colored by nutrient values
###########################################################################################################################

sp.coord <- read_csv(here("Data","FE_4D_coord_species.csv"))

## relative abundance in ggplot
fig2.fd.sgd <- sp.coord %>%
  rename(Species = '...1')

fig2b_species_biplot <- fig2.fd.sgd %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(shape = 21,
             fill = paste0(mypalette[14],"70"),
             color = "black",
             size = 4) + # shape of a fillable circle. lets us fill with alpha values
  labs(x = "PCoA1", y = "PCoA2") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_text_repel(aes(label = Species), size = 2.5)

fig2b_species_biplot

ggsave(here("Output", "Species_biplot_pcoa.png"), fig2b_species_biplot, width = 6, height = 6)








## Can use the three values above (SpR, FER, Vol4D), and also community composition: either relative abundance or presence-absence
## then can do a permanova / nMDS of community comp with the volume / FErichness



### relative abundance
FE_nmds_data <- myspecies %>%
  filter(CowTagID != "VSEEP") %>%  # remove seep for nMDS for now
  pivot_longer(cols = Turf:'Caulerpa racemosa', names_to = "Taxa", values_to = "pCover") %>%
  filter(pCover > 0) %>%
  left_join(as_tibble(rownames_to_column(species_entities, var = "Taxa"))) %>%
  group_by(CowTagID, FE) %>% # get relative abundance of FE (pCvoer is already percent, so just add percentages of FE)
  mutate(pCoverFE = sum(pCover)) %>%
  distinct(CowTagID, FE, pCoverFE) %>%
  drop_na(FE) %>%
  pivot_wider(names_from = FE, values_from = pCoverFE) %>% # longform for the nmds and will establish absence through NAs
  mutate_at(vars(2:ncol(.)), .funs = ~if_else(is.na(.), 0, .)) # zero for NA's 1's for presence
# will cbind cowtags later
# set levels as numerical order of plates
CTlevels <- c('V1','V2','V3','V4','V5','V6','V7','V8','V9','V10','V11','V12','V14','V15','V16','V17','V18','V19','V20')
FE_nmds_data$CowTagID <- factor(FE_nmds_data$CowTagID, levels = CTlevels)
# arrange by cowtag and then remove for nmds
FE_nmds_data <- FE_nmds_data %>%
  arrange(CowTagID) %>%
  ungroup() %>%
  select(-CowTagID)


ord1 <- metaMDS(FE_nmds_data, k=2, distance='jaccard') # jaccard for P/A

# stress with k=2 dimensions. Is it < 0.3?
ord1$stress

# stress plot - want to minimize scatter
stressplot(ord1)

#param_mds <- nMDS_species(ord1) # MDS1 and MDS2 for FEs
# get points for species
Group <- rownames(ord1$species) # get characteristic names
MDS1 <- c(ord1$species[,1]) # MDS1 for characteristics
MDS2 <- c(ord1$species[,2]) # MDS2 for characteristics
Data <- as_tibble(cbind(Group, MDS1, MDS2)) %>%  # bind all cols into tibble
  mutate(MDS1 = as.numeric(MDS1), # as numeric
         MDS2 = as.numeric(MDS2)) %>%
  #mutate(Taxon_Group = if_else(Taxa == "Hard Substrate", "Abiotic", Taxon_Group)) %>%
  select(MDS1, MDS2, Group)


#param_mds_cat <- nMDS_points(ord1, meta, c('CowTagID', 'dist_to_seep_m')) # MDS1 and MDS2 for CowTagID
Groupb <- as.character(CTlevels) # assign CowTagID
MDS1b <- ord1$points[,1] # MDS1 for CowTagID
MDS2b <- ord1$points[,2] # MDS2 for CowTagID
Datab <- as_tibble(cbind(Groupb, MDS1b, MDS2b)) %>%  # bind all cols into tibble
  mutate(MDS1b = as.numeric(MDS1b), # as numeric
         MDS2b = as.numeric(MDS2b)) %>%
  rename(CowTagID = Groupb)

joinDF <- chem %>%
  select(CowTagID, Phosphate_umolL)

Datab <- Datab %>%
  left_join(joinDF)


## plot
nMDSplot <- ggplot(data = Data,
                   aes(x = MDS1,
                       y = MDS2)) +
  geom_point(color = "black") +
  geom_point(data = Datab,
             aes(x = MDS1b,
                 y = MDS2b,
                 color = (Phosphate_umolL)),
             size = 3) +
  geom_text_repel(data = Data, # site characteristics
                  aes(x = MDS1,
                      y = MDS2,
                      label = Group),
                  size = 3) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "right") +
  geom_label_repel(data = Datab, # site characteristics
                   aes(x = MDS1b,
                       y = MDS2b,
                       label = CowTagID),
                   size = 5,
                   max.overlaps = 16) + # increase from 10 because too many FEs overlapping
  scale_color_gradient(low = "red", high = "yellow")
nMDSplot

ggsave(here("Output", "FE_nmds_plot.png"), nMDSplot, width = 10, height = 10)

### PERMANOVA
richPermFull <- cbind(Groupb, FE_nmds_data) %>%  # bind cowTagIDs
  rename(CowTagID = Groupb) %>%
  left_join(joinDF) %>%
  mutate(rel = if_else(Phosphate_umolL > 0.15, "High", # D, C, N, H
                       if_else(Phosphate_umolL < 0.06, "Low", # Q, E, F, B
                               "Moderate")))



permanovamodel<-adonis2(richPermFull[,2:24]~rel, richPermFull, permutations = 999,
                        method="bray") # should change out cowtagid with some grouping name
permanovamodel

# #If we are to trust the results of the permanova, then we have to assume that the dispersion among
# #data is the same in each group. We can test with assumption with a PermDisp test:
# disper<-vegdist(richPermFull[,2:25])
# betadisper(disper, richPermFull$relDist)
# #Look at the Average distance to median...these numbers should be reasonably similar
# #A rule of thumb is that one number should not be twice as high as any other
#
# pairwise.adonis(richPermFull[2:25], richPermFull$relDist, perm=999)
#
# #Get coefficients to see which species are most important in explaining site differences:
# #permanovamodel$coefficients
#
#
#
