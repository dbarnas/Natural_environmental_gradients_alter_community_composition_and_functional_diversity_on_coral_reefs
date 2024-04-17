### Comparing functional presence and volume across each quadrat

### This script creates and exports a dataframe of Species Richness,
### Functional Entity Richness, and Functional Volume for each CowTagID

## Calculate convex hull (modified Teixido script)
#### Calculate convex hull

### Created by Danielle Barnas
### Created November 2022
### Modified February 26, 2023

##### LOAD LIBRARIES #####

library(tidyverse)
library(here)
library(FD)
library(tripack) # Triangulation of Irregularly Spaced Data
library(geometry) # Mesh Generation and Surface Tessellation
library(matrixStats) # Functions that Apply to Rows and Columns of Matrices (and to Vectors)
library(patchwork)
library(PNWColors)


##### READ IN DATA #####

traits <- read_csv(here("Data", "Surveys","Distinct_Taxa.csv"))
myspecies <- read_csv(here("Data", "Surveys", "Species_Composition_2022.csv"))
meta <- read_csv(here("Data", "Full_Metadata.csv"))
chem <- read_csv(here("Data","Biogeochem", "Nutrients_Processed_All.csv")) %>%
  filter(Season == "Dry") %>%
  filter(Location == "Varari",
         #CowTagID != "VSEEP" &
           CowTagID != "V13") %>%
  select(CowTagID, Parameters, CVSeasonal) %>%
  pivot_wider(names_from = Parameters, values_from = CVSeasonal)

myspecies <- myspecies %>%
  filter(Location == "Varari", # only analyze varari for now
         #CowTagID != "VSEEP",
         CowTagID != "V13")

# only cowtag ID's
quad.label <- myspecies %>%
  select(CowTagID) %>%
  distinct()



# create attribute for Functional Entity of each species
traits <- traits %>%
  full_join(myspecies) %>% # apply traits to species composition df
  filter(Location == "Varari",
         #CowTagID != "VSEEP",
         CowTagID != "V13",
         Identified == "yes",
         Taxon_Group != "Hard Substrate" &
           Taxon_Group != "Sand") %>%
  select(Taxa,
         Taxon_Group,
         Morph2,
         #Symb,
         Calc,
         #Life_Span,
         #MS_cat,
         #GR_cat,
         ER
         #FM
         ) %>%
  unite(col = "FE", Taxon_Group:ER, sep = ",", remove = F) %>%
  relocate(FE, .after = ER) %>%
  distinct()

# write csv for all distinct FE combinations
fes_traits.sgd <- traits %>%
  select(FE = FE, Taxon_Group:ER) %>%
  distinct()

write_csv(fes_traits.sgd, here("Data", "Distinct_FE.csv"))




# species abundance as wide format species
myspecies <- myspecies %>%
  filter(Taxa != "Bare Rock" & # only include biological data
           Taxa != "Sand" &
           Taxa != "Rubble") %>%
  select(CowTagID, Taxa, SpeciesCounts) %>%
  group_by(CowTagID) %>%
  mutate(totalCount = sum(SpeciesCounts), # calculate percent cover of taxa
         cover = SpeciesCounts / totalCount * 100) %>%
  ungroup() %>%
  select(-c(SpeciesCounts,totalCount)) %>%
  group_by(CowTagID, Taxa) %>% # make sure all distinct taxa are summed together
  mutate(cover = sum(cover)) %>%
  distinct() %>%  # include only one of each taxa category. may have some repeats after ID'ing previously unknown organisms
  pivot_wider(names_from = Taxa, values_from = cover) %>%
  mutate_all(.funs = ~if_else(is.na(.), 0, .)) %>%  # zero percent for any NA values
  ungroup()

# write csv for species abundances
write_csv(myspecies, here("Data", "Species_Abundances_wide.csv"))



# select only for species and functional entities
species_entities <- traits %>%
  select(Taxa, FE)

# write csv for species entities tied to FE
write_csv(species_entities, here("Data", "Species_FE.csv"))




# df with unique functional entities for each row
entity <- traits %>%
  select(-Taxa) %>%
  distinct() %>%
  relocate(FE, .before = Taxon_Group) %>%
  mutate_all(.funs = as_factor) # groups need to be factors to run quality_funct_space()
entity <- column_to_rownames(entity, var = "FE")

#  CowTagIDs as factors, to be used as relative SGD identifiers along gradient
relative.sgd <- quad.label$CowTagID

## Computing multidimensional functional space

# source function from Teixido: function for computing the quality of functional dendrogramm and multidimensional functional spaces
source(here("Scripts", "Teixido", "quality_funct_space.R"))

qfs <- quality_funct_space(mat_funct = entity, # distinct functional trait entities (NAs not allowed. Traits can be different types: numeric, ordinal, nominal)
                           traits_weights = NULL, # default = same weight for all traits
                           nbdim = 14, # default = 7, max number of dimensions
                           metric = "Gower", # other option is "Euclidean" for cluster::daisy dissimilarity matrix calculation; when using Gower's distance, number of dimensions tested should be lower than number of species
                           dendro = FALSE, # whether the best functional dendrogram should be looked for. default = FALSE
                           plot = "DMB_quality_funct_space") # set the name of the jpeg file for plots illustrating the quality of functional space. NA means no plot

# low meanSD = high quality space
round(qfs$meanSD, 4)

# keep coordinates on 4 dimensions, where meanSD < 0.004
# WHY < 0.004??
fd.coord.sgd <- qfs$details_funct_space$mat_coord[,1:4] #%>% # keeps PC1-4

write.csv(fd.coord.sgd, here("Data","FE_4D_coord_dmb.csv"))

# see variance explained by the PCoA axes
gower <- qfs$details_funct_space$mat_dissim

fit <- cmdscale(gower, eig = TRUE, k = 4) # PCoA

# variance explained by the axes
cumsum(fit$eig[fit$eig >= 0]) / sum(fit$eig[fit$eig > 0])


# Functional Richness

## Calculate species and function entity richness

# Species Richness within quadrats
sprich <- myspecies %>%
  pivot_longer(cols = Turf:ncol(.), names_to = "Taxa", values_to = "cover") %>% # reset to pivot by SGD
  distinct() %>%
  filter(cover > 0) %>%
  mutate(CowTagID = as_factor(CowTagID)) %>%
  dplyr::count(CowTagID, Taxa) %>%
  group_by(CowTagID) %>%
  summarise(counts = sum(n)) %>%
  left_join(quad.label)


## Percent Cover
# Set categories of relative SGD (based on clustering or PERMANOVA from biogeochemistry)
# OR Set categories of CowTagIDs (will sort later based on biogeochemistr or nutrients)
sgd.sp <- myspecies %>%
  left_join(quad.label) %>%
  pivot_longer(cols = Turf:ncol(.), names_to = "Taxa", values_to = "cover") %>% # reset to pivot by SGD
  group_by(CowTagID, Taxa) %>% # make sure cover is totaled by taxa
  mutate(cover = sum(cover)) %>%
  ungroup() %>%
  distinct() %>% # remove duplicate values after summing up percent cover
  group_by(CowTagID) %>%
  pivot_wider(names_from = Taxa, values_from = cover) %>%   # pivot species wide by relative SGD
  # remove seep from convex hull calculations
  filter(CowTagID != "VSEEP")
relative.sgd <- relative.sgd[1:19]

# CowTagIDs as rownames
sgd.sp <- column_to_rownames(.data = sgd.sp, var = "CowTagID")
sgd.sp <- as.data.frame(sgd.sp)

# species names as rownames
species_entities <- column_to_rownames(.data = species_entities, var = "Taxa")
species_entities <- as.data.frame(species_entities)


## Calculate convex hull (modified Teixido script)

#### Calculate convex hull (Teixido script after I gave up at convhulln)

Fric <- lapply(relative.sgd, function (x) {

  species.sgd <- colnames(sgd.sp)[which(sgd.sp[x,] > 0)]

  fes_cond.sgd <- species_entities[rownames(species_entities) %in% species.sgd, ]

  m.sgd <- fd.coord.sgd[rownames(fd.coord.sgd) %in% fes_cond.sgd, ]

  ch.sgd <- convhulln(m.sgd, options = "FA") # FA: generalized areas and volumes

  chg.sgd <- convhulln(fd.coord.sgd, options = "FA")

  c(length(species.sgd), length(species.sgd)/ncol(sgd.sp)*100, dim(m.sgd)[1], dim(m.sgd)[1]/dim(fd.coord.sgd)[1]*100, ch.sgd$vol/chg.sgd$vol*100)
  #  72 is Teixido's total number of species, so I am dividing by my total number of species (51 at Varari)

})#eo lapply


names(Fric) = relative.sgd

# Fric contains the number of species(NbSp) and FEs (NbFEs), relative percentages (NbSpP,NbFEsP ) , and the volume along the survey sites
Fric <- do.call(rbind, Fric)

colnames(Fric) <- c("NbSp", "NbSpP", "NbFEs","NbFEsP", "Vol8D")
# Number of species per condition
# percentage of species present per condition
# number of functional entities present per condition
# percentage of functional entities present per condition
# percentage of functional space volume occupied

Fric <- rownames_to_column(as.data.frame(Fric), var = "CowTagID")


## CALCULATE VALUES FOR VSEEP AND RBIND
Fric
SeepFric <- tibble(CowTagID = "VSEEP",
                   NbSp = 1,
                   NbSpP = 1/51*100,
                   NbFEs = 1,
                   NbFEsP = 1/28*100,
                   Vol8D = 0)
Fric <- Fric %>%
  rbind(SeepFric)

write_csv(Fric, here("Data", "Sp_FE_Vol.csv"))


## CALCULATE RESIDUALS AND SAVE CSV

resFric <- Fric %>%
  left_join(meta) %>%
  mutate(meanRugosity = if_else(CowTagID == "VSEEP", 0.97, meanRugosity)) %>%
  arrange(CowTagID) %>%
  mutate(meanRugosity = 1 - meanRugosity) # such that higher rugosity relates to more structural heterogeneity

# calculate residuals and join together
resSp <- residuals(lm(data = resFric, NbSp ~ meanRugosity))
resSpp <- residuals(lm(data = resFric, NbSpP ~ meanRugosity))
resFE <- residuals(lm(data = resFric, NbFEs ~ meanRugosity))
resFEp <- residuals(lm(data = resFric, NbFEsP ~ meanRugosity))
resVol <- residuals(lm(data = resFric, Vol8D ~ meanRugosity))
resFric <- resFric %>%
  cbind(resSp) %>%
  cbind(resSpp) %>%
  cbind(resFE) %>%
  cbind(resFEp) %>%
  cbind(resVol) %>%
  select(CowTagID, NbSp, NbSpP, NbFEs, NbFEsP, Vol8D,
         resSp, resSpp, resFE, resFEp, resVol)



write_csv(resFric, here("Data", "Sp_FE_Vol_res.csv"))


