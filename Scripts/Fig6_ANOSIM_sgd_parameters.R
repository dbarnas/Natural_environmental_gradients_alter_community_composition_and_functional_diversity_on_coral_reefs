#### Analysis of Similarity (ANOSIM) comparing species or functional entity community composition along SGD variability gradient

### Created by Danielle Barnas
### Created on July 1, 2023


# ANOSIM

#############################
# LOAD LIBRARIES
#############################
library(tidyverse)
library(here)
library(vegan)
library(patchwork)
# library(geosphere)





#############################
# LOAD DATA
#############################
chem <- read_csv(here("Data","Biogeochem","Nutrients_Processed_All.csv")) %>%
  filter(CowTagID != "V13",
         CowTagID != "VSEEP")

alphatag <- read_csv(here("Data", "CowTag_to_AlphaTag.csv"))

myorder <- chem %>%
  select(CowTagID, Phosphate_umolL) %>%
  left_join(alphatag) %>%
  arrange(Phosphate_umolL)
myorder <- myorder$AlphaTag



ab.sgd <- read_csv(here("Data", "Species_Abundances_wide.csv")) %>%
  filter(CowTagID != "VSEEP")
ab.sgd <- ab.sgd %>%
  left_join(alphatag) %>%
  mutate(AlphaTag = factor(AlphaTag, levels = myorder)) %>%
  arrange(AlphaTag) %>%
  dplyr::select(-AlphaTag)

mychem <- chem %>%
  select(CowTagID,
         NN_umolL,
         Phosphate_umolL,
         Salinity,
         Silicate_umolL,
         Temperature,
         pH) %>%
  left_join(alphatag) %>%
  mutate(AlphaTag = factor(AlphaTag, levels = myorder)) %>%
  arrange(AlphaTag) %>%
  dplyr::select(-AlphaTag)

# move tag names to rownames and make data.frame class as needed
mychem <- column_to_rownames(mychem,var = "CowTagID")
ab.sgd <- as.data.frame(column_to_rownames(ab.sgd, var = 'CowTagID'))

# color plot by % SR
resFric <- read_csv(here("Data", "Sp_FE_Vol_res.csv")) %>%
  filter(CowTagID != "VSEEP") %>%
  left_join(alphatag) %>%
  mutate(AlphaTag = factor(AlphaTag, levels = myorder)) %>%
  arrange(AlphaTag)


set.seed(7)




#############################
# Normalize before scaling
#############################
hist(log(mychem$Salinity + 0.0001))
hist(mychem$Silicate_umolL)
hist(mychem$NN_umolL)
hist(mychem$Phosphate_umolL)
hist(mychem$Temperature)
hist(log(mychem$pH + 0.0001))

mychem <- mychem %>%
  mutate(Salinity = log(Salinity + 0.0001),
         pH = log(pH + 0.0001))


#############################
# ANOSIM - grouped by AlphaTag
#############################

#anosim(ab.sgd, mygroup$rel, distance = "bray", permutations = 9999)




#############################
# ANOSIM - along Phosphate
#############################
# along continuous parameter: Mantel Test

### Now we have to convert these subsets into distance matrices.

#abundance data frame - bray-curtis dissimilarity
dist.abund = vegdist(ab.sgd, method = "bray")

#environmental vector - euclidean distance
# need to scale chem data
dist.chem = dist(scale(x = mychem, scale = T, center = T), method = "euclidean")

df.resFric<-as.data.frame(column_to_rownames(resFric, var = 'CowTagID')) %>%
  select(NbSpP)
dist.var = dist(df.resFric, method = "euclidean") # if I want to color points by resSpp

### Now we can run the mantel command:

#diversity vs environmental
abund_chemSR = mantel(dist.abund, dist.chem, method = "spearman", permutations = 9999, na.rm = TRUE)
abund_chemSR # p > 0.05

#### The point of this figure will be to visualize the correlation between two corresponding matrices of data. Each point in these pairwise scatter plots will represent the difference between two samples

# Next, I need to convert these distance matrices (multiple columns) into vectors (one column), and then combine these matrices into one new data frame mat

aa = as.vector(dist.abund)
tt = as.vector(dist.chem)
gg = as.vector(dist.var)

#new data frame with vectorized distance matrices
matSR = data.frame(aa,tt,gg)


#First is the pairwise comparison of community dissimilarity and differences in temperature:

#abundance vs temperature
mmSR = ggplot(matSR, aes(y = aa, x = tt)) +
  geom_point(size = 3, alpha = 0.5, shape = 21, aes(fill = gg)) +
  labs(x = expression(Delta*" SGD"), y = "Bray-Curtis Dissimilarity", fill = "% SR Dissimilarity") + #fill = expression("PO"[4]^"3-")) +
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12),
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"),
         axis.title= element_text(face = "bold", size = 14, colour = "black"),
         panel.background = element_blank(),
         panel.border = element_rect(fill = NA, colour = "black"),
         legend.position = "top",
         legend.text = element_text(size = 10, face = "bold"),
         legend.title = element_text(size = 11, face = "bold")) +
  scale_fill_continuous(high = "purple", low = "lightpink")
mmSR




####### SAME TEST BUT NOW WITH FE
#############################
# LOAD DATA
#############################
spe_fes.sgd <- as.data.frame(read_csv(here("Data", "Species_FE.csv")))

ab.sgd <- read_csv(here("Data", "Species_Abundances_wide.csv")) %>% filter(CowTagID != "VSEEP")

# join species abundance with FE to get Fe abundance
ab.sgd <- ab.sgd %>%
  pivot_longer(cols = Turf:ncol(.), names_to = "Taxa", values_to = "pCover") %>%
  left_join(spe_fes.sgd) %>%
  group_by(CowTagID, FE) %>%
  summarise(pCover = sum(pCover)) %>%
  pivot_wider(names_from = FE, values_from = pCover)

# order cowtagid's
ab.sgd <- ab.sgd %>%
  left_join(alphatag) %>%
  mutate(AlphaTag = factor(AlphaTag, levels = myorder)) %>%
  arrange(AlphaTag) %>%
  dplyr::select(-AlphaTag)

# environmental matrix for changing SGD
mychem # same as above

# move tag names to rownames and make data.frame class
ab.sgd <- as.data.frame(column_to_rownames(ab.sgd, var = 'CowTagID'))
#mychem <- as.data.frame(column_to_rownames(mychem, var = 'CowTagID'))


set.seed(7)




#############################
# ANOSIM - grouped by AlphaTag
#############################

# ref: anosim(m_com, pc$Time, distance = "bray", permutations = 9999)
#anosim(ab.sgd, mygroup$rel, distance = "bray", permutations = 9999)




#############################
# ANOSIM - along Phosphate
#############################
# along continuous parameter: Mantel Test

### Now we have to convert these subsets into distance matrices.

#abundance data frame - bray-curtis dissimilarity
dist.abund = vegdist(ab.sgd, method = "bray")

#environmental vector - euclidean distance
# need to scale chem data
#dist.chem = dist(scale(x = mychem, scale = T, center = T), method = "euclidean")
#dist.chem = dist(mychem, method = "euclidean")


# third vector variable
df.resFric <- as.data.frame(column_to_rownames(resFric, 'CowTagID')) %>% select(NbFEs)
dist.var = dist(df.resFric, method = "euclidean")

### Now we can run the mantel command:

#abundance vs environmental matrix
# spearman correlation is a nonparametric test, not relying on any particular pattern assumptions
abund_chemFER = mantel(dist.abund, dist.chem, method = "spearman", permutations = 9999, na.rm = TRUE)
abund_chemFER # p < 0.05



# Pairwise scatterplot

#### The point of this figure will be to visualize the correlation between two corresponding matrices of data. Each point in these pairwise scatter plots will represent the difference between two samples

# Next, I need to convert these distance matrices (multiple columns) into vectors (one column), and then combine these matrices into one new data frame mat

aa = as.vector(dist.abund)
tt = as.vector(dist.chem)
gg = as.vector(dist.var)

#new data frame with vectorized distance matrices
matFER = data.frame(aa,tt,gg)


#First is the pairwise comparison of community dissimilarity and differences in temperature:

#abundance vs temperature
mmFER = ggplot(matFER, aes(y = aa, x = tt)) +
  geom_point(size = 3, alpha = 0.5, shape = 21, aes(fill = gg)) +
  labs(x = expression(Delta*" SGD"), y = "Bray-Curtis Dissimilarity", fill = "% FER Dissimilarity") + #fill = expression("PO"[4]^"3-")) +
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12),
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"),
         axis.title= element_text(face = "bold", size = 14, colour = "black"),
         panel.background = element_blank(),
         panel.border = element_rect(fill = NA, colour = "black"),
         legend.position = "top",
         legend.text = element_text(size = 10, face = "bold"),
         legend.title = element_text(size = 11, face = "bold")) +
  scale_fill_continuous(high = "navy", low = "skyblue") +
  geom_smooth(method = "lm", color= "black")
mmFER



## PATCH PLOTS
Bray_Curtis_Plot <- mmSR / mmFER #+
  #plot_annotation(tag_levels = 'A')
ggsave(here("Output", "PaperFigures", "Fig6_Composition_SGD.png"), Bray_Curtis_Plot, device = "png", height = 6, width = 6)

Bray_Curtis_Plot


