### Multivariate GLM analysis of environmental parameter influence on taxonomic and functional entity diversity

### Created by Danielle Barnas
### Created on January 13, 2025
### Modified on February 23, 2025

###############################
# LOAD LIBRARIES
###############################
library(tidyverse)
library(here)
library(stats) #MANOVA
library(vegan) #nMDS
library(pairwiseAdonis) #post-hoc permanova
library(PNWColors)
library(ggordiplots)
library(BiodiversityR)
library(mgcv)
library(patchwork)

###############################
# READ IN DATA
###############################
comp <- read_csv(here("Data", "Species_Abundances_wide.csv"))
meta <- read_csv(here("Data", "Full_Metadata.csv"))
envi <- read_csv(here("Data", "Biogeochem", "Nutrients_Processed_All.csv"))
fe.traits <- read_csv(here("Data", "Distinct_Taxa_FEtraits.csv")) %>% relocate(FE, .after = Taxa)

###############################
# PROCESS AND JOIN DATA
###############################
## metadata
meta <- meta %>%
  select(Location, CowTagID, AlphaTag, complexity)

## environmental data
envi <- envi %>%
  select(Location:NN_umolL) %>%
  select(-TA)

## join
envi <- envi %>%
  left_join(meta) %>%
  filter(CowTagID != "V13",
         CowTagID != "VSEEP")

## compositional data
comp <- comp %>%
  left_join(meta[,2:3]) %>%
  filter(CowTagID != "V13",
         CowTagID != "VSEEP") %>%
  relocate(AlphaTag, .after = CowTagID) %>%
  arrange(AlphaTag)

## functional entity data
fe.traits <- fe.traits %>%
  mutate(broad_group = case_when(Taxon_Group == "Rhodophyta" ~ "Macroalgae",
                                 Taxon_Group == "Chlorophyta" ~ "Macroalgae",
                                 Taxon_Group == "Phaeophyta" ~ "Macroalgae",
                                 Taxon_Group == "Cnidaria" ~ "Coral",
                                 Taxon_Group == "Turf" ~ "Turf",
                                 Taxon_Group == "Porifera" ~ "Porifera",
                                 Taxon_Group == "Cyanobacteria" ~ "Cyanobacteria"),
         broad_group = if_else(Taxa == "Lithophyllum kotschyanum", "CCA",
                               if_else(Taxa == "Crustose Corallines", "CCA",
                                       if_else(Taxa == "Heteractis magnifica", "Other cnidaria",
                                               if_else(Taxa == "Discosoma nummiforme", "Other cnidaria", broad_group)))))

## categorize species by coral, macroalgae, and other
comp.long <- comp %>%
  pivot_longer(cols = 3:ncol(.),
               names_to = "Taxa",
               values_to = "cover") %>%
  left_join(fe.traits)

coral.data <- comp.long %>%
  filter(Taxon_Group == "Cnidaria") %>%
  filter(Taxa != "Heteractis magnifica" &
           Taxa != "Discosoma nummiforme") # story coral only
# mutate(Group = "Coral")

ma.data <- comp.long %>%
  filter(Taxon_Group == "Phaeophyta" |
           Taxon_Group == "Rhodophyta" |
           Taxon_Group == "Chlorophyta") %>%
  filter(Taxa != "Crustose Corallines" &
           Taxa != "Lithophyllum kotschyanum" &
           Taxa != "Turf")  # fleshy macroalgae only
# mutate(Group = "Macroalgae")

comp.long <- comp.long %>%
  select(CowTagID, AlphaTag, broad_group, Group = Taxa, cover) %>%
  arrange(AlphaTag)

cma.data <- rbind(coral.data, ma.data)




###############################
# PREPEARE FOR MULTIVARIATE
###############################
# scale to starndardize and calculate z-scores for all parameters
envi <- envi %>%
  relocate(AlphaTag, .after = CowTagID) %>%
  arrange(AlphaTag) #%>%
# mutate_at(.vars = vars(Salinity:complexity), .funs = ~scale(., center = TRUE, scale = TRUE))

envi.meta <- envi %>%
  select(Location:AlphaTag)
envi.only <- envi %>%
  select(-c(Location, CowTagID)) %>%
  rename(`Phosphate` = Phosphate_umolL,
         `Silicate` = Silicate_umolL,
         `N+N` = NN_umolL,
         `Complexity` = complexity)

df <- full_join(comp, envi.only)

envi.only <- column_to_rownames(envi.only, var = "AlphaTag")



###############################
# nMDS of community composition
###############################
# nMDS
comp$AlphaTag <- factor(comp$AlphaTag)
CTlevels <- unique(comp$AlphaTag)

ord1<-metaMDS(comp[,3:ncol(comp)],k=2, distance='bray')
ord1$stress # < 0.3
stressplot(ord1)



# get points for species
Group <- rownames(ord1$species) # get characteristic names
MDS1 <- c(ord1$species[,1]) # MDS1 for characteristics
MDS2 <- c(ord1$species[,2]) # MDS2 for characteristics
Data <- as_tibble(cbind(Group, MDS1, MDS2)) %>%  # bind all cols into tibble
  mutate(MDS1 = as.numeric(MDS1), # as numeric
         MDS2 = as.numeric(MDS2)) %>%
  #mutate(Taxon_Group = if_else(Taxa == "Hard Substrate", "Abiotic", Taxon_Group)) %>%
  select(MDS1, MDS2, Group) %>%
  left_join(comp.long[,3:4]) %>%
  distinct() %>%
  mutate(broad_group = factor(broad_group, levels = c("Coral", "CCA", "Macroalgae", "Turf", "Porifera", "Cyanobacteria", "Other cnidaria")))

# get points cowtagid's
Groupb <- as.character(CTlevels) # assign CowTagID
MDS1b <- ord1$points[,1] # MDS1 for CowTagID
MDS2b <- ord1$points[,2] # MDS2 for CowTagID
Datab <- as_tibble(cbind(Groupb, MDS1b, MDS2b)) %>%  # bind all cols into tibble
  mutate(MDS1b = as.numeric(MDS1b), # as numeric
         MDS2b = as.numeric(MDS2b)) %>%
  rename(AlphaTag = Groupb)

joinDF <- envi %>%
  left_join(meta) %>%
  select(CowTagID, AlphaTag, complexity, Salinity:NN_umolL)

Datab <- Datab %>%
  left_join(joinDF)

### palette
pal <- pnw_palette("Bay", n = length(unique(Data$broad_group)))


# Add loadings from environmental data to nMDS plot
envi.only.red <- envi.only[,] # remove TA 3
en <- envfit(ord1, envi.only.red, permutations = 999, na.rm = TRUE) # envfit: linear relationships only
en.red <- envfit(ord1, envi.only.red[,-ncol(envi.only.red)], permutations = 999, na.rm = TRUE) # removing complexity
en
en.red

en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)
en_coord_cont2 = as.data.frame(scores(en.red, "vectors")) * ordiArrowMul(en)

sp.mds.pcover <- ggplot(data = Data,
                        aes(x = MDS1,
                            y = MDS2)) +
  geom_point(aes(color = broad_group),
             size = 2) +
  geom_segment(data = en_coord_cont, aes(x = 0, y = 0,
                                         xend = NMDS1*2, yend = NMDS2*2),
               size =1, alpha = 0.5,
               colour = "grey30",
               arrow = arrow(angle = 30, length = unit(0.1, "inches"))) +
  geom_text(data = en_coord_cont, aes(x = NMDS1*2, y = NMDS2*2), colour = "grey30",
            fontface = "bold", label = row.names(en_coord_cont)) +
  coord_cartesian(xlim = c(-2, 2), ylim = c(-1, 1)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_color_manual(values = pal) +
  labs(color = "Category")

# ggsave(here("Output", "sp.mds.pcover.pdf"), sp.mds.pcover, device = "pdf", width = 7, height = 6)



# ordiplot theme elements
BioR.theme <- theme(
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_line("gray25"),
  text = element_text(size = 12),
  axis.text = element_text(size = 10, colour = "gray25"),
  axis.title = element_text(size = 14, colour = "gray25"),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 14),
  legend.key = element_blank())



## ggordisurf
nnord.surface <- ordisurf(ord1, y = envi.only.red$`N+N`)
summary.gam(nnord.surface)
# N+N, Silicate, Temperature

nnord.grid <- ordisurfgrid.long(nnord.surface)

plot2 <- ordiplot(ord1, choices = c(1,2))
sites.long2 <- sites.long(plot2, env.data = envi.only.red)
species.long2 <- species.long(plot2)
species.long2 <- rownames_to_column(species.long2, var = "Taxa") %>%
  left_join(fe.traits) %>%
  # modify species text
  mutate(Taxa = if_else(Taxa == "Crustose Corallines", " CCA",
                        if_else(Taxa == "Turf", " Turf",
                                if_else(Taxa == "Cyanobacteria", " Cyanobacteria", Taxa)))) %>%
  separate(Taxa, into = c("Genus", "Species"), sep = " ") %>%
  mutate(Genus = substr(Genus, start = 1, stop = 1)) %>%
  mutate(Genus = paste0(Genus,".")) %>%
  unite(col = "Taxa", Genus, Species, sep = "") %>%
  mutate(Taxa = if_else(Taxa == ".CCA", "CCA",
                        if_else(Taxa == ".Turf", "Turf",
                                if_else(Taxa == ".Cyanobacteria", "Cyanobacteria",
                                        if_else(Taxa == "D.sp.", "Dictyosphaeria sp.",
                                                if_else(Taxa == "V.sp.", "Verongida sp.", Taxa)))))) %>%
  mutate(broad_group = factor(Taxon_Group, levels = c("Cnidaria", "Porifera", "Chlorophyta", "Phaeophyta", "Rhodophyta", "Cyanobacteria", "Turf")))
axis.long2 <- axis.long(ord1, choice=c(1,2))

# prepare Taxa in species.long2 as italic species names
species.long2 <- species.long2 %>%
  mutate(pre_taxa = "italic('",
         post_taxa = "')") %>%
  unite(col = "taxa_ital", pre_taxa, Taxa, sep = "") %>%
  unite(col = "taxa_ital2", taxa_ital, post_taxa, sep = "") %>%
  rename(Taxa = taxa_ital2) %>%
  mutate(Taxa =
           if_else(Taxa == "italic('Turf')", "Turf",
                   if_else(Taxa == "italic('CCA')", "CCA",
                           if_else(Taxa == "italic('Dictyosphaeria sp.')", "italic('Dictyosphaeria')~sp.",
                                   if_else(Taxa == "italic('Verongida sp.')", "italic('Verongida')~sp.",
                                           if_else(Taxa == "italic('P.unknown')", "italic('Porifera')~unk.",
                                                   Taxa))))))


sp.ordi.nn.site <- ggplot() +
  geom_contour_filled(data = nnord.grid,
                      aes(x = x, y = y, z = z)) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  xlab(axis.long2[1, "label"]) +
  ylab(axis.long2[2, "label"]) +
  scale_x_continuous(sec.axis = dup_axis(labels = NULL, name = NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels = NULL, name = NULL)) +
  geom_point(data = sites.long2,
             aes(x = axis1, y = axis2),
             colour = "black",
             size = 4) +
  BioR.theme +
  scale_fill_viridis_d(alpha = 0.7) +
  labs(fill = "N+N") +
  coord_fixed(ratio = 1)

sp.ordi.nn.species <- ggplot() +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  xlab(axis.long2[1, "label"]) +
  ylab(axis.long2[2, "label"]) +
  scale_x_continuous(sec.axis = dup_axis(labels = NULL, name = NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels = NULL, name = NULL)) +
  geom_point(data = species.long2,
             aes(x = axis1, y = axis2,
                 color = Taxon_Group),
             # colour = "red",
             size = 1.5) + #3
  ggrepel::geom_text_repel(data = species.long2,
                           aes(x = axis1, y = axis2, label = Taxa),
                           max.overlaps = 15,
                           size = 3,
                           parse = TRUE) +
  BioR.theme +
  scale_color_manual(values = pal) +
  labs(color = "Phyla") +
  # coord_fixed(ratio = 1)
coord_cartesian(xlim = c(-1.5, 1.4), ylim = c(-1,0.8)) # use when saving separately

ggsave(here("Output", "Fig5a_sp.ordi.nn.site.png"), sp.ordi.nn.site, device = "png", width = 7, height = 7)
ggsave(here("Output", "Fig5b_sp.ordi.nn.species.png"), sp.ordi.nn.species, device = "png", width = 7, height = 5) # height=5

# layout <- '
# AB
# '
# p1 <- sp.ordi.nn.site + sp.ordi.nn.species +
#   plot_layout(design = layout)
# ggsave(here("Output", "sp.ordi.nn.png"), p1, device = "png", width = 12, height = 7)

