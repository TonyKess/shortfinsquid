setwd("~/Desktop/Projects/SQUID_popgen/")
library(tidyverse)
library(data.table)
library(qvalue)
library(vegan)
library(ggman)
library(marmap)
library(pcadapt)


#manage Squid meta
Squidmeta <- fread("Squid_Meta_Sept2024_Simple") 
Squidmetanames <- colnames(Squidmeta)
Squidfam <- fread("SQUIDKING/SQUID.tkpipeline.80geno.HWE.fam") %>%  mutate(ID = str_replace(V1, "NS.*i5.", "")) %>% 
  mutate(ID = str_replace(ID, ".real.*", "")) %>%   
  inner_join(Squidmeta) 


#map
bathydata <- getNOAA.bathy(-76,-43,35,57, res=3,keep=T)
plot(bathydata)
map=autoplot.bathy(bathydata, geom=c("r", "c")) +
  scale_fill_etopo(guide = FALSE) +
  geom_contour(aes(z=z), breaks=c(-100, -200, -500, -1000, -2000, -4000), colour="grey90", size=0.01) +
  xlab("Degrees Longitude") +
  ylab("Degrees Latitude") +
  theme(legend.position = "right")

Squidsimplemta <- Squidfam %>%  select(Lat, Lon, Zone) %>%  distinct()

map + geom_point(data = Squidsimplemta, aes(x = Lon, y  = Lat, colour = Zone, size = 3))

#PCA
squidpcaobject <- read.pcadapt(input = "SQUIDKING/SQUID.tkpipeline.80geno.HWE.bed", type = "bed")
SquidPCA <- pcadapt(squidpcaobject, K = 10)

#var explained
100 * (SquidPCA$singular.values ^ 2)

#PCA scores for mapping
PCAscores <- SquidPCA$scores
colnames(PCAscores) <- paste0("PC", rep(1:10))

#and umap
Squidumap <- umap::umap(PCAscores)
Squidumap <-Squidumap$layout
colnames(Squidumap) <- c("UMAP1", "UMAP2")

SquidPCAmeta <- bind_cols(Squidfam, PCAscores, Squidumap)

#get depths too
SquidPCAmeta <- fread("SQUID.tkpipeline.idepth") %>% 
  mutate(V1 = INDV) %>% 
  inner_join(SquidPCAmeta)

#plots - PCA
ggplot() + geom_point(data = SquidPCAmeta, aes(x = PC1, y = PC2, colour = Zone)) + theme_classic() + labs(x = "PC1 - 1.28%", y = "PC2 - 0.40%")
ggplot() + geom_point(data = SquidPCAmeta, aes(x = PC1, y = PC2, colour = as.factor(Year))) + theme_classic()
ggplot() + geom_point(data = SquidPCAmeta, aes(x = PC3, y = PC4, colour = Zone)) + theme_classic()

ggplot() + geom_point(data = SquidPCAmeta, aes(x = PC1, y = PC2, colour = MEAN_DEPTH)) + theme_classic() + theme_classic() + scale_colour_gradient(high = "red", low = "blue")
ggplot() + geom_point(data = SquidPCAmeta, aes(x = PC1, y = PC3, colour = MEAN_DEPTH)) + theme_classic() + scale_colour_gradient(high = "red", low = "blue")
ggplot() + geom_point(data = SquidPCAmeta, aes(x = PC1, y = PC2, colour =Year )) + theme_classic() 

#plots umap
ggplot() + geom_point(data = SquidPCAmeta, aes(x = UMAP1, y = UMAP2, colour = Zone)) + theme_classic()
ggplot() + geom_point(data = SquidPCAmeta, aes(x = UMAP1, y = UMAP2, colour = MEAN_DEPTH)) + theme_classic()  + theme_classic() + scale_colour_gradient(high = "red", low = "blue")

#depth by zones
ggplot() + geom_point(data = SquidPCAmeta, aes(x = Zone, y = PC1, colour = MEAN_DEPTH)) + theme_classic() + scale_colour_gradient(high = "red", low = "blue")

#admix data
SQUIDadmix <- fread("admix/SQUID.tkpipeline.80geno.pruned.2.Q") %>%  mutate(Q1 = V1, Q2 = V2) %>% select(Q1, Q2)

SquidPCAmeta.admix <- bind_cols(SquidPCAmeta, SQUIDadmix)

ggplot() + geom_point(data = SquidPCAmeta.admix, aes(x = PC1, y = Q1, colour = MEAN_DEPTH))  + theme_classic() + scale_colour_gradient(high = "red", low = "blue")
ggplot() + geom_point(data = SquidPCAmeta.admix, aes(x = PC1, y = Q2, colour = MEAN_DEPTH))  + theme_classic() + scale_colour_gradient(high = "red", low = "blue")

#Diversity
SQUIDdiv  <-fread("SQUD.tkpipeline.80geno.het") 
min(SQUIDdiv$F)
SQUIDhet <- fread("SQUD.tkpipeline.80geno.hwe")
mean(SQUIDhet$`O(HET)`)
mean(SQUIDhet$`E(HET)`)

