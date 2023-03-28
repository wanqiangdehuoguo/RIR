## This scripts build for PCA analysis of small sample physicochemical properties .
## An example of an analysis is provided here. How to group depends on the design of the study.
## author:sangyimeng  
## time:2023-2-14

# 1.Loading R packages----
rm(list=ls())
library(tidyverse)
library(ggplot2)
library(patchwork)
library(vegan)
library(FactoMineR)

# 2.Data collation----
phy <- read.csv("sample_data/Physical_and_chemical_properties.csv") %>% 
  filter(mainVariable2 == "A") %>%
  mutate(group_name = c("vfd1","vfd2","vfd3","vfd4","vfd5","vfd6","vfd7","vfd8","vfd9","vfw1","vfw2","vfw3","vfw4","vfw5","vfw6","vfw7","vfw8","vfw9",
                        "spw1","spw2","spw3","spw4","spw5","spw6","spw7","spw8","spw9","spd1","spd2","spd3","spd4","spd5","spd6","spd7","spd8","spd9")) %>%
  select(group_name,pH,Salinity,TOC,Sulfide,TN,TP)
rownames(phy) <- phy$group_name
phy <- phy[,-1]
phy <- phy %>% mutate(Group1 = c("VF","VF","VF","VF","VF","VF","VF","VF","VF","VF","VF","VF","VF","VF","VF","VF","VF","VF",
                                  "SP","SP","SP","SP","SP","SP","SP","SP","SP","SP","SP","SP","SP","SP","SP","SP","SP","SP"),
                      Group2 = c("Dry","Dry","Dry","Dry","Dry","Dry","Dry","Dry","Dry","Wet","Wet","Wet","Wet","Wet","Wet","Wet","Wet","Wet",
                                  "Wet","Wet","Wet","Wet","Wet","Wet","Wet","Wet","Wet","Dry","Dry","Dry","Dry","Dry","Dry","Dry","Dry","Dry"))
phy.wet <- phy %>% filter(Group2 == "Wet")
phy.dry <- phy %>% filter(Group2 == "Dry")
phy.sp <- phy %>% filter(Group1 == "SP")
phy.vf <- phy %>% filter(Group1 == "VF")

# 3.ANOSIM intergroup differences----
distance.wet <- vegdist(phy.wet[,c(-7,-8)],method = 'bray')
anosim.wet <- anosim(distance.wet,phy.wet$Group1,permutations = 999)
distance.dry <- vegdist(phy.dry[,c(-7,-8)],method = 'bray')
anosim.dry <- anosim(distance.dry,phy.dry$Group1,permutations = 999)
distance.sp <- vegdist(phy.sp[,c(-7,-8)],method = 'bray')
anosim.sp <- anosim(distance.sp,phy.sp$Group2,permutations = 999)
distance.vf <- vegdist(phy.vf[,c(-7,-8)],method = 'bray')
anosim.vf <- anosim(distance.vf,phy.vf$Group2,permutations = 999)

# 4. PCA and visualization----
pca.phy.wet <- PCA(phy.wet[,c(-7,-8)], graph = FALSE)
bioplot.wet <- fviz_pca_biplot(pca.phy.wet,
                col.ind = phy.wet$Group1, palette = "jco", # Observe the measurement color
                addEllipses = TRUE, label = "var", ellipse.type = "confidence", # Confidence interval ellipse
                col.var = "black", repel = TRUE, # Line color
                legend.title = "Location",
                title = "R=0.992  P=0.001")
screen.wet = fviz_screeplot(pca.phy.wet,addlabels = TRUE,title="")

pca.phy.dry <- PCA(phy.dry[,c(-7,-8)], graph = FALSE)
bioplot.dry <- fviz_pca_biplot(pca.phy.dry,
                col.ind = phy.dry$Group1, palette = "jco",
                addEllipses = TRUE, label = "var", ellipse.type = "confidence", 
                col.var = "black", repel = TRUE, 
                legend.title = "Location",
                title = "R=0.916   P=0.001")
screen.dry = fviz_screeplot(pca.phy.dry,addlabels = TRUE,title="")

pca.phy.sp  <- PCA(phy.sp[,c(-7,-8)], graph = FALSE)
bioplot.sp <- fviz_pca_biplot(pca.phy.sp,
                col.ind = phy.sp$Group2, palette = c("#3CB371", "#CD5C5C"), 
                addEllipses = TRUE, label = "var", ellipse.type = "confidence", 
                col.var = "black", repel = TRUE, 
                legend.title = "Season",
                title = "R=0.1591  P=0.048")
screen.sp = fviz_screeplot(pca.phy.sp,addlabels = TRUE,title="")

pca.phy.vf  <- PCA(phy.vf[,c(-7,-8)], graph = FALSE)
bioplot.vf <- fviz_pca_biplot(pca.phy.vf,
                col.ind = phy.vf$Group2, palette = c("#3CB371", "#CD5C5C"), 
                addEllipses = TRUE, label = "var", ellipse.type = "confidence", 
                col.var = "black", repel = TRUE, 
                legend.title = "Season",
                title = "R=0.07922  P=0.122")
screen.vf = fviz_screeplot(pca.phy.vf,addlabels = TRUE,title="")

p <- (bioplot.wet+screen.wet)/(bioplot.dry+screen.dry)/(bioplot.sp+screen.sp)/(bioplot.vf+screen.vf)

# 5.save picture
ggsave(filename = "picture/physical.pca.png",height = 18,width = 9,dpi = 300)
ggsave(filename = "picture/pca_phy/bioplot.wet.png",plot = bioplot.wet,dpi = 600,width = 4.8,height = 4.125)
ggsave(filename = "picture/pca_phy/bioplot.dry.png",plot = bioplot.dry,dpi = 600,width = 4.8,height = 4.125)
ggsave(filename = "picture/pca_phy/bioplot.sp.png",plot = bioplot.sp,dpi = 600,width = 4.8,height = 4.125)
ggsave(filename = "picture/pca_phy/bioplot.vf.png",plot = bioplot.vf,dpi = 600,width = 4.8,height = 4.125)
ggsave(filename = "picture/pca_phy/screen.wet.png",plot = screen.wet,dpi = 600,width = 4.8,height = 4.125)
ggsave(filename = "picture/pca_phy/screen.dry.png",plot = screen.dry,dpi = 600,width = 4.8,height = 4.125)
ggsave(filename = "picture/pca_phy/screen.sp.png",plot = screen.sp,dpi = 600,width = 4.8,height = 4.125)
ggsave(filename = "picture/pca_phy/screen.vf.png",plot = screen.vf,dpi = 600,width = 4.8,height = 4.125)

