# FIG.6: PCA metabolomics  ######################
## Load libraries####
library(FactoMineR)
library(factoextra)
library(upstartr)
library(vegan)
library(pairwiseAdonis)
## Load datasets and pca analysis ####
meta <- read.csv("~/eab-paper-ny/metabo_pcaFF.csv", row.names = 1);dim(meta) 
env <- meta[, 1:3];dim(env)
meta.hel <- decostand(meta[,-(1:3)], "hellinger")
res.pca <- PCA(meta.hel, scale.unit=TRUE, ncp=5, graph=F)
contri <- res.pca$var$contrib
cols <- c("darkgreen","#7fbc41","#AFFF2A",
          "black","#727472","darkgrey",
          "#ec7014","#D2781E","#fcbba1")
Species  <-  as.factor(meta$Habitat)
Fig.6 <- fviz_pca_biplot(res.pca, label = c("var"), habillage = Species,
                           select.var = list(contrib = 10),#select.ind = list(contrib = 30),
                           palette = c("gray","black","lightgreen","#00a200","#FFE239","#FF8735"),
                           addEllipses=T, ellipse.level=0.95,ellipse.type = "confidence",
                           #                                title="Top10 contibuting to the model",
                           subtitle="F = 2.1252  p = 1e-04",
                           title = "",
                           legend.title = "Sample type",mean.point = FALSE,
                           labelsize=3,repel = T,ggtheme = theme_bw())+
  theme(plot.subtitle = element_text(color="black", size=10, face="italic",hjust = 1),
        plot.title = element_text(color="black", size=12, face="bold",hjust = 0),
        legend.title = element_text(color="black", size=12, face="bold"),
        legend.text = element_text(color="black", size=10, face="plain"),
        axis.title.x = element_text(color="black", size=10, face="bold"),
        axis.title.y = element_text(color="black", size=10, face="bold")); Fig.6_2#Fig.5A.PCA+
## test the PCA model ####
meta.dist <- vegdist(meta.hel, method = "euclidean") 
bdisp <- betadisper(meta.dist, meta$Habitat, "centroid") 
bdisp
anova(bdisp) 
PERMA = adonis2(meta.dist ~ Habitat, data = meta, 
                permutations = 9999, method = "euclidean")
PERMA
posthoc <- pairwiseAdonis::pairwise.adonis(meta.dist,factors=meta$Habitat, 
                                           p.adjust.m = "BH", perm = 99999);posthoc
