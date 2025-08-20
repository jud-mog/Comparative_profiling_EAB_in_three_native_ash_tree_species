# FIG.4: BETADIVERSITY ######################
## Load packages all packages at once ####
pcoa_packages <- c('upstartr', 'ggplot2', 'cowplot', 
                   'ggpubr', 'ade4', 'ggforce','pairwiseAdonis',
                   'pairwiseAdonis', 'tidyverse', 'patchwork',
                   'ape','phyloseq')
lapply(pcoa_packages, require, character.only=TRUE)
## Phloem samples community/bacterial community ####
phy<- readRDS("~/eab-paper-ny/phy_its.RData");phy
phyF <- subset_samples(phy, sample_data(phy)$Tissue=="Phloem");phyF
spe <- otu_table(phyF);dim(spe)
env <- as.data.frame(sample_data(phyF));dim(env)
spe <- as.data.frame(t(spe));dim(spe)
envT <- cbind(env,spe);dim(envT)
spe.hel <- decostand(envT[,13:1110], method="hellinger")
is.euclid(dist(spe.hel)) #check before moving forward
spe.dbc <- vegdist(spe.hel,"bray") # bray curtis
is.euclid(dist(spe.dbc))
spe.djac <- vegdist(spe.hel,"jac", binary = TRUE) # jaccard
is.euclid(dist(spe.djac))
#test the homogeneity of the data
#for group of interest (species or Habitat)
bdisper <- betadisper(spe.dbc, envT$Habitat, "centroid", bias.adjust=TRUE) 
bdisper 
anova(bdisper)
#test the model at the species level
PERMA = adonis2(spe.hel ~ Habitat, data = envT, permutations = 9999, method = "bray")
PERMA
pairwise.adonis(spe, envT$Habitat, p.adjust.m = "fdr", 
                         perm = 99999, sim.method = "bray")
#pcoa phloem samples
spe.pcoa <- cmdscale(spe.dbc,  k=sqrt(nrow(spe)-1), eig=TRUE)
pc1 <- round(spe.pcoa$eig[1]/sum(spe.pcoa$eig)*100);pc1
pc2 <- round(spe.pcoa$eig[2]/sum(spe.pcoa$eig)*100);pc2
pc3 <- round(spe.pcoa$eig[3]/sum(spe.pcoa$eig)*100);pc3
pc4 <- round(spe.pcoa$eig[4]/sum(spe.pcoa$eig)*100);pc4
spe.pcoa <- as.data.frame(spe.pcoa$points); dim(spe.pcoa)
spe.pcoa <- spe.pcoa[,1:3]
dim(spe.pcoa)
Species <- envT$Species
Habitat <- envT$Sample_Type
final <- cbind(spe.pcoa, Species, Habitat)
#now plot the pcoa graph
cols <- c("#727472","#7fbc41","#D2781E")
pcoa_phloem_ITS <- ggplot(final, aes(x = V1, y = V3, color = Species)) + 
  geom_point(size = 0.1) +
  geom_point(aes(shape = Habitat, stroke = 1), size=4)+  #
  scale_shape_manual(values = c(10,15,16))+ 
  scale_color_manual(values = cols)+
  stat_ellipse( aes(group=Species),level = 0.90) + 
  theme_bw() +
  guides(shape = guide_legend(title = "Sample type",size = 10,
                              theme = theme(legend.title = element_text(hjust = 1, size = 12)),
                              order = 1),colour = guide_legend(title = "Species",size = 10,
                                                               theme = theme(legend.title = element_text(hjust = 0, size = 12))))+
  xlim(-0.7, 0.7) + 
  ylim(-0.5, 0.5)+
  geom_hline(yintercept=0) + geom_vline(xintercept=0)+
  theme(strip.background = element_rect(color="black", fill="white", size=10, linetype="solid"), 
        strip.text.x = element_text(size = 10, color = "black", face = "bold.italic"))+
  labs(x=paste("PCoA 1 (", pc1, "%)", sep = ""), y=paste("PCoA 3 (", pc3, "%)", sep = ""))+
  ggtitle("F = 1.3503, p = 0.0001")+
  theme(plot.title = element_text(color="black", size=8, face="bold.italic",hjust = 1),
        axis.title.x = element_text(color="black", size=10, face="bold"),
        axis.title.y = element_text(color="black", size=10, face="bold"))
pcoa_phloem_ITS

## Larva samples community/bacterial community ####
phy<- readRDS("~/eab-paper-ny/phy_its.RData");phy
phyF1 <- subset_samples(phy, sample_data(phy)$Tissue=="Larva");phyF1
spe <- otu_table(phyF1);dim(spe)
env <- as.data.frame(sample_data(phyF1));dim(env)
spe <- as.data.frame(t(spe));dim(spe)
envT <- cbind(env,spe);dim(envT)
spe.hel <- decostand(envT[,13:1110], method="hellinger")
is.euclid(dist(spe.hel))
spe.dbc <- vegdist(spe.hel,"bray") # bray curtis
is.euclid(dist(spe.dbc))
spe.djac <- vegdist(spe.hel,"jac", binary = TRUE) # jaccard
is.euclid(dist(spe.djac))
#test the homogeneity of the data
bdisper <- betadisper(spe.dbc, envT$Species, "centroid", bias.adjust=TRUE) 
bdisper 
anova(bdisper)
#test the model at the species level
PERMA = adonis2(spe.hel ~ Species, data = envT, permutations = 9999, method = "jaccard")
PERMA
#pcoa
spe.pcoa <- cmdscale(spe.djac,  k=sqrt(nrow(spe)-1), eig=TRUE)
pc1 <- round(spe.pcoa$eig[1]/sum(spe.pcoa$eig)*100);pc1
pc2 <- round(spe.pcoa$eig[2]/sum(spe.pcoa$eig)*100);pc2
pc3 <- round(spe.pcoa$eig[3]/sum(spe.pcoa$eig)*100);pc3
pc4 <- round(spe.pcoa$eig[4]/sum(spe.pcoa$eig)*100);pc4
spe.pcoa <- as.data.frame(spe.pcoa$points); dim(spe.pcoa)
spe.pcoa <- spe.pcoa[,1:3]
dim(spe.pcoa)
Species <- envT$Species
final <- cbind(spe.pcoa, Species)
#now plot the pcoa graph
cols <- c("#727472","#7fbc41","#D2781E")
pcoa_larva_ITS <- ggplot(final, aes(x = V1, y = V2, color = Species)) + 
  geom_point(size = 0.1) +
  geom_point(aes(stroke = 1), size=4)+  #
  scale_shape_manual(values = c(1,1,1))+ 
  scale_color_manual(values = cols)+
  theme_bw() +
  guides(shape = guide_legend(title = "Sample type",size = 10,
                              theme = theme(legend.title = element_text(hjust = 1, size = 12)),
                              order = 1),colour = guide_legend(title = "Species",size = 10,
                                                               theme = theme(legend.title = element_text(hjust = 0, size = 12))))+
  xlim(-0.6, 0.6) + 
  ylim(-0.6, 0.6)+
  geom_hline(yintercept=0) + geom_vline(xintercept=0)+
  theme(strip.background = element_rect(color="black", fill="white", size=10, linetype="solid"), 
        strip.text.x = element_text(size = 10, color = "black", face = "bold.italic"))+
  labs(x=paste("PCoA 1 (", pc1, "%)", sep = ""), y=paste("PCoA 2 (", pc2, "%)", sep = ""))+
  ggtitle("F = 0.7817  p = 0.8444")+
  theme(plot.title = element_text(color="black", size=8, face="bold.italic",hjust = 1),
        axis.title.x = element_text(color="black", size=10, face="bold"),
        axis.title.y = element_text(color="black", size=10, face="bold"))
pcoa_larva_ITS
