# FIG.3: TAXONOMIC PROFILES ######################
##load libraries ####
library(microbiomeutilities)
library(phyloseq)
library(vegan)
library(forcats)
##Load data fungal/bacterial community: Phloem ####
phy<- readRDS("~/eab-paper-ny/phy_its.RData");phy
phyP <- subset_samples(phy, sample_data(phy)$Tissue=="Phloem");phyP
### Top20 (phloem, fungi)####
phyloFg.gen <- phyloseq::tax_glom(phyP, "Genus")
phyloseq::tax_table(phyloFg.gen)
#Check the names
taxa_names(phyloFg.gen)[1:5]
#set this the Taxonomic names in italics
guide_italics <- guides(fill = guide_legend(label.theme = element_text(
  size = 11,face = "italic", 
  colour = "Black", angle = 0))) 
#fix colors
mycolors <- c("#88CCEE", "#ED0DFD", "#DDCC77", "#888888", "#332288", 
              "#00CBA7","#00C2F9","#FFB935","#FFACC6", "#FF8735",
              "grey", "#0079FA", "#D55E00", "#CC79A7","#000000",
              "#00D302","#661100","#999933","#DE0D2E","#009503",
              "#D80D7B","#86FFDE","grey","#0079FA","#8E06CD")#,
#select the top20 genus 
phylo.genus <- aggregate_top_taxa2(phyloFg.gen, "Genus", top = 20)
rank_names(phylo.genus)
phyloseq::tax_table(phylo.genus)
phylo.genus.rel <- microbiome::transform(phylo.genus, "compositional")
plot.composition.relAbun <- plot_composition(phylo.genus.rel,
                                             sample.sort = "Species", 
                                             average_by = "Habitat",
                                             x.label = "Sample_Type") 
data.com <- plot.composition.relAbun$data
colnames(data.com)
colnames(data.com) <- c("Tax","Sample","Abundance","Habitat")
levels(data.com$Habitat)
data.com$Habitat <- fct_relevel(data.com$Habitat, 
                                "Green-P","Black-P","White-P",
                                "Green-Pa","Black-Pa","White-Pa",
                                "Green-G","Black-G","White-G")
top20 <- ggplot(data.com, aes(x = Habitat, y = Abundance, fill = Tax));top20
top20 <- top20 + geom_bar(position = "stack", stat = "identity");top20
top20 <- top20 + scale_x_discrete(labels = data.com$Habitat, breaks = data.com$Habitat);top20
top20 <- top20 + scale_y_continuous(labels = scales::percent)+
  theme(strip.background = element_rect(color="black", fill="white", linewidth=1, linetype="solid"), 
        strip.text.x = element_text(size = 11, color = "black", face = "bold.italic"));top20
top20 <- top20 + theme_bw()+scale_fill_manual("Genus", values = mycolors);top20
top20 <- top20 + theme(axis.text.x=element_text(angle=90, vjust=0.5,face = "plain",size = 12));top20
top20 <- top20 + labs(y= "Relative abundance (%)", x = "");top20
top20 <- top20 + theme(plot.title = element_text(family = "Tahoma", face = "bold", size = (15)),
                       legend.title = element_text(colour = "black", face = "bold", family = "Tahoma", size = (11)),
                       legend.text = element_text(face = "italic", colour = "black", family = "Tahoma", size = (4)),
                       axis.title = element_text(family = "Tahoma", size = (10), colour = "black", face = "bold"),
                       axis.text = element_text(family = "Tahoma", colour = "black", size = (9)));top20
phloemtop20.ITS <- top20 + ggtitle(" ") + guide_italics + 
  theme(legend.title = element_text(size = 11, face = "bold"));phloemtop20.ITS

### Phylum level (phloem, fungi)####
ph1 <- tax_glom(phyP, taxrank = 'Phylum') 
(ph2 = merge_samples(ph1, "Habitat")) 
ph3 <- transform_sample_counts(ph2, function(x)100* x/sum(x)) 
ph4 <- psmelt(ph3) 
ph4$Phylum <- as.character(ph4$Phylum) 
mycolors <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", 
              "#DE0D2E","#AFFF2A","#FFB935","#00F407", "#999933", "#FF8735",
              "#009FFA","#FF4235","#00B408","#DE0D2E","#009503",
              "#B40AFC","#009175","#AB0D61","#00735C","#810D49",
              "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#000000")
ph4$Sample <- fct_relevel(ph4$Sample, 
                          "Green-P","Black-P","White-P",
                          "Green-Pa","Black-Pa","White-Pa",
                          "Green-G","Black-G","White-G")
#plot
guide_text <- guides(fill = guide_legend(label.theme = element_text(
  size = 9, face = "plain", 
  colour = "Black", angle = 0))) 
phylum <- ggplot(data=ph4, aes(x=Sample, y=Abundance, fill=Phylum))+
  geom_bar(aes(), stat="identity", position="stack") + theme_bw()+
  scale_fill_manual(values=mycolors)+
  theme(legend.position="right") + guides(fill=guide_legend(nrow=5))+
  theme(axis.text.x=element_text(family = "Tahoma",angle=90, vjust=0.5,hjust=0.5,face = "plain",size = 9))+
  theme(plot.title = element_text(family = "Tahoma", face = "bold", size = (12)),
        legend.title = element_text(colour = "black", face = "bold", family = "Tahoma", size = 9),
        legend.text = element_text(colour = "black", family = "Tahoma", size = (9)),#face = "italic", 
        axis.title = element_text(family = "Tahoma", size = (9), colour = "black", face = "bold"),
        axis.text = element_text(family = "Tahoma", colour = "black", size = (9)));phylum
phloemphylumITS <- phylum + ggtitle(" ") + guide_text + 
  theme(legend.title = element_text(size = 9, face = "bold"));phloemphylumITS

### Class level (phloem, fungi)####
cla1 <- tax_glom(phyP, taxrank = 'Class') 
(cla2 = merge_samples(cla1, "Habitat")) 
cla3 <- transform_sample_counts(cla2, function(x)100* x/sum(x)) 
cla4 <- psmelt(cla3) 
cla4$Class <- as.character(cla4$Class) 
cla4$Class[cla4$Abundance < 0.1] <- "Class < 0.1% abund." 
#set color palette to accommodate the number of genera
mycolors <- c("#88CCEE", "#CC6677", "#DDCC77","#000000","#117733", "#332288", 
              "grey","#AFFF2A","#FFB935","#00F407", "#999933", "#FF8735",
              "#00C2F9","#888888","#DE0D2E","#FFACC6","#6699CC",
              "#D80D7B","#009175","#00E5F8","#00735C","#810D49",
              "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#86FFDE",
              "#0079FA","#8E06CD","#005FCC","#6B069F","#00489E",
              "#460B70","#00306F","#FFD7E1","#86FFDE","#FFACC6")

cla4$Sample <- fct_relevel(cla4$Sample, 
                           "Green-P","Black-P","White-P",
                           "Green-Pa","Black-Pa","White-Pa",
                           "Green-G","Black-G","White-G")
#plot
guide_text <- guides(fill = guide_legend(label.theme = element_text(
  size = 12, face = "plain", 
  colour = "Black", angle = 0))) 
class <- ggplot(data=cla4, aes(x=Sample, y=Abundance, fill=Class))+
  geom_bar(aes(), stat="identity", position="stack") + theme_bw()+ 
  scale_fill_manual(values=mycolors) + 
  theme(legend.position="right") + guides(fill=guide_legend(nrow=5))+
  theme(axis.text.x=element_text(family = "Tahoma",angle=90, 
                                 vjust=0.5,hjust=0.5,
                                 face = "plain",size = 10))+
  ggtitle("Class level ")+
  labs(y= "Relative abundance (%)", x = "")+
  theme(plot.title = element_text(family = "Tahoma", face = "bold", size = (15)),
        legend.title = element_text(colour = "black", face = "bold", family = "Tahoma", size = (10)),
        legend.text = element_text(colour = "black", family = "Tahoma", size = (10)),#face = "italic", 
        axis.title = element_text(family = "Tahoma", size = (10), colour = "black", face = "bold"),
        axis.text = element_text(family = "Tahoma", colour = "black", size = (9)));class
phloemclassITS <- class + ggtitle(" ") + guide_text + 
  theme(legend.title = element_text(size = 12, face = "bold"));phloemclassITS

##Load data fungal/bacterial community: Larva ####
### Top20 (larva, fungi)####
phy<- readRDS("~/eab-paper-ny/phy_its.RData");phy
phyL <- subset_samples(phy, sample_data(phy)$Tissue=="Larva");phyL
phyloFg.gen <- phyloseq::tax_glom(phyL, "Genus") 
#Check the names.
taxa_names(phyloFg.gen)[1:5]
#set this the Taxonomic names in italics
guide_italics <- guides(fill = guide_legend(label.theme = element_text(
  size = 11,face = "italic", 
  colour = "Black", angle = 0))) 
#fix colors
mycolors <- c("#88CCEE", "#ED0DFD", "#DDCC77", "#888888", "#332288", 
              "#00CBA7","#00C2F9","#FFB935","#FFACC6", "#FF8735",
              "grey", "#0079FA", "#D55E00", "#000000","#CC79A7",
              "#00D302","#661100","#999933","#DE0D2E","#009503",
              "#D80D7B","#86FFDE","grey","#0079FA","#8E06CD")#,
#select the top20 genus 
phylo.genus <- aggregate_top_taxa2(phyloFg.gen, "Genus", top = 20)
rank_names(phylo.genus)
phyloseq::tax_table(phylo.genus)
phylo.genus.rel <- microbiome::transform(phylo.genus, "compositional")
plot.composition.relAbun <- plot_composition(phylo.genus.rel,
                                             sample.sort = "Species", average_by = "Habitat",x.label = "Sample_Type") 
data.com <- plot.composition.relAbun$data
colnames(data.com)
colnames(data.com) <- c("Tax","Sample","Abundance","Habitat")
data.com$Habitat <- fct_relevel(data.com$Habitat, 
                                "Green-La","Black-La","White-La")
top20 <- ggplot(data.com, aes(x = Habitat, y = Abundance, fill = Tax));top20
top20 <- top20 + geom_bar(position = "stack", stat = "identity");top20
top20 <- top20 + scale_x_discrete(labels = data.com$Habitat, breaks = data.com$Habitat);top20
top20 <- top20 + scale_y_continuous(labels = scales::percent)+ 
  theme(strip.background = element_rect(color="black", fill="white", linewidth=1, linetype="solid"), 
        strip.text.x = element_text(size = 11, color = "black", face = "bold.italic"));top20
top20 <- top20 + theme_bw()+scale_fill_manual("Genus", values = mycolors);top20
top20 <- top20 + theme(axis.text.x=element_text(angle=90, vjust=0.5,face = "plain",size = 12));top20
top20 <- top20 + labs(y= "Relative abundance (%)", x = "");top20
top20 <- top20 + theme(plot.title = element_text(family = "Tahoma", face = "bold", size = (15)),
                       legend.title = element_text(colour = "black", face = "bold", family = "Tahoma", size = (11)),
                       legend.text = element_text(face = "italic", colour = "black", family = "Tahoma", size = (4)),
                       axis.title = element_text(family = "Tahoma", size = (10), colour = "black", face = "bold"),
                       axis.text = element_text(family = "Tahoma", colour = "black", size = (9)));top20
larvatop20.ITS <- top20 + ggtitle(" ") + guide_italics + 
  theme(legend.title = element_text(size = 11, face = "bold"));larvatop20.ITS

### Phylum level (larva, fungi)####
phyL <- subset_samples(phy, sample_data(phy)$Tissue=="Larva");phyL
ph1 <- tax_glom(phyL, taxrank = 'Phylum') 
(ph2 = merge_samples(ph1, "Habitat"))
ph3 <- transform_sample_counts(ph2, function(x)100* x/sum(x)) 
ph4 <- psmelt(ph3) 
ph4$Phylum <- as.character(ph4$Phylum) 
mycolors <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", 
              "#DE0D2E","#AFFF2A","#FFB935","#00F407", "#999933", "#FF8735",
              "#009FFA","#FF4235","#00B408","#DE0D2E","#009503",
              "#B40AFC","#009175","#AB0D61","#00735C","#810D49",
              "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#000000")
ph4$Sample <- fct_relevel(ph4$Sample, 
                          "Green-La","Black-La","White-La")
#plot
guide_text <- guides(fill = guide_legend(label.theme = element_text(
  size = 9, face = "plain", 
  colour = "Black", angle = 0))) 
phylum <- ggplot(data=ph4, aes(x=Sample, y=Abundance, fill=Phylum))+
  geom_bar(aes(), stat="identity", position="stack") + theme_bw()+
  scale_fill_manual(values=mycolors)+
  theme(legend.position="right") + guides(fill=guide_legend(nrow=5))+
  theme(axis.text.x=element_text(family = "Tahoma",angle=90, vjust=0.5,hjust=0.5,face = "plain",size = 9))+
  ggtitle("Phylum level ")+
  labs(y= "Relative abundance (%)", x = "")+
  theme(plot.title = element_text(family = "Tahoma", face = "bold", size = (12)),
        legend.title = element_text(colour = "black", face = "bold", family = "Tahoma", size = 9),
        legend.text = element_text(colour = "black", family = "Tahoma", size = (9)),#face = "italic", 
        axis.title = element_text(family = "Tahoma", size = (9), colour = "black", face = "bold"),
        axis.text = element_text(family = "Tahoma", colour = "black", size = (9)));phylum

larvaphylumITS <- phylum + ggtitle(" ") + guide_text + 
  theme(legend.title = element_text(size = 9, face = "bold"));larvaphylumITS

### Class level (larva, fungi)####
cla1 <- tax_glom(phyL, taxrank = 'Class') 
(cla2 = merge_samples(cla1, "Habitat")) 
cla3 <- transform_sample_counts(cla2, function(x)100* x/sum(x)) 
cla4 <- psmelt(cla3) 
cla4$Class <- as.character(cla4$Class) 
cla4$Class[cla4$Abundance < 0.1] <- "Class < 0.1% abund." 

#set color palette to accommodate the number of genera
mycolors <- c("#88CCEE","#000000","#005FCC","grey","#AFFF2A","#FFB935","#8E06CD",
              "#FFD7E1","#00F407","#DE0D2E")

cla4$Sample <- fct_relevel(cla4$Sample, 
                           "Green-La","Black-La","White-La")
#plot
guide_text <- guides(fill = guide_legend(label.theme = element_text(
  size = 12, face = "plain", 
  colour = "Black", angle = 0))) 
class <- ggplot(data=cla4, aes(x=Sample, y=Abundance, fill=Class))+
  geom_bar(aes(), stat="identity", position="stack") + theme_bw()+ 
  scale_fill_manual(values=mycolors) + 
  theme(legend.position="right") + guides(fill=guide_legend(nrow=5))+
  theme(axis.text.x=element_text(family = "Tahoma",angle=90, vjust=0.5,
                                 hjust=0.5,face = "plain",size = 10))+
  ggtitle("Class level ")+
  labs(y= "Relative abundance (%)", x = "")+
  theme(plot.title = element_text(family = "Tahoma", face = "bold", size = (15)),
        legend.title = element_text(colour = "black", face = "bold", family = "Tahoma", size = (10)),
        legend.text = element_text(colour = "black", family = "Tahoma", size = (10)),#face = "italic", 
        axis.title = element_text(family = "Tahoma", size = (10), colour = "black", face = "bold"),
        axis.text = element_text(family = "Tahoma", colour = "black", size = (9)));class

larvaclassITS <- class + ggtitle(" ") + guide_text + 
  theme(legend.title = element_text(size = 12, face = "bold"));larvaclassITS
