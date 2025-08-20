# FIG.5: DISCRIMINANT ANALYSES ######################
## Load libraries ####
library(microbiomeMarker)
library(phyloseq)
## identifying markers fungal/bacterial community ####
phy<- readRDS("~/eab-paper-ny/phy_its.RData");phy
tax_table(phy) <- tax_table(phy)[,1:6]
phyFP <- subset_samples(phy, sample_data(phy)$Tissue=="Phloem"); phyFP
#remove the PhloemA as no marker was found related to that group
phyFP1 <- subset_samples(phyFP, !sample_data(phyFP)$Sample_Type=="PhloemA"); phyFP1
### markers for the black ash ####
phyfungFPB <- subset_samples(phyFP1, sample_data(phyFP1)$Species=="Black");phyfungFPB
summarize_taxa(phyfungFPB)
sample_data(phyfungFPB)
black_edgerITS <- run_edger(phyfungFPB,norm="TMM",transform = "identity",
                            group = "Sample_Type",#contrast = c("Phloem","Gallery"),
                            pvalue_cutoff = 0.01,
                            p_adjust = "none")
black_edgerITS
marker_table(black_edgerITS)
black_da_barITS <- plot_ef_bar(black_edgerITS);black_da_barITS
### markers for the green ash ####
phyfungFPG <- subset_samples(phyFP1, sample_data(phyFP1)$Species=="Green");phyfungFPG
summarize_taxa(phyfungFPG)
sample_data(phyfungFPG)
green_edgerITS <- run_edger(phyfungFPG,norm="TMM",transform = "identity",
                            group = "Sample_Type",#contrast = c("Phloem","Gallery"),
                            pvalue_cutoff = 0.01,
                            p_adjust = "none")
green_edgerITS
marker_table(green_edgerITS)
green_da_barITS <- plot_ef_bar(green_edgerITS);green_da_barITS
### markers for the white ash ####
phyfungFPW <- subset_samples(phyFP1, sample_data(phyFP1)$Species=="White");phyfungFPW
summarize_taxa(phyfungFPW)
sample_data(phyfungFPW)
white_edgerITS <- run_edger(phyfungFPW,norm="TMM",transform = "identity",
                            group = "Sample_Type",#contrast = c("Phloem","Gallery"),
                            pvalue_cutoff = 0.01,
                            p_adjust = "none")
white_edgerITS
marker_table(white_edgerITS)
white_da_barITS <- plot_ef_bar(white_edgerITS);white_da_barITS
