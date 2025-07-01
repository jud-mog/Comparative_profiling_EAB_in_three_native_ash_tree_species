# FIGURE 1 DIVERSITY INDEXES ######################

##figure 1: load libraries ####
library(vegan)
library(forcats)
library(rstatix)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(upstartr)
##fungi phloem ITS diversity ####
dat <- read.csv("~/eab-paper-ny/div_its.csv", row.names = 1);dim(dat)
dat_phloem <- subset(dat,dat$Tissue=="Phloem");dim(dat_phloem)
dat_phloem$Habitat <- fct_relevel(dat_phloem$Habitat, 
                                  "Green-P","Green-Pa","Green-G",
                                  "Black-P","Black-Pa","Black-G",
                                  "White-P","White-Pa","White-G")
#Chao diversity (Hill numbers)
chao_dunn <- dunn_test(Chao_obs ~ Habitat, data=dat_phloem, p.adjust.method="none")
chao_dunn <- chao_dunn %>% arrange(p.adj)
chao_dunn <- chao_dunn %>% add_xy_position(x = "Habitat")
chao_ITS_phloem <- ggboxplot(dat_phloem, x = "Habitat", y = "Chao_obs",
                             color = "Species", fill = "#f5f3f4", palette = c("#727472","#7fbc41","#D2781E"),
                             add = "jitter")+ 
  scale_y_continuous(n.breaks = 6) +
  theme_light()+
  stat_pvalue_manual(chao_dunn[1,], label = "p.adj.signif",
                     #                     y.position = c(4.50,3.90,4.40,3.70), 
                     tip.length = 0.015,hide.ns = TRUE)+
  labs(y = "Chao index (16S)", x = "")+
  theme(axis.text.x = element_text(family = "Tahoma", size = 9, face = "plain",angle = 45,hjust = 1)) +
  theme(axis.text.y = element_text(family = "Tahoma", size = 9, face = "plain")) +
  theme(axis.title.x = element_text(family = "Tahoma", size = 9, face = "plain"))+
  theme(axis.title.y = element_text(family = "Tahoma", size = 9, face = "bold"))+
  theme(plot.title = element_text(hjust = 0.5, size = 12, family = "Tahoma", face = "bold"));chao_ITS_phloem
write_plot(chao_ITS_phloem, width=4, height=3, format = "tiff",dpi=300) 
#Shannon diversity (Hill numbers)
shan_dunn <- dunn_test(Shannon_obs ~ Habitat, data=dat_phloem, p.adjust.method="none")
shan_dunn <- shan_dunn %>% arrange(p.adj)
shan_dunn <- shan_dunn %>% add_xy_position(x = "Habitat")
Shannon_ITS_phloem <- ggboxplot(dat_phloem, x = "Habitat", y = "Shannon_obs",
                                color = "Species", fill = "#f5f3f4", palette = c("#727472","#7fbc41","#D2781E"),
                                add = "jitter")+ 
  scale_y_continuous(n.breaks = 6) +
  theme_light()+
  stat_pvalue_manual(shan_dunn[4,], label = "p.adj.signif",
                     y.position = c(4.50), tip.length=0.015,#,3.90,4.40,3.70
                     hide.ns = TRUE)+
  labs(y = "Shannon diversity", x = "")+
  theme(axis.text.x = element_text(family = "Tahoma", size = 9, face = "plain",angle = 45,hjust = 1)) +
  theme(axis.text.y = element_text(family = "Tahoma", size = 9, face = "plain")) +
  theme(axis.title.x = element_text(family = "Tahoma", size = 9, face = "plain"))+
  theme(axis.title.y = element_text(family = "Tahoma", size = 9, face = "bold"))+
  theme(plot.title = element_text(hjust = 0.5, size = 12, family = "Tahoma", face = "bold"));Shannon_ITS_phloem

#Simpson diversity (Hill numbers)
sim_dunn <- dunn_test(Simpson_obs ~ Habitat, data=dat_phloem, p.adjust.method="none")
sim_dunn <- sim_dunn %>% arrange(p.adj)
sim_dunn <- sim_dunn %>% add_xy_position(x = "Habitat")
Simpson_ITS_phloem <- ggboxplot(dat_phloem, x = "Habitat", y = "Simpson_obs",
                                color = "Species", fill = "#f5f3f4", palette = c("#727472","#7fbc41","#D2781E"),
                                add = "jitter")+ 
  scale_y_continuous(n.breaks = 6) +
  theme_light()+
  stat_pvalue_manual(sim_dunn[c(3),], label = "p.adj.signif",
                     y.position = c(0.98),#1.02,1.04,0.98,1.0,1.06
                     tip.length = 0.015, hide.ns = TRUE)+
  labs(y = "Simpson diversity", x = "")+
  theme(axis.text.x = element_text(family = "Tahoma", size = 9, face = "plain",angle = 45,hjust = 1)) +
  theme(axis.text.y = element_text(family = "Tahoma", size = 9, face = "plain")) +
  theme(axis.title.x = element_text(family = "Tahoma", size = 9, face = "plain"))+
  theme(axis.title.y = element_text(family = "Tahoma", size = 9, face = "bold"))+
  theme(plot.title = element_text(hjust = 0.5, size = 12, family = "Tahoma", face = "bold"));Simpson_ITS_phloem
write.csv(as.matrix(rbind(chao_dunn, shan_dunn,sim_dunn)), "div_its_indexes.csv")
##fungi larva ITS diversity ####
dat <- read.csv("~/eab-paper-ny/div_its.csv", row.names = 1);dim(dat)
dat_larva <- subset(dat,dat$Tissue=="Larva")
dat_larva$Habitat <- fct_relevel(dat_larva$Habitat, 
                                 "Green-La","Black-La","White-La")
#chao index
chao_dunn <- dunn_test(Chao_obs ~ Habitat, data=dat_larva, p.adjust.method="none")
chao_dunn <- chao_dunn %>% arrange(p.adj)
chao_dunn <- chao_dunn %>% add_xy_position(x = "Habitat")
chao_ITS_larva <- ggboxplot(dat_larva, x = "Habitat", y = "Chao_obs",
                            color = "Species", fill = "#f5f3f4", palette = c("#727472","#7fbc41","#D2781E"),
                            add = "jitter")+ 
  scale_y_continuous(n.breaks = 6) +
  theme_light()+
  stat_pvalue_manual(chao_dunn, label = "p.adj.signif",
                     #                     y.position = c(4.50,3.90,4.40,3.70), 
                     tip.length = 0.015,hide.ns = TRUE)+
  labs(y = "Chao index (ITS)", x = "")+
  theme(axis.text.x = element_text(family = "Tahoma", size = 9, face = "plain",angle = 45,hjust = 1)) +
  theme(axis.text.y = element_text(family = "Tahoma", size = 9, face = "plain")) +
  theme(axis.title.x = element_text(family = "Tahoma", size = 9, face = "plain"))+
  theme(axis.title.y = element_text(family = "Tahoma", size = 9, face = "bold"))+
  theme(plot.title = element_text(hjust = 0.5, size = 12, family = "Tahoma", face = "bold"));chao_ITS_larva
write_plot(chao_ITS_larva, width=4, height=3, format = "tiff",dpi=300) 
#Shannon diversity (Hill numbers)
shan_dunn <- dunn_test(Shannon_obs ~ Habitat, data=dat_larva, p.adjust.method="none")
shan_dunn <- shan_dunn %>% arrange(p.adj)
shan_dunn <- shan_dunn %>% add_xy_position(x = "Habitat")
Shannon_ITS_larva <- ggboxplot(dat_larva, x = "Habitat", y = "Shannon_obs",
                               color = "Species", fill = "#f5f3f4", palette = c("#727472","#7fbc41","#D2781E"),
                               add = "jitter")+ 
  scale_y_continuous(n.breaks = 6) +
  theme_light()+
  stat_pvalue_manual(shan_dunn, label = "p.adj.signif",
                     #                     y.position = c(4.50,3.90,4.40,3.70), 
                     tip.length = 0.015,hide.ns = TRUE)+
  labs(y = "Shannon diversity", x = "")+
  theme(axis.text.x = element_text(family = "Tahoma", size = 9, face = "plain",angle = 45,hjust = 1)) +
  theme(axis.text.y = element_text(family = "Tahoma", size = 9, face = "plain")) +
  theme(axis.title.x = element_text(family = "Tahoma", size = 9, face = "plain"))+
  theme(axis.title.y = element_text(family = "Tahoma", size = 9, face = "bold"))+
  theme(plot.title = element_text(hjust = 0.5, size = 12, family = "Tahoma", face = "bold"));Shannon_ITS_larva
#Simpson diversity (Hill numbers)
sim_dunn <- dunn_test(Simpson_obs ~ Habitat, data=dat_larva, p.adjust.method="none")
sim_dunn <- sim_dunn %>% arrange(p.adj)
sim_dunn <- sim_dunn %>% add_xy_position(x = "Habitat")
Simpson_ITS_larva <- ggboxplot(dat_larva, x = "Habitat", y = "Simpson_obs",
                               color = "Species", fill = "#f5f3f4", palette = c("#727472","#7fbc41","#D2781E"),
                               add = "jitter")+ 
  scale_y_continuous(n.breaks = 6) +
  theme_light()+
  stat_pvalue_manual(sim_dunn, label = "p.adj.signif",
                     #                     y.position = c(1.02,1.04,0.98,1.0,1.06),
                     tip.length = 0.015, hide.ns = TRUE)+
  labs(y = "Simpson diversity", x = "")+
  theme(axis.text.x = element_text(family = "Tahoma", size = 9, face = "plain",angle = 45,hjust = 1)) +
  theme(axis.text.y = element_text(family = "Tahoma", size = 9, face = "plain")) +
  theme(axis.title.x = element_text(family = "Tahoma", size = 9, face = "plain"))+
  theme(axis.title.y = element_text(family = "Tahoma", size = 9, face = "bold"))+
  theme(plot.title = element_text(hjust = 0.5, size = 12, family = "Tahoma", face = "bold"));Simpson_ITS_larva
write.csv(as.matrix(rbind(chao_dunn, shan_dunn,sim_dunn)), "div_its_indexes_larva.csv")

##bacteria phloem diversity####
dat <- read.csv("~/eab-paper-ny/div_16s.csv", row.names = 1);dim(dat)
dat_phloem <- subset(dat,dat$Tissue=="Phloem");dim(dat_phloem)
dat_phloem$Habitat <- fct_relevel(dat_phloem$Habitat, 
                                  "Green-P","Green-Pa","Green-G",
                                  "Black-P","Black-Pa","Black-G",
                                  "White-P","White-Pa","White-G")
#chao index
chao_dunn <- dunn_test(Chao_obs ~ Habitat, data=dat_phloem, p.adjust.method="none")
chao_dunn <- chao_dunn %>% arrange(p.adj)
chao_dunn <- chao_dunn %>% add_xy_position(x = "Habitat")
chao_16S_phloem <- ggboxplot(dat_phloem, x = "Habitat", y = "Chao_obs",
                             color = "Species", fill = "#f5f3f4", palette = c("#727472","#7fbc41","#D2781E"),
                             add = "jitter")+ 
  scale_y_continuous(n.breaks = 6) +
  theme_light()+
  stat_pvalue_manual(chao_dunn[1,], label = "p.adj.signif",
                     #                     y.position = c(4.50,3.90,4.40,3.70), 
                     tip.length = 0.015,hide.ns = TRUE)+
  labs(y = "Chao index (16S)", x = "")+
  theme(axis.text.x = element_text(family = "Tahoma", size = 9, face = "plain",angle = 45,hjust = 1)) +
  theme(axis.text.y = element_text(family = "Tahoma", size = 9, face = "plain")) +
  theme(axis.title.x = element_text(family = "Tahoma", size = 9, face = "plain"))+
  theme(axis.title.y = element_text(family = "Tahoma", size = 9, face = "bold"))+
  theme(plot.title = element_text(hjust = 0.5, size = 12, family = "Tahoma", face = "bold"));chao_16S_phloem
write_plot(chao_16S_phloem, width=4, height=3, format = "tiff",dpi=300) 

#Shannon diversity (Hill numbers)
shan_dunn <- dunn_test(Shannon_obs ~ Habitat, data=dat_phloem, p.adjust.method="none")
shan_dunn <- shan_dunn %>% arrange(p.adj)
shan_dunn <- shan_dunn %>% add_xy_position(x = "Habitat")
Shannon_16S_phloem <- ggboxplot(dat_phloem, x = "Habitat", y = "Shannon_obs",
                                color = "Species", fill = "#f5f3f4", palette = c("#727472","#7fbc41","#D2781E"),
                                add = "jitter")+ 
  scale_y_continuous(n.breaks = 6) +
  theme_light()+
  stat_pvalue_manual(shan_dunn[1,], label = "p.adj.signif", 
                     tip.length = 0.015,hide.ns = TRUE)+
  labs(y = "Shannon diversity", x = "")+
  theme(axis.text.x = element_text(family = "Tahoma", size = 9, face = "plain",angle = 45,hjust = 1)) +
  theme(axis.text.y = element_text(family = "Tahoma", size = 9, face = "plain")) +
  theme(axis.title.x = element_text(family = "Tahoma", size = 9, face = "plain"))+
  theme(axis.title.y = element_text(family = "Tahoma", size = 9, face = "bold"))+
  theme(plot.title = element_text(hjust = 0.5, size = 12, family = "Tahoma", face = "bold"));Shannon_16S_phloem
#Simpson diversity (Hill numbers)
sim_dunn <- dunn_test(Simpson_obs ~ Habitat, data=dat_phloem, p.adjust.method="none")
sim_dunn <- sim_dunn %>% arrange(p.adj)
sim_dunn <- sim_dunn %>% add_xy_position(x = "Habitat")
Simpson_16S_phloem <- ggboxplot(dat_phloem, x = "Habitat", y = "Simpson_obs",
                                color = "Species", fill = "#f5f3f4", palette = c("#727472","#7fbc41","#D2781E"),
                                add = "jitter")+ 
  scale_y_continuous(n.breaks = 6) +
  theme_light()+
  stat_pvalue_manual(sim_dunn[c(2,4,8),], label = "p.adj.signif",#[c(3,6),]
                     y.position = c(1.04,1.08,1.04),#1.07,1.02,1.01,1.09, 1.04,1.05
                     #                     y.position = c(0.01,0.01,1.01,0.01,0.01, 0.01),
                     tip.length = 0.015, hide.ns = TRUE)+
  labs(y = "Simpson diversity", x = "")+
  theme(axis.text.x = element_text(family = "Tahoma", size = 9, face = "plain",angle = 45,hjust = 1)) +
  theme(axis.text.y = element_text(family = "Tahoma", size = 9, face = "plain")) +
  theme(axis.title.x = element_text(family = "Tahoma", size = 9, face = "plain"))+
  theme(axis.title.y = element_text(family = "Tahoma", size = 9, face = "bold"))+
  theme(plot.title = element_text(hjust = 0.5, size = 12, family = "Tahoma", face = "bold"));Simpson_16S_phloem
write.csv(as.matrix(rbind(chao_dunn, shan_dunn,sim_dunn)), "div_16S_indexes_phloem.csv")
##bacteria larva diversity ####
dat <- read.csv("~/eab-paper-ny/div_16s.csv", row.names = 1);dim(dat)
dat_larva <- subset(dat,dat$Tissue=="Larva")
dat_larva$Habitat <- fct_relevel(dat_larva$Habitat, 
                                 "Green-La","Black-La","White-La")
#chao index
chao_dunn <- dunn_test(Chao_obs ~ Habitat, data=dat_larva, p.adjust.method="none")
chao_dunn <- chao_dunn %>% arrange(p.adj)
chao_dunn <- chao_dunn %>% add_xy_position(x = "Habitat")
chao_16S_larva <- ggboxplot(dat_larva, x = "Habitat", y = "Chao_obs",
                            color = "Species", fill = "#f5f3f4", palette = c("#727472","#7fbc41","#D2781E"),
                            add = "jitter")+ 
  scale_y_continuous(n.breaks = 6) +
  theme_light()+
  stat_pvalue_manual(chao_dunn, label = "p.adj.signif",
                     #                     y.position = c(4.50,3.90,4.40,3.70), 
                     tip.length = 0.015,hide.ns = TRUE)+
  labs(y = "Chao index (16S)", x = "")+
  theme(axis.text.x = element_text(family = "Tahoma", size = 9, face = "plain",angle = 45,hjust = 1)) +
  theme(axis.text.y = element_text(family = "Tahoma", size = 9, face = "plain")) +
  theme(axis.title.x = element_text(family = "Tahoma", size = 9, face = "plain"))+
  theme(axis.title.y = element_text(family = "Tahoma", size = 9, face = "bold"))+
  theme(plot.title = element_text(hjust = 0.5, size = 12, family = "Tahoma", face = "bold"));chao_16S_larva
write_plot(chao_16S_larva, width=4, height=3, format = "tiff",dpi=300) 
#Shannon diversity (Hill numbers)
shan_dunn <- dunn_test(Shannon_obs ~ Habitat, data=dat_larva, p.adjust.method="none")
shan_dunn <- shan_dunn %>% arrange(p.adj)
shan_dunn <- shan_dunn %>% add_xy_position(x = "Habitat")
Shannon_16S_larva <- ggboxplot(dat_larva, x = "Habitat", y = "Shannon_obs",
                               color = "Species", fill = "#f5f3f4", palette = c("#727472","#7fbc41","#D2781E"),
                               add = "jitter")+ 
  scale_y_continuous(n.breaks = 6) +
  theme_light()+
  stat_pvalue_manual(shan_dunn, 
                     label = "p.adj.signif", 
                     tip.length = 0.015,hide.ns = TRUE)+
  labs(y = "Shannon diversity", x = "")+
  theme(axis.text.x = element_text(family = "Tahoma", size = 9, face = "plain",angle = 45,hjust = 1)) +
  theme(axis.text.y = element_text(family = "Tahoma", size = 9, face = "plain")) +
  theme(axis.title.x = element_text(family = "Tahoma", size = 9, face = "plain"))+
  theme(axis.title.y = element_text(family = "Tahoma", size = 9, face = "bold"))+
  theme(plot.title = element_text(hjust = 0.5, size = 12, family = "Tahoma", face = "bold"));Shannon_16S_larva
#Simpson diversity (Hill numbers)
sim_dunn <- dunn_test(Simpson_obs ~ Habitat, data=dat_larva, p.adjust.method="none")
sim_dunn <- sim_dunn %>% arrange(p.adj)
sim_dunn <- sim_dunn %>% add_xy_position(x = "Habitat")
Simpson_16S_larva <- ggboxplot(dat_larva, x = "Habitat", y = "Simpson_obs",
                               color = "Species", fill = "#f5f3f4", palette = c("#727472","#7fbc41","#D2781E"),
                               add = "jitter")+ 
  scale_y_continuous(n.breaks = 6) +
  theme_light()+
  stat_pvalue_manual(sim_dunn, label = "p.adj.signif",
                     #                     y.position = c(1.055,1.04,1.01,1.02,1.095,1.08),
                     tip.length = 0.015, hide.ns = TRUE)+
  labs(y = "Simpson diversity", x = "")+
  theme(axis.text.x = element_text(family = "Tahoma", size = 9, face = "plain",angle = 45,hjust = 1)) +
  theme(axis.text.y = element_text(family = "Tahoma", size = 9, face = "plain")) +
  theme(axis.title.x = element_text(family = "Tahoma", size = 9, face = "plain"))+
  theme(axis.title.y = element_text(family = "Tahoma", size = 9, face = "bold"))+
  theme(plot.title = element_text(hjust = 0.5, size = 12, family = "Tahoma", face = "bold"));Simpson_16S_larva
write.csv(as.matrix(rbind(chao_dunn, shan_dunn,sim_dunn)), "div_16S_indexes_larva.csv")

##figure 1: final figure ####
Fig.1 <- ggarrange(ggarrange(Shannon_ITS_phloem,"",Simpson_ITS_phloem,"",
                                    Shannon_ITS_larva,"",Simpson_ITS_larva,
                                    legend="none", ncol = 7, nrow = 1,
                                    widths = c(2,0.05,2,0.4,2,0.05,2),heights = 5,
                                    labels = c("A","","","","B","","")), "",
                          ggarrange(Shannon_16S_phloem,"",Simpson_16S_phloem,"",
                                    Shannon_16S_larva,"",Simpson_16S_larva,
                                    legend="none", ncol = 7, nrow = 1,
                                    widths = c(2,0.05,2,0.4,2,0.05,2),heights = 5,
                                    labels = c("C","","","","D","","")), 
                          nrow=3,ncol = 1,widths = 10,heights = c(5,0.2,5));Fig.1
write_plot(Fig.1, width=10, height=8, format = "tiff",dpi=300) 
#write_plot(Fig.1, width=15, height=12, format = "pdf",dpi=300) 
