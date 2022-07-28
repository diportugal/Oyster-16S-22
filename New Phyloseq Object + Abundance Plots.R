
library(vegan)
#install.packages("taxa")
library(taxa)
library(ggplot2)
#install.packages("ggpubr")
library(ggpubr)

#Original phyloseq object ####
physeq = phyloseq(OTU, TAX, SAMP) 
physeq


    
# CREATING NEW PHYLOSEQ OBJECTS ####
##peacrabs.f_1_vs_0
PEApa_comp <- read.csv("Reduced Tax Tables/PresAbs_PEA_COMP.csv") 
PEApa_comp <- as.matrix(PEApa_comp) 
colnames(PEApa_comp) <- c("X", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus.x","Genus.y","Species")
rownames(PEApa_comp)
PEApa_comp <- tax_table(PEApa_comp)
rownames(PEApa_comp)
sample_names(SAMP) 
sample_names(OTU)
sample_names(PEApa_comp) 
phyloPEA_PA <- phyloseq(OTU, PEApa_comp, SAMP) 
phyloPEA_PA


##RFTM_pa1.peacrabs.f1 
RFTMPEApa_comp <- read.csv("Reduced Tax Tables/PresAbs_RFTMPEA_COMP.csv") 
RFTMPEApa_comp <- as.matrix(RFTMPEApa_comp) 
colnames(RFTMPEApa_comp) <- c("X", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus.x","Genus.y","Species")
rownames(RFTMPEApa_comp)
RFTMPEApa_comp <- tax_table(RFTMPEApa_comp)
rownames(RFTMPEApa_comp)
sample_names(SAMP) 
sample_names(OTU)
sample_names(RFTMPEApa_comp) 
phyloRFTMPEA_PA <- phyloseq(OTU, RFTMPEApa_comp, SAMP) 
phyloRFTMPEA_PA


## RFTM_pa_1_vs_0
#Making a new phyloseq object for RFTM_PA
RFTMpa_comp <- read.csv("Reduced Tax Tables/PresAbs_RFTM_COMP.csv") 
#Setting new tax table 
RFTMpa_comp <- as.matrix(RFTMpa_comp) 
colnames(RFTMpa_comp) <- c("X", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus.x","Genus.y","Species")
rownames(RFTMpa_comp)
#SEQUENCE 
RFTMpa_comp <- tax_table(RFTMpa_comp)
rownames(RFTMpa_comp) #SEQUENCE
#INSPECTING SAMPLE NAMES
sample_names(SAMP) #UniqueID
sample_names(OTU) #UniqueID
sample_names(RFTMpa_comp) #NULL
phyloRFTM_PA <- phyloseq(OTU, RFTMpa_comp, SAMP) 
phyloRFTM_PA


##Site.x_OY_vs_NW
SITEOYNWpa_comp <- read.csv("Reduced Tax Tables/PresAbs_SITEOYNW_COMP.csv") 
SITEOYNWpa_comp <- as.matrix(SITEOYNWpa_comp) 
colnames(SITEOYNWpa_comp) <- c("X", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus.x","Genus.y","Species")
rownames(SITEOYNWpa_comp)
SITEOYNWpa_comp <- tax_table(SITEOYNWpa_comp)
rownames(SITEOYNWpa_comp)
sample_names(SAMP) 
sample_names(OTU)
sample_names(SITEOYNWpa_comp) 
phyloSITEOYNW_PA <- phyloseq(OTU, SITEOYNWpa_comp, SAMP) 
phyloSITEOYNW_PA


##Site.x_SW_vs_NW
SITESWNWpa_comp <- read.csv("Reduced Tax Tables/PresAbs_SITESWNW_COMP.csv") 
SITESWNWpa_comp <- as.matrix(SITESWNWpa_comp) 
colnames(SITESWNWpa_comp) <- c("X", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus.x","Genus.y","Species")
rownames(SITESWNWpa_comp)
SITESWNWpa_comp <- tax_table(SITESWNWpa_comp)
rownames(SITESWNWpa_comp)
sample_names(SAMP) 
sample_names(OTU)
sample_names(SITESWNWpa_comp) 
phyloSITESWNW_PA <- phyloseq(OTU, SITESWNWpa_comp, SAMP) 
phyloSITESWNW_PA


##RFTM_pa1.Site.xOY
RFTMOYpa_comp <- read.csv("Reduced Tax Tables/PresAbs_RFTMOY_COMP.csv") 
RFTMOYpa_comp <- as.matrix(RFTMOYpa_comp) 
colnames(RFTMOYpa_comp) <- c("X", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus.x","Genus.y","Species")
rownames(RFTMOYpa_comp)
RFTMOYpa_comp <- tax_table(RFTMOYpa_comp)
rownames(RFTMOYpa_comp)
sample_names(SAMP) 
sample_names(OTU)
sample_names(RFTMOYpa_comp) 
phyloRFTMOY_PA <- phyloseq(OTU, RFTMOYpa_comp, SAMP) 
phyloRFTMOY_PA


##RFTM_pa1.Site.xSW
RFTMSWpa_comp <- read.csv("Reduced Tax Tables/PresAbs_RFTMSW_COMP.csv") 
RFTMSWpa_comp <- as.matrix(RFTMSWpa_comp) 
colnames(RFTMSWpa_comp) <- c("X", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus.x","Genus.y","Species")
rownames(RFTMSWpa_comp)
RFTMSWpa_comp <- tax_table(RFTMSWpa_comp)
rownames(RFTMSWpa_comp)
sample_names(SAMP) 
sample_names(OTU)
sample_names(RFTMSWpa_comp) 
phyloRFTMSW_PA <- phyloseq(OTU, RFTMSWpa_comp, SAMP) 
phyloRFTMSW_PA


##peacrabs.f1.Site.xOY
PEAOYpa_comp <- read.csv("Reduced Tax Tables/PresAbs_PEAOY_COMP.csv") 
PEAOYpa_comp <- as.matrix(PEAOYpa_comp) 
colnames(PEAOYpa_comp) <- c("X", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus.x","Genus.y","Species")
rownames(PEAOYpa_comp)
PEAOYpa_comp <- tax_table(PEAOYpa_comp)
rownames(PEAOYpa_comp)
sample_names(SAMP) 
sample_names(OTU)
sample_names(PEAOYpa_comp) 
phyloPEAOY_PA <- phyloseq(OTU, PEAOYpa_comp, SAMP) 
phyloPEAOY_PA


##peacrabs.f1.Site.xSW
PEASWpa_comp <- read.csv("Reduced Tax Tables/PresAbs_PEASW_COMP.csv") 
PEASWpa_comp <- as.matrix(PEASWpa_comp) 
colnames(PEASWpa_comp) <- c("X", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus.x","Genus.y","Species")
rownames(PEASWpa_comp)
PEASWpa_comp <- tax_table(PEASWpa_comp)
rownames(PEASWpa_comp)
sample_names(SAMP) 
sample_names(OTU)
sample_names(PEASWpa_comp) 
phyloPEASW_PA <- phyloseq(OTU, PEASWpa_comp, SAMP) 
phyloPEASW_PA


# PHYLOSEQ OBJECTS ####
#phyloPEA_PA
#phyloRFTMPEA_PA
#phyloRFTM_PA
#phyloSITEOYNW_PA
#phyloSITESWNW_PA
#phyloRFTMOY_PA
#phyloRFTMSW_PA
#phyloPEAOY_PA
#phyloPEASW_PA


# TAXONOMY TABLE ID ####
#PEApa_comp
#RFTMPEApa_comp
#RFTMpa_comp
#SITEOYNWpa_comp
#SITESWNWpa_comp
#RFTMOYpa_comp
#RFTMSWpa_comp
#PEAOYpa_comp
#PEASWpa_comp



##### Figuring out the number of class groupings for each ID ####

#PEApa_comp
PEApa_tax <- read.csv("Reduced Tax Tables/PresAbs_PEA_COMP.csv") 
PEApa_tax_tab <- table(PEApa_tax$Class, useNA = "ifany")
view(PEApa_tax_tab)

PEApa_tax <- read.csv("Reduced Tax Tables/PresAbs_PEA_COMP.csv") 
PEApa_tax_ord <- table(PEApa_tax$Order, useNA = "ifany")
view(PEApa_tax_ord)



#RFTMPEApa_comp
RFTMPEApa_tax <- read.csv("Reduced Tax Tables/PresAbs_RFTMPEA_COMP.csv") 
RFTMPEApa_tax_tab <- table(RFTMPEApa_tax$Class, useNA = "ifany")
view(RFTMPEApa_tax_tab)

#RFTMpa_comp
RFTMpa_tax <- read.csv("Reduced Tax Tables/PresAbs_RFTM_COMP.csv") 
RFTMpa_tax_tab <- table(RFTMpa_tax$Class, useNA = "ifany")
view(RFTMpa_tax_tab)


#SITEOYNWpa_comp
SITEOYNWpa_tax <- read.csv("Reduced Tax Tables/PresAbs_SITEOYNW_COMP.csv") 
SITEOYNWpa_tax_tab <- table(SITEOYNWpa_tax$Class, useNA = "ifany")
view(SITEOYNWpa_tax_tab)

#SITESWNWpa_comp
SITESWNWpa_tax <- read.csv("Reduced Tax Tables/PresAbs_SITESWNW_COMP.csv") 
SITESWNWpa_tax_tab <- table(SITESWNWpa_tax$Class, useNA = "ifany")
view(SITESWNWpa_tax_tab)

#RFTMOYpa_comp
RFTMOYpa_tax <- read.csv("Reduced Tax Tables/PresAbs_RFTMOY_COMP.csv") 
RFTMOYpa_tax_tab <- table(RFTMOYpa_tax$Class, useNA = "ifany")
view(RFTMOYpa_tax_tab)

#RFTMSWpa_comp
RFTMSWpa_tax <- read.csv("Reduced Tax Tables/PresAbs_RFTMSW_COMP.csv") 
RFTMSWpa_tax_tab <- table(RFTMSWpa_tax$Class, useNA = "ifany")
view(RFTMSWpa_tax_tab)

#PEAOYpa_comp
PEAOYpa_tax <- read.csv("Reduced Tax Tables/PresAbs_PEAOY_COMP.csv") 
PEAOYpa_tax_tab <- table(PEAOYpa_tax$Class, useNA = "ifany")
view(PEAOYpa_tax_tab)


#PEASWpa_comp
PEASWpa_tax <- read.csv("Reduced Tax Tables/PresAbs_PEASW_COMP.csv") 
PEASWpa_tax_tab <- table(PEASWpa_tax$Class, useNA = "ifany")
view(PEASWpa_tax_tab)


##### Abundance Plot for RFTM_PA ####

#phyloPEA_PA
plot_bar(phyloPEA_PA, fill="Class")+
  theme(legend.position="right", legend.text=element_text(size=8),
        axis.ticks.x = element_blank(), 
        axis.line = element_line(color="black"), text = element_text(size=10), 
        plot.title = element_text(family = "sans", face = "bold", hjust = 0.5, size = 15, colour = "black", margin=margin(15, 15, 0, 15)), 
        plot.subtitle = element_text(family = "sans", face = "bold", hjust = 0.5, color = "black", margin=margin(15, 15, 15, 15), ),
        plot.caption = element_text(color = "black", face = "italic", margin=margin(15, 15, 15, 15)),
        axis.title.x = element_text(family = "sans", face = "bold", size = 15, margin=margin(15, 15, 15, 15)), 
        axis.title.y = element_text(family = "sans", face = "bold", size = 15, margin=margin(15, 15, 15, 15)))+
  labs(title = "Microbial Abundance for Pea Crab Prescence Absence ", 
       subtitle = "Taxonomic Rank: Class",
       caption = "Data source: Oyster 16s 2017",
       x = "Sample ID",
       y = "Abundance")
ggsave(filename = "Abun_PEApa_class.jpeg", plot=last_plot(), path ="PEA_RFTM/G_RFTM_pa/", width = 15, height = 10)


#phyloRFTMPEA_PA
plot_bar(phyloRFTMPEA_PA, fill="Class")+
  theme(legend.position="right", legend.text=element_text(size=8),
        axis.ticks.x = element_blank(), 
        axis.line = element_line(color="black"), text = element_text(size=10), 
        plot.title = element_text(family = "sans", face = "bold", hjust = 0.5, size = 15, colour = "black", margin=margin(15, 15, 0, 15)), 
        plot.subtitle = element_text(family = "sans", face = "bold", hjust = 0.5, color = "black", margin=margin(15, 15, 15, 15), ),
        plot.caption = element_text(color = "black", face = "italic", margin=margin(15, 15, 15, 15)),
        axis.title.x = element_text(family = "sans", face = "bold", size = 15, margin=margin(15, 15, 15, 15)), 
        axis.title.y = element_text(family = "sans", face = "bold", size = 15, margin=margin(15, 15, 15, 15)))+
  labs(title = "Microbial Abundance for Pea Crab and RFTM Prescence Absence ", 
       subtitle = "Taxonomic Rank: Class",
       caption = "Data source: Oyster 16s 2017",
       x = "Sample ID",
       y = "Abundance")
ggsave(filename = "Abun_RFTMPEApa_class.jpeg", plot=last_plot(), path ="PEA_RFTM/G_RFTM_pa/", width = 15, height = 10)


#phyloRFTM_PA
plot_bar(phyloRFTM_PA, fill="Class")+
  theme(legend.position="right", legend.text=element_text(size=8),
        axis.ticks.x = element_blank(), 
        axis.line = element_line(color="black"), text = element_text(size=10), 
        plot.title = element_text(family = "sans", face = "bold", hjust = 0.5, size = 15, colour = "black", margin=margin(15, 15, 0, 15)), 
        plot.subtitle = element_text(family = "sans", face = "bold", hjust = 0.5, color = "black", margin=margin(15, 15, 15, 15), ),
        plot.caption = element_text(color = "black", face = "italic", margin=margin(15, 15, 15, 15)),
        axis.title.x = element_text(family = "sans", face = "bold", size = 15, margin=margin(15, 15, 15, 15)), 
        axis.title.y = element_text(family = "sans", face = "bold", size = 15, margin=margin(15, 15, 15, 15)))+
  labs(title = "Microbial Abundance for RFTM Prescence Absence ", 
       subtitle = "Taxonomic Rank: Class",
       caption = "Data source: Oyster 16s 2017",
       x = "Sample ID",
       y = "Abundance")
ggsave(filename = "Abun_RFTMpa_class.jpeg", plot=last_plot(), path ="PEA_RFTM/G_RFTM_pa/", width = 15, height = 10)


#phyloSITEOYNW_PA
plot_bar(phyloSITEOYNW_PA, fill="Class")+
  theme(legend.position="right", legend.text=element_text(size=8),
        axis.ticks.x = element_blank(), 
        axis.line = element_line(color="black"), text = element_text(size=10), 
        plot.title = element_text(family = "sans", face = "bold", hjust = 0.5, size = 15, colour = "black", margin=margin(15, 15, 0, 15)), 
        plot.subtitle = element_text(family = "sans", face = "bold", hjust = 0.5, color = "black", margin=margin(15, 15, 15, 15), ),
        plot.caption = element_text(color = "black", face = "italic", margin=margin(15, 15, 15, 15)),
        axis.title.x = element_text(family = "sans", face = "bold", size = 15, margin=margin(15, 15, 15, 15)), 
        axis.title.y = element_text(family = "sans", face = "bold", size = 15, margin=margin(15, 15, 15, 15)))+
  labs(title = "Microbial Abundance for Sites OY and NW Prescence Absence ", 
       subtitle = "Taxonomic Rank: Class",
       caption = "Data source: Oyster 16s 2017",
       x = "Sample ID",
       y = "Abundance")
ggsave(filename = "Abun_SITEOYNWpa_class.jpeg", plot=last_plot(), path ="PEA_RFTM/G_RFTM_pa/", width = 15, height = 10)


#phyloSITESWNW_PA
plot_bar(phyloSITESWNW_PA, fill="Class")+
  theme(legend.position="right", legend.text=element_text(size=8),
        axis.ticks.x = element_blank(), 
        axis.line = element_line(color="black"), text = element_text(size=10), 
        plot.title = element_text(family = "sans", face = "bold", hjust = 0.5, size = 15, colour = "black", margin=margin(15, 15, 0, 15)), 
        plot.subtitle = element_text(family = "sans", face = "bold", hjust = 0.5, color = "black", margin=margin(15, 15, 15, 15), ),
        plot.caption = element_text(color = "black", face = "italic", margin=margin(15, 15, 15, 15)),
        axis.title.x = element_text(family = "sans", face = "bold", size = 15, margin=margin(15, 15, 15, 15)), 
        axis.title.y = element_text(family = "sans", face = "bold", size = 15, margin=margin(15, 15, 15, 15)))+
  labs(title = "Microbial Abundance for Sites SW and NW Prescence Absence ", 
       subtitle = "Taxonomic Rank: Class",
       caption = "Data source: Oyster 16s 2017",
       x = "Sample ID",
       y = "Abundance")
ggsave(filename = "Abun_SITESWNWpa_class.jpeg", plot=last_plot(), path ="PEA_RFTM/G_RFTM_pa/", width = 15, height = 10)


#phyloRFTMOY_PA
plot_bar(phyloRFTMOY_PA, fill="Class")+
  theme(legend.position="right", legend.text=element_text(size=8),
        axis.ticks.x = element_blank(), 
        axis.line = element_line(color="black"), text = element_text(size=10), 
        plot.title = element_text(family = "sans", face = "bold", hjust = 0.5, size = 15, colour = "black", margin=margin(15, 15, 0, 15)), 
        plot.subtitle = element_text(family = "sans", face = "bold", hjust = 0.5, color = "black", margin=margin(15, 15, 15, 15), ),
        plot.caption = element_text(color = "black", face = "italic", margin=margin(15, 15, 15, 15)),
        axis.title.x = element_text(family = "sans", face = "bold", size = 15, margin=margin(15, 15, 15, 15)), 
        axis.title.y = element_text(family = "sans", face = "bold", size = 15, margin=margin(15, 15, 15, 15)))+
  labs(title = "Microbial Abundance for RFTM and Site OY Prescence Absence ", 
       subtitle = "Taxonomic Rank: Class",
       caption = "Data source: Oyster 16s 2017",
       x = "Sample ID",
       y = "Abundance")
ggsave(filename = "Abun_RFTMOYpa_class.jpeg", plot=last_plot(), path ="PEA_RFTM/G_RFTM_pa/", width = 15, height = 10)


#phyloRFTMSW_PA
plot_bar(phyloRFTMSW_PA, fill="Class")+
  theme(legend.position="right", legend.text=element_text(size=8),
        axis.ticks.x = element_blank(), 
        axis.line = element_line(color="black"), text = element_text(size=10), 
        plot.title = element_text(family = "sans", face = "bold", hjust = 0.5, size = 15, colour = "black", margin=margin(15, 15, 0, 15)), 
        plot.subtitle = element_text(family = "sans", face = "bold", hjust = 0.5, color = "black", margin=margin(15, 15, 15, 15), ),
        plot.caption = element_text(color = "black", face = "italic", margin=margin(15, 15, 15, 15)),
        axis.title.x = element_text(family = "sans", face = "bold", size = 15, margin=margin(15, 15, 15, 15)), 
        axis.title.y = element_text(family = "sans", face = "bold", size = 15, margin=margin(15, 15, 15, 15)))+
  labs(title = "Microbial Abundance for RFTM and Site SW Prescence Absence ", 
       subtitle = "Taxonomic Rank: Class",
       caption = "Data source: Oyster 16s 2017",
       x = "Sample ID",
       y = "Abundance")
ggsave(filename = "Abun_RFTMSWpa_class.jpeg", plot=last_plot(), path ="PEA_RFTM/G_RFTM_pa/", width = 15, height = 10)


#phyloPEAOY_PA
plot_bar(phyloPEAOY_PA, fill="Class")+
  theme(legend.position="right", legend.text=element_text(size=8),
        axis.ticks.x = element_blank(), 
        axis.line = element_line(color="black"), text = element_text(size=10), 
        plot.title = element_text(family = "sans", face = "bold", hjust = 0.5, size = 15, colour = "black", margin=margin(15, 15, 0, 15)), 
        plot.subtitle = element_text(family = "sans", face = "bold", hjust = 0.5, color = "black", margin=margin(15, 15, 15, 15), ),
        plot.caption = element_text(color = "black", face = "italic", margin=margin(15, 15, 15, 15)),
        axis.title.x = element_text(family = "sans", face = "bold", size = 15, margin=margin(15, 15, 15, 15)), 
        axis.title.y = element_text(family = "sans", face = "bold", size = 15, margin=margin(15, 15, 15, 15)))+
  labs(title = "Microbial Abundance for Pea Crabs and Site OY Prescence Absence ", 
       subtitle = "Taxonomic Rank: Class",
       caption = "Data source: Oyster 16s 2017",
       x = "Sample ID",
       y = "Abundance")
ggsave(filename = "Abun_PEAOYpa_class.jpeg", plot=last_plot(), path ="PEA_RFTM/G_RFTM_pa/", width = 15, height = 10)


#phyloPEASW_PA
plot_bar(phyloPEASW_PA, fill="Class")+
  theme(legend.position="right", legend.text=element_text(size=8),
        axis.ticks.x = element_blank(), 
        axis.line = element_line(color="black"), text = element_text(size=10), 
        plot.title = element_text(family = "sans", face = "bold", hjust = 0.5, size = 15, colour = "black", margin=margin(15, 15, 0, 15)), 
        plot.subtitle = element_text(family = "sans", face = "bold", hjust = 0.5, color = "black", margin=margin(15, 15, 15, 15), ),
        plot.caption = element_text(color = "black", face = "italic", margin=margin(15, 15, 15, 15)),
        axis.title.x = element_text(family = "sans", face = "bold", size = 15, margin=margin(15, 15, 15, 15)), 
        axis.title.y = element_text(family = "sans", face = "bold", size = 15, margin=margin(15, 15, 15, 15)))+
  labs(title = "Microbial Abundance for Pea Crabs and Site SW Prescence Absence ", 
       subtitle = "Taxonomic Rank: Class",
       caption = "Data source: Oyster 16s 2017",
       x = "Sample ID",
       y = "Abundance")
ggsave(filename = "Abun_PEASWpa_class.jpeg", plot=last_plot(), path ="PEA_RFTM/G_RFTM_pa/", width = 15, height = 10)






