# Other plot types 


#PHYLOSEQ OBJECTS####
#phyloPEA_PA
#phyloRFTMPEA_PA
#phyloRFTM_PA
#phyloSITEOYNW_PA
#phyloSITESWNW_PA
#phyloRFTMOY_PA
#phyloRFTMSW_PA
#phyloPEAOY_PA
#phyloPEASW_PA


##### Plot Richness for RFTM_PA####
plot_richness(phyloRFTM_PA, x="RFTM_pa", measures = c("Observed", "Chao", "Simpson", "Shannon"))+
  geom_boxplot()+
  stat_compare_means(method = "t.test", label.y = 0)+
  theme(legend.position="right", legend.text=element_text(size=8),
        axis.ticks.x = element_blank(), 
        axis.line = element_line(color="black"), text = element_text(size=10), 
        plot.title = element_text(family = "sans", face = "bold", hjust = 0.5, size = 15, colour = "black", margin=margin(15, 15, 0, 15)), 
        plot.subtitle = element_text(family = "sans", face = "bold", hjust = 0.5, color = "black", margin=margin(15, 15, 15, 15), ),
        plot.caption = element_text(color = "black", face = "italic", margin=margin(15, 15, 15, 15)),
        axis.title.x = element_text(family = "sans", size = 15, margin=margin(15, 15, 15, 15)), 
        axis.title.y = element_text(family = "sans", size = 15, margin=margin(15, 15, 15, 15)))+
  labs(title = "Microbial Richness ", 
       subtitle = "RFTM Prescence Absence",
       caption = "Data source: Oyster 16s 2017",
       x = "RFTM Prescence Absence",
       y = "Richness")
ggsave(filename = "Rich_RFTMpa.jpeg", plot=last_plot(), path ="PEA_RFTM/G_RFTM_pa/", width = 15, height = 10)






























