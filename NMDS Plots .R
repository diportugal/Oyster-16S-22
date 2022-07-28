# Other plot types 


# PHYLOSEQ OBJECTS ####
## TAX TABS ####

#PEApa_comp
##PEApa_tax_tab

#RFTMPEApa_comp
##RFTMPEApa_tax_tab

#RFTMpa_comp
##RFTMpa_tax_tab

#SITEOYNWpa_comp
##SITEOYNWpa_tax_tab

#SITESWNWpa_comp
##SITESWNWpa_tax_tab

#RFTMOYpa_comp
##RFTMOYpa_tax_tab

#RFTMSWpa_comp
##RFTMSWpa_tax_tab

#PEAOYpa_comp
##PEAOYpa_tax_tab

#PEASWpa_comp
##PEASWpa_tax_tab


#######################################
physeq

physeq = physeq
physeqwh0 = genefilter_sample(physeq, filterfun_sample(function(x) x > 5), A=0.5*nsamples(physeq))
physeq1 = prune_taxa(physeqwh0, physeq)
physeq1 = transform_sample_counts(physeq1, function(x) 1E6 * x/sum(x))

class.sum = tapply(taxa_sums(physeq), tax_table(physeq)[, "Class"], sum, na.rm=TRUE)
top5class_nmds = names(sort(class.sum, TRUE))[1:5]
physeq_class_5class = prune_taxa((tax_table(physeq)[, "Class"] %in% top5class_nmds), physeq)
physeq_class5_no_na <- subset_samples(physeq_class_5class, "Class" != "NA", )
data17.ord_class5_no_na <- ordinate(physeq_class5_no_na, "NMDS", "bray")

plot_ordination(physeq_class_5class, data17.ord_class5_no_na, type="taxa", color="Class")


############################################























