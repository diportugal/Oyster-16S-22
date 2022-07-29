#Adding in the Log2FoldChange Column


#Testing with the main object
physeqlog <- physeq

physeqlog_dds <- phyloseq_to_deseq2(physeqlog, ~RFTM_pa) 
physeqlog_dds <- DESeq(physeqlog_dds, test="Wald", fitType="parametric")

physeqlog_res= results(physeqlog_dds, cooksCutoff = FALSE)
alpha = 0.05
physeqlog_sig = physeqlog_res[which(physeqlog_res$padj < alpha), ]
physeqlog_sig = cbind(as(physeqlog_sig, "data.frame"), as(tax_table(physeqlog)[rownames(physeqlog_sig), ], "matrix"))

#Simplifying the columns and adding 
view(physeqlog_sig)
physeqlog_sig <- select(physeqlog_sig, -baseMean,-lfcSE, -stat, -pvalue, -padj, -Genus.y)
names(physeqlog_sig)[names(physeqlog_sig) == "Genus.x"] <- "Genus"
view(physeqlog_sig)



---------------------------------------------------------------------------------------------------

ddsPA_rftmpeaM <- phyloseq_to_deseq2(physeq, ~ RFTM_pa * peacrabs.f)
ddsPA_rftmpeaM <- DESeq(ddsPA_rftmpeaM, test="Wald", fitType="parametric")
resultsNames(ddsPA_rftmpeaM)

###RFTM_pa_1_vs_0####
resPA_rftm = results(ddsPA_rftmpeaM, cooksCutoff = FALSE, name="RFTM_pa_1_vs_0")
alpha = 0.05
sigPA_rftm = resPA_rftm[which(resPA_rftm$padj < alpha), ]
sigPA_rftm = cbind(as(sigPA_rftm, "data.frame"), as(tax_table(physeq)[rownames(sigPA_rftm), ], "matrix"))
view(sigPA_rftm)
sigPA_rftm <- select(sigPA_rftm, -baseMean, -lfcSE, -X, -stat, -pvalue, -padj, -Genus.y)
names(sigPA_rftm)[names(sigPA_rftm) == "Genus.x"] <- "Genus"
view(sigPA_rftm)
dim(sigPA_rftm)
#458 8
write.table(sigPA_rftm, file="Log2Fold/RFTMpa_compLog.csv", quote=FALSE,sep = ",", col.names=T)



###peacrabs.f_1_vs_0####
resPA_pea = results(ddsPA_rftmpeaM, cooksCutoff = FALSE, name="peacrabs.f_1_vs_0")
alpha = 0.05
sigPA_pea = resPA_pea[which(resPA_pea$padj < alpha), ]
sigPA_pea = cbind(as(sigPA_pea, "data.frame"), as(tax_table(physeq)[rownames(sigPA_pea), ], "matrix"))
view(sigPA_pea)
sigPA_pea <- select(sigPA_pea, -baseMean, -lfcSE, -X, -stat, -pvalue, -padj, -Genus.y)
names(sigPA_pea)[names(sigPA_pea) == "Genus.x"] <- "Genus"
view(sigPA_pea)
dim(sigPA_pea)
#265 8 


###RFTM_pa1.peacrabs.f1####
resPA_rftmpea = results(ddsPA_rftmpeaM, cooksCutoff = FALSE, name="RFTM_pa1.peacrabs.f1")
alpha = 0.05
sigPA_rftmpea = resPA_rftmpea[which(resPA_rftmpea$padj < alpha), ]
sigPA_rftmpea = cbind(as(sigPA_rftmpea, "data.frame"), as(tax_table(physeq)[rownames(sigPA_rftmpea), ], "matrix"))
view(sigPA_rftmpea)
sigPA_rftmpea <- select(sigPA_rftmpea, -baseMean, -lfcSE, -X, -stat, -pvalue, -padj, -Genus.y)
names(sigPA_rftmpea)[names(sigPA_rftmpea) == "Genus.x"] <- "Genus"
view(sigPA_rftmpea)
dim(sigPA_rftmpea)
#35 8

-------------------------------------------------------------------------------------------------------------------------

ddsPA_rftmsiteM <- phyloseq_to_deseq2(physeq, ~ RFTM_pa * Site.x)
ddsPA_rftmsiteM <- DESeq(ddsPA_rftmsiteM, test="Wald", fitType="parametric")
resultsNames(ddsPA_rftmsiteM)


###Site.x_OY_vs_NW####
resPA_siteOYNW = results(ddsPA_rftmsiteM, cooksCutoff = FALSE, name="Site.x_OY_vs_NW")
alpha = 0.05
sigPA_siteOYNW = resPA_siteOYNW[which(resPA_siteOYNW$padj < alpha), ]
sigPA_siteOYNW = cbind(as(sigPA_siteOYNW, "data.frame"), as(tax_table(physeq)[rownames(sigPA_siteOYNW), ], "matrix"))
view(sigPA_siteOYNW)
sigPA_siteOYNW <- select(sigPA_siteOYNW, -baseMean, -lfcSE, -X, -stat, -pvalue, -padj, -Genus.y)
names(sigPA_siteOYNW)[names(sigPA_siteOYNW) == "Genus.x"] <- "Genus"
view(sigPA_siteOYNW)
dim(sigPA_siteOYNW)
#313 8


###sigPA_siteSWNW####
resPA_siteSWNW = results(ddsPA_rftmsiteM, cooksCutoff = FALSE, name="Site.x_SW_vs_NW")
alpha = 0.05
sigPA_siteSWNW = resPA_siteSWNW[which(resPA_siteSWNW$padj < alpha), ]
sigPA_siteSWNW = cbind(as(sigPA_siteSWNW, "data.frame"), as(tax_table(physeq)[rownames(sigPA_siteSWNW), ], "matrix"))
view(sigPA_siteSWNW)
sigPA_siteSWNW <- select(sigPA_siteSWNW, -baseMean, -lfcSE, -X, -stat, -pvalue, -padj, -Genus.y)
names(sigPA_siteSWNW)[names(sigPA_siteSWNW) == "Genus.x"] <- "Genus"
view(sigPA_siteSWNW)
dim(sigPA_siteSWNW)
#295 8


###RFTM_pa1.Site.xOY####
resPA_RFTMOY = results(ddsPA_rftmsiteM, cooksCutoff = FALSE, name="RFTM_pa1.Site.xOY")
alpha = 0.05
sigPA_RFTMOY = resPA_RFTMOY[which(resPA_RFTMOY$padj < alpha), ]
sigPA_RFTMOY = cbind(as(sigPA_RFTMOY, "data.frame"), as(tax_table(physeq)[rownames(sigPA_RFTMOY), ], "matrix"))
view(sigPA_RFTMOY)
sigPA_RFTMOY <- select(sigPA_RFTMOY, -baseMean, -lfcSE, -X, -stat, -pvalue, -padj, -Genus.y)
names(sigPA_RFTMOY)[names(sigPA_RFTMOY) == "Genus.x"] <- "Genus"
view(sigPA_RFTMOY)
dim(sigPA_RFTMOY)
#184 8


###RFTM_pa1.Site.xSW####
resPA_RFTMSW = results(ddsPA_rftmsiteM, cooksCutoff = FALSE, name="RFTM_pa1.Site.xSW")
alpha = 0.05
sigPA_RFTMSW = resPA_RFTMSW[which(resPA_RFTMSW$padj < alpha), ]
sigPA_RFTMSW = cbind(as(sigPA_RFTMSW, "data.frame"), as(tax_table(physeq)[rownames(sigPA_RFTMSW), ], "matrix"))
view(sigPA_RFTMSW)
sigPA_RFTMSW <- select(sigPA_RFTMSW, -baseMean, -lfcSE, -X, -stat, -pvalue, -padj, -Genus.y)
names(sigPA_RFTMSW)[names(sigPA_RFTMSW) == "Genus.x"] <- "Genus"
view(sigPA_RFTMSW)
dim(sigPA_RFTMSW)
#209 8


-------------------------------------------------------------------------------------------------------------------------
  
ddsPA_PeasiteM <- phyloseq_to_deseq2(physeq, ~ peacrabs.f * Site.x)
ddsPA_PeasiteM <- DESeq(ddsPA_PeasiteM, test="Wald", fitType="parametric")
resultsNames(ddsPA_PeasiteM)

###peacrabs.f1.Site.xOY####
resPA_PeaOY = results(ddsPA_PeasiteM, cooksCutoff = FALSE, name="peacrabs.f1.Site.xOY")
alpha = 0.05
sigPA_PeaOY = resPA_PeaOY[which(resPA_PeaOY$padj < alpha), ]
sigPA_PeaOY = cbind(as(sigPA_PeaOY, "data.frame"), as(tax_table(physeq)[rownames(sigPA_PeaOY), ], "matrix"))
view(sigPA_PeaOY)
sigPA_PeaOY <- select(sigPA_PeaOY, -baseMean, -lfcSE, -X, -stat, -pvalue, -padj, -Genus.y)
names(sigPA_PeaOY)[names(sigPA_PeaOY) == "Genus.x"] <- "Genus"
view(sigPA_PeaOY)
dim(sigPA_PeaOY)
#50 8

###peacrabs.f1.Site.xSW####
resPA_PeaSW = results(ddsPA_PeasiteM, cooksCutoff = FALSE, name="peacrabs.f1.Site.xSW")
alpha = 0.05
sigPA_PeaSW = resPA_PeaSW[which(resPA_PeaSW$padj < alpha), ]
sigPA_PeaSW = cbind(as(sigPA_PeaSW, "data.frame"), as(tax_table(physeq)[rownames(sigPA_PeaSW), ], "matrix"))
view(sigPA_PeaSW)
sigPA_PeaSW <- select(sigPA_PeaSW, -baseMean, -lfcSE, -X, -stat, -pvalue, -padj, -Genus.y)
names(sigPA_PeaSW)[names(sigPA_PeaSW) == "Genus.x"] <- "Genus"
view(sigPA_PeaSW)
dim(sigPA_PeaSW)
#37 8



