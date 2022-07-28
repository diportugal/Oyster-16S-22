#Adding in the Log2FoldChange Column

physeq_t_log <- physeq







view(SAMP) 
view(OTU) 
view(TAX) 

#Start here ####

#Analysis of RFTM Scores ####
#Make RFTM_score.x not a factor - Make it a numerical value 
RFTM_dds17 <- phyloseq_to_deseq2(physeq_class, ~RFTM_score.asnum) 
RFTM_dds17 <- DESeq(RFTM_dds17, test="Wald", fitType="parametric")

RFTM_res17= results(RFTM_dds17, cooksCutoff = FALSE)
alpha = 0.05
RFTM_sig17 = RFTM_res17[which(RFTM_res17$padj < alpha), ]
RFTM_sig17 = cbind(as(RFTM_sig17, "data.frame"), as(tax_table(physeq_class)[rownames(RFTM_sig17), ], "matrix"))

head(RFTM_sig17)
dim(RFTM_sig17)