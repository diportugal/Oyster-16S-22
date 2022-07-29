#Full Script - Modified 6/28/22

#Libraries ####
#install.packages("DESeq2")
#install.packages("genefilter")
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("plyr"); packageVersion("plyr")
library("data.table")
library("DESeq2")
library("GenomicRanges")
library("GenomeInfoDb") 
library("dplyr")
library("tidyverse")
require("genefilter")
theme_set(theme_bw())


#creating the PHYLOSEQ object ####

c_meta17data <- read.csv("Oyster_data_raw/cleanmetadata17")
  
asvtable17 <- fread("Oyster_data_raw/asvtable_de17.csv")

run23 <- read.csv("Oyster_data_raw/Run123_taxa_complete.csv")

#***************************************************************************************************************************************

## CHANGING ROW NAMES FOR EACH DATA SET 
rownames(c_meta17data) = c_meta17data$X
c_meta17data$X=NULL
rownames(c_meta17data)
view(c_meta17data)
#ROW NAMES ARE THE UNIQUE IDs 

rownames(asvtable17) = asvtable17$V1
rownames(asvtable17)
asvtable17$V1= NULL
#ROW NAMES ARE THE UNIQUE IDs 

rownames(run23) = run23$Row.names
run23$Row.names = NULL        
rownames(run23)
#ROW NAMES ARE THE SEQUENCE 


## CONVERTING TO MATRICIES 
meta17_matrix <- as.data.frame(c_meta17data, rownames("X"))
rownames(meta17_matrix)
#STILL UNIQUE ID

otumat_matrix <- as.matrix(asvtable17, rownames=rownames(asvtable17))
rownames(otumat_matrix)
#STILL UNIQUE ID 


rownames(run23)

taxmat_matrix <- as.matrix(run23) 
colnames(taxmat_matrix) <- c("X", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus.x","Genus.y","Species")
rownames(taxmat_matrix)
#STILL SEQUENCE 


## SETTING OTU, TAX, SAMP 
OTU <- otu_table(otumat_matrix, taxa_are_rows = FALSE)
rownames(OTU) #UniqueID

TAX <- tax_table(taxmat_matrix)
rownames(TAX) #Sequence

SAMP <- sample_data(meta17_matrix)
rownames(SAMP)#UniqueID


## INSPECTING SAMPLE NAMES
sample_names(SAMP) #UniqueID
sample_names(OTU) #UniqueID
sample_names(TAX) #NULL


## EVENING OUT THE DATA
OTU=transform_sample_counts(OTU, function(x) 1E6 * x/sum(x))

## READING THROUGH PHYLOSEQ
physeq = phyloseq(OTU, TAX, SAMP) 
physeq

#Final Object (Original)
physeq


#***************************************************************************************************************************************

#Make peacrabs into a factor
SAMP$peacrabs.f <- factor(SAMP$peacrabs.x)
typeof(SAMP$peacrabs.f)#Remains "integer"
is.factor(SAMP$peacrabs.f) #True

#Make RFTM_score into a factor
SAMP$RFTM_score.f <- factor(SAMP$RFTM_score.x)
typeof(SAMP$RFTM_score.f) #Remains "integer"
is.factor(SAMP$RFTM_score.f)

typeof(SAMP$peacrabs.x) #Integer
typeof(SAMP$RFTM_score.x) #Double


#make a simpler rftm factor without 0.5 and combining 4 and 5
#changing 0.5 to 1 and 5 to 4, we are left with 1,2,3,4
SAMP$RFTM_simp <- SAMP$RFTM_score.f
SAMP$RFTM_simp <- sub("0.5", "1", SAMP$RFTM_simp)
SAMP$RFTM_simp <- sub("5", "4", SAMP$RFTM_simp)


#make RFTM presence absence
#Anything equal to O is absence, everything else=(1,2,3,4) is presence
SAMP$RFTM_pa <- ifelse(SAMP$RFTM_score.x=="0", 0, 1)
SAMP$RFTM_pa <- factor(SAMP$RFTM_pa)


physeq = phyloseq(OTU, TAX, SAMP) 

view(SAMP)

#***************************************************************************************************************************************

## Intersection of RFTM and Pea */+ ####
ddsPA_rftmpeaP <- phyloseq_to_deseq2(physeq, ~ RFTM_pa + peacrabs.f)
ddsPA_rftmpeaP <- DESeq(ddsPA_rftmpeaP, test="Wald", fitType="parametric")
resultsNames(ddsPA_rftmpeaP)

#Result Names
#RFTM_pa_1_vs_0
#peacrabs.f_1_vs_0


ddsPA_rftmpeaM <- phyloseq_to_deseq2(physeq, ~ RFTM_pa * peacrabs.f)
ddsPA_rftmpeaM <- DESeq(ddsPA_rftmpeaM, test="Wald", fitType="parametric")
resultsNames(ddsPA_rftmpeaM)

#Result Names
#RFTM_pa_1_vs_0
#peacrabs.f_1_vs_0
#RFTM_pa1.peacrabs.f1"



###RFTM_pa_1_vs_0####
resPA_rftm <- results(ddsPA_rftmpeaM, name="RFTM_pa_1_vs_0")
sigPA_rftm <- resPA_rftm[which(resPA_rftm$padj < 0.05), ]
dim(sigPA_rftm)
view(sigPA_rftm)
#148 6

st_sigPA_rftm <- subset_taxa(prune_taxa(rownames(sigPA_rftm), physeq))
st_sigPA_rftm

seq_sigPA_rftm <- as.data.frame(tax_table(st_sigPA_rftm))

#Creating Reduced Tax Tables for this variable
write.table(seq_sigPA_rftm, file="Reduced Tax Tables/PresAbs_RFTM_COMP.csv", quote=FALSE,sep = ",", col.names=T)
            

#Getting the +/- Log of this variable

#POSITIVE
sigPA_rftm_pos <- sigPA_rftm[sigPA_rftm$log2FoldChange>0,]
dim(sigPA_rftm_pos)
# 77 6

st_sigPA_rftm_pos <- subset_taxa(prune_taxa(rownames(sigPA_rftm_pos), physeq))
st_sigPA_rftm_pos

seq_sigPA_rftm_pos <- as.data.frame(tax_table(st_sigPA_rftm_pos))

#Creating Reduced Tax Tables for this variable
write.table(seq_sigPA_rftm_pos, file="Reduced Tax Tables/PresAbs_RFTM_POS.csv", quote=FALSE,sep = ",", col.names=T)


#NEGATIVE
sigPA_rftm_neg <- sigPA_rftm[sigPA_rftm$log2FoldChange<0,]
dim(sigPA_rftm_neg)
# 71 6

st_sigPA_rftm_neg <- subset_taxa(prune_taxa(rownames(sigPA_rftm_neg), physeq))
st_sigPA_rftm_neg

seq_sigPA_rftm_neg <- as.data.frame(tax_table(st_sigPA_rftm_neg))

write.table(seq_sigPA_rftm_neg, file="Reduced Tax Tables/PresAbs_RFTM_NEG.csv", quote=FALSE,sep = ",", col.names=T)


###peacrabs.f_1_vs_0####
resPA_pea <- results(ddsPA_rftmpeaM, name="peacrabs.f_1_vs_0")
sigPA_pea <- resPA_pea[which(resPA_pea$padj < 0.05), ]
dim(sigPA_pea)
# 66 6

st_sigPA_pea <- subset_taxa(prune_taxa(rownames(sigPA_pea), physeq))
st_sigPA_pea

seq_sigPA_pea <- as.data.frame(tax_table(st_sigPA_pea))

#Creating Reduced Tax Tables for this variable
write.table(seq_sigPA_pea, file="Reduced Tax Tables/PresAbs_PEA_COMP.csv", quote=FALSE,sep = ",", col.names=T)


#Getting the +/- Log of this variable

#POSITIVE - DOES NOT WORK BECAUSE THERE ARE NO POSITIVE VALUES

#sigPA_pea
#sigPA_pea_pos <- sigPA_pea[sigPA_pea$log2FoldChange>0,]
#dim(sigPA_pea_pos)
# 0 6
#st_sigPA_pea_pos <- subset_taxa(prune_taxa(rownames(sigPA_pea_pos), physeq))
#st_sigPA_pea_pos
#seq_sigPA_pea_pos <- as.data.frame(tax_table(st_sigPA_pea_pos))
#Creating Reduced Tax Tables for this variable
#write.table(seq_sigPA_pea_pos, file="Reduced Tax Tables/PresAbs_PEA_POS.csv", quote=FALSE,sep = ",", col.names=FALSE)


#NEGATIVE - WILL BE THE SAME AS THE COMPLETE DATA SET 

sigPA_pea_neg <- sigPA_pea[sigPA_pea$log2FoldChange<0,]
dim(sigPA_pea_neg)
# 66 6

st_sigPA_pea_neg <- subset_taxa(prune_taxa(rownames(sigPA_pea_neg), physeq))
st_sigPA_pea_neg

seq_sigPA_pea_neg <- as.data.frame(tax_table(st_sigPA_pea_neg))

write.table(seq_sigPA_pea_neg, file="Reduced Tax Tables/PresAbs_PEA_NEG.csv", quote=FALSE,sep = ",", col.names=T)



###RFTM_pa1.peacrabs.f1####
resPA_rftmpea <- results(ddsPA_rftmpeaM, name="RFTM_pa1.peacrabs.f1")
sigPA_rftmpea <- resPA_rftmpea[which(resPA_rftmpea$padj < 0.05), ]
dim(sigPA_rftmpea)
#33 6

st_sigPA_rftmpea <- subset_taxa(prune_taxa(rownames(sigPA_rftmpea), physeq))
st_sigPA_rftmpea

seq_sigPA_rftmpea <- as.data.frame(tax_table(st_sigPA_rftmpea))

#Creating Reduced Tax Tables for this variable
write.table(seq_sigPA_rftmpea, file="Reduced Tax Tables/PresAbs_RFTMPEA_COMP.csv", quote=FALSE,sep = ",", col.names=T)


#Getting the +/- Log of this variable

#POSITIVE
sigPA_rftmpea_pos <- sigPA_rftmpea[sigPA_rftmpea$log2FoldChange>0,]
dim(sigPA_rftmpea_pos)
# 25 6

st_sigPA_rftmpea_pos <- subset_taxa(prune_taxa(rownames(sigPA_rftmpea_pos), physeq))
st_sigPA_rftmpea_pos

seq_sigPA_rftmpea_pos <- as.data.frame(tax_table(st_sigPA_rftmpea_pos))

#Creating Reduced Tax Tables for this variable
write.table(seq_sigPA_rftmpea_pos, file="Reduced Tax Tables/PresAbs_RFTMPEA_POS.csv", quote=FALSE,sep = ",", col.names=T)


#NEGATIVE
sigPA_rftmpea_neg <- sigPA_rftmpea[sigPA_rftmpea$log2FoldChange<0,]
dim(sigPA_rftmpea_neg)
# 8 6

st_sigPA_rftmpea_neg <- subset_taxa(prune_taxa(rownames(sigPA_rftmpea_neg), physeq))
st_sigPA_rftmpea_neg

seq_sigPA_rftmpea_neg <- as.data.frame(tax_table(st_sigPA_rftmpea_neg))

write.table(seq_sigPA_rftmpea_neg, file="Reduced Tax Tables/PresAbs_RFTMPEA_NEG.csv", quote=FALSE,sep = ",", col.names=T)


#***************************************************************************************************************************************

#Evaluating Site */+ Peacrabs and RFTM ####

ddsPA_rftmsiteM <- phyloseq_to_deseq2(physeq, ~ RFTM_pa * Site.x)
ddsPA_rftmsiteM <- DESeq(ddsPA_rftmsiteM, test="Wald", fitType="parametric")
resultsNames(ddsPA_rftmsiteM)
#RFTM_pa_1_vs_0
#Site.x_OY_vs_NW
#Site.x_SW_vs_NW
#RFTM_pa1.Site.xOY
#RFTM_pa1.Site.xSW

ddsPA_rftmsiteP <- phyloseq_to_deseq2(physeq, ~ RFTM_pa + Site.x)
ddsPA_rftmsiteP <- DESeq(ddsPA_rftmsiteP, test="Wald", fitType="parametric")
resultsNames(ddsPA_rftmsiteP)
#RFTM_pa_1_vs_0
#Site.x_OY_vs_NW
#Site.x_SW_vs_NW


ddsPA_PeasiteP <- phyloseq_to_deseq2(physeq, ~ peacrabs.f + Site.x)
ddsPA_PeasiteP <- DESeq(ddsPA_PeasiteP, test="Wald", fitType="parametric")
resultsNames(ddsPA_PeasiteP)
#peacrabs.f_1_vs_0
#Site.x_OY_vs_NW
#Site.x_SW_vs_NW

RFTM_pa1.peacrabs.f1


ddsPA_PeasiteM <- phyloseq_to_deseq2(physeq, ~ peacrabs.f * Site.x)
ddsPA_PeasiteM <- DESeq(ddsPA_PeasiteM, test="Wald", fitType="parametric")
resultsNames(ddsPA_PeasiteM)
#peacrabs.f_1_vs_0
#Site.x_OY_vs_NW
#Site.x_SW_vs_NW
#peacrabs.f1.Site.xOY
#peacrabs.f1.Site.xSW"



ddsPA_rftmsiteM
##Site.x_OY_vs_NW
##Site.x_SW_vs_NW
##RFTM_pa1.Site.xOY
##RFTM_pa1.Site.xSW
ddsPA_PeasiteM
##Site.x_OY_vs_NW -- 
##Site.x_SW_vs_NW --- 
##peacrabs.f1.Site.xOY
##peacrabs.f1.Site.xSW




###Site.x_OY_vs_NW####
resPA_siteOYNW <- results(ddsPA_rftmsiteM, name="Site.x_OY_vs_NW")
sigPA_siteOYNW <- resPA_siteOYNW[which(resPA_siteOYNW$padj < 0.05), ]
dim(sigPA_siteOYNW)
# 313 6
st_sigPA_siteOYNW <- subset_taxa(prune_taxa(rownames(sigPA_siteOYNW), physeq))
seq_sigPA_siteOYNW <- as.data.frame(tax_table(st_sigPA_siteOYNW))
write.table(seq_sigPA_siteOYNW, file="Reduced Tax Tables/PresAbs_SITEOYNW_COMP.csv", quote=FALSE,sep = ",", col.names=T)

#POSITIVE
sigPA_siteOYNW_pos <- sigPA_siteOYNW[sigPA_siteOYNW$log2FoldChange>0,]
dim(sigPA_siteOYNW_pos)
# 176 6
st_sigPA_siteOYNW_pos <- subset_taxa(prune_taxa(rownames(sigPA_siteOYNW_pos), physeq))
seq_sigPA_siteOYNW_pos <- as.data.frame(tax_table(st_sigPA_siteOYNW_pos))
write.table(seq_sigPA_siteOYNW_pos, file="Reduced Tax Tables/PresAbs_SITEOYNW_POS.csv", quote=FALSE,sep = ",", col.names=T)

#NEGATIVE
sigPA_siteOYNW_neg <- sigPA_siteOYNW[sigPA_siteOYNW$log2FoldChange<0,]
dim(sigPA_siteOYNW_neg)
# 137 6
st_sigPA_siteOYNW_neg <- subset_taxa(prune_taxa(rownames(sigPA_siteOYNW_neg), physeq))
seq_sigPA_siteOYNW_neg <- as.data.frame(tax_table(st_sigPA_siteOYNW_neg))
write.table(seq_sigPA_siteOYNW_neg, file="Reduced Tax Tables/PresAbs_SITEOYNW_NEG.csv", quote=FALSE,sep = ",", col.names=T)




###Site.x_SW_vs_NW####
resPA_siteSWNW <- results(ddsPA_rftmsiteM, name="Site.x_SW_vs_NW")
sigPA_siteSWNW <- resPA_siteSWNW[which(resPA_siteSWNW$padj < 0.05), ]
dim(sigPA_siteSWNW)
# 295 6
st_sigPA_siteSWNW <- subset_taxa(prune_taxa(rownames(sigPA_siteSWNW), physeq))
seq_sigPA_siteSWNW <- as.data.frame(tax_table(st_sigPA_siteSWNW))
write.table(seq_sigPA_siteSWNW, file="Reduced Tax Tables/PresAbs_SITESWNW_COMP.csv", quote=FALSE,sep = ",", col.names=T)

#POSITIVE
sigPA_siteSWNW_pos <- sigPA_siteSWNW[sigPA_siteSWNW$log2FoldChange>0,]
dim(sigPA_siteSWNW_pos)
# 185 6
st_sigPA_siteSWNW_pos <- subset_taxa(prune_taxa(rownames(sigPA_siteSWNW_pos), physeq))
seq_sigPA_siteSWNW_pos <- as.data.frame(tax_table(st_sigPA_siteSWNW_pos))
write.table(seq_sigPA_siteSWNW_pos, file="Reduced Tax Tables/PresAbs_SITESWNW_POS.csv", quote=FALSE,sep = ",", col.names=T)

#NEGATIVE
sigPA_siteSWNW_neg <- sigPA_siteSWNW[sigPA_siteSWNW$log2FoldChange<0,]
dim(sigPA_siteSWNW_neg)
# 110 6
st_sigPA_siteSWNW_neg <- subset_taxa(prune_taxa(rownames(sigPA_siteSWNW_neg), physeq))
seq_sigPA_siteSWNW_neg <- as.data.frame(tax_table(st_sigPA_siteSWNW_neg))
write.table(seq_sigPA_siteSWNW_neg, file="Reduced Tax Tables/PresAbs_SITESWNW_NEG.csv", quote=FALSE,sep = ",", col.names=T)



###RFTM_pa1.Site.xOY####
resPA_RFTMOY <- results(ddsPA_rftmsiteM, name="RFTM_pa1.Site.xOY")
sigPA_RFTMOY <- resPA_RFTMOY[which(resPA_RFTMOY$padj < 0.05), ]
dim(sigPA_RFTMOY)
# 184 6
st_sigPA_RFTMOY <- subset_taxa(prune_taxa(rownames(sigPA_RFTMOY), physeq))
seq_sigPA_RFTMOY <- as.data.frame(tax_table(st_sigPA_RFTMOY))
write.table(seq_sigPA_RFTMOY, file="Reduced Tax Tables/PresAbs_RFTMOY_COMP.csv", quote=FALSE,sep = ",", col.names=T)

#POSITIVE
sigPA_RFTMOY_pos <- sigPA_RFTMOY[sigPA_RFTMOY$log2FoldChange>0,]
dim(sigPA_RFTMOY_pos)
# 111 6
st_sigPA_RFTMOY_pos <- subset_taxa(prune_taxa(rownames(sigPA_RFTMOY_pos), physeq))
seq_sigPA_RFTMOY_pos <- as.data.frame(tax_table(st_sigPA_RFTMOY_pos))
write.table(seq_sigPA_RFTMOY_pos, file="Reduced Tax Tables/PresAbs_RFTMOY_POS.csv", quote=FALSE,sep = ",", col.names=T)

#NEGATIVE
sigPA_RFTMOY_neg <- sigPA_RFTMOY[sigPA_RFTMOY$log2FoldChange<0,]
dim(sigPA_RFTMOY_neg)
# 73 6
st_sigPA_RFTMOY_neg <- subset_taxa(prune_taxa(rownames(sigPA_RFTMOY_neg), physeq))
seq_sigPA_RFTMOY_neg <- as.data.frame(tax_table(st_sigPA_RFTMOY_neg))
write.table(seq_sigPA_RFTMOY_neg, file="Reduced Tax Tables/PresAbs_RFTMOY_NEG.csv", quote=FALSE,sep = ",", col.names=T)




###RFTM_pa1.Site.xSW####
resPA_RFTMSW <- results(ddsPA_rftmsiteM, name="RFTM_pa1.Site.xSW")
sigPA_RFTMSW <- resPA_RFTMSW[which(resPA_RFTMSW$padj < 0.05), ]
dim(sigPA_RFTMSW)
# 209 6
st_sigPA_RFTMSW <- subset_taxa(prune_taxa(rownames(sigPA_RFTMSW), physeq))
seq_sigPA_RFTMSW <- as.data.frame(tax_table(st_sigPA_RFTMSW))
write.table(seq_sigPA_RFTMSW, file="Reduced Tax Tables/PresAbs_RFTMSW_COMP.csv", quote=FALSE,sep = ",", col.names=T)

#POSITIVE
sigPA_RFTMSW_pos <- sigPA_RFTMSW[sigPA_RFTMSW$log2FoldChange>0,]
dim(sigPA_RFTMSW_pos)
# 124 6
st_sigPA_RFTMSW_pos <- subset_taxa(prune_taxa(rownames(sigPA_RFTMSW_pos), physeq))
seq_sigPA_RFTMSW_pos <- as.data.frame(tax_table(st_sigPA_RFTMSW_pos))
write.table(seq_sigPA_RFTMSW_pos, file="Reduced Tax Tables/PresAbs_RFTMSW_POS.csv", quote=FALSE,sep = ",", col.names=T)

#NEGATIVE
sigPA_RFTMSW_neg <- sigPA_RFTMSW[sigPA_RFTMSW$log2FoldChange<0,]
dim(sigPA_RFTMSW_neg)
# 85 6
st_sigPA_RFTMSW_neg <- subset_taxa(prune_taxa(rownames(sigPA_RFTMOY_neg), physeq))
seq_sigPA_RFTMSW_neg <- as.data.frame(tax_table(st_sigPA_RFTMSW_neg))
write.table(seq_sigPA_RFTMSW_neg, file="Reduced Tax Tables/PresAbs_RFTMSW_NEG.csv", quote=FALSE,sep = ",", col.names=T)


###peacrabs.f1.Site.xOY####
resPA_PeaOY <- results(ddsPA_PeasiteM, name="peacrabs.f1.Site.xOY")
sigPA_PeaOY <- resPA_PeaOY[which(resPA_PeaOY$padj < 0.05), ]
dim(sigPA_PeaOY)
# 52 6 
st_sigPA_PeaOY <- subset_taxa(prune_taxa(rownames(sigPA_PeaOY), physeq))
seq_sigPA_PeaOY <- as.data.frame(tax_table(st_sigPA_PeaOY))
write.table(seq_sigPA_PeaOY, file="Reduced Tax Tables/PresAbs_PEAOY_COMP.csv", quote=FALSE,sep = ",", col.names=T)

#POSITIVE
sigPA_PeaOY_pos <- sigPA_PeaOY[sigPA_PeaOY$log2FoldChange>0,]
dim(sigPA_PeaOY_pos)
# 29 6
st_sigPA_PeaOY_pos <- subset_taxa(prune_taxa(rownames(sigPA_PeaOY_pos), physeq))
seq_sigPA_PeaOY_pos <- as.data.frame(tax_table(st_sigPA_PeaOY_pos))
write.table(seq_sigPA_PeaOY_pos, file="Reduced Tax Tables/PresAbs_PEAOY_POS.csv", quote=FALSE,sep = ",", col.names=T)

#NEGATIVE
sigPA_PeaOY_neg <- sigPA_PeaOY[sigPA_PeaOY$log2FoldChange<0,]
dim(sigPA_PeaOY_neg)
# 23 6
st_sigPA_PeaOY_neg <- subset_taxa(prune_taxa(rownames(sigPA_PeaOY_neg), physeq))
seq_sigPA_PeaOY_neg <- as.data.frame(tax_table(st_sigPA_PeaOY_neg))
write.table(seq_sigPA_PeaOY_neg, file="Reduced Tax Tables/PresAbs_PEAOY_NEG.csv", quote=FALSE,sep = ",", col.names=T)





###peacrabs.f1.Site.xSW####
resPA_PeaSW <- results(ddsPA_PeasiteM, name="peacrabs.f1.Site.xSW")
sigPA_PeaSW <- resPA_PeaSW[which(resPA_PeaSW$padj < 0.05), ]
dim(sigPA_PeaSW)
# 37 6
st_sigPA_PeaSW <- subset_taxa(prune_taxa(rownames(sigPA_PeaSW), physeq))
seq_sigPA_PeaSW <- as.data.frame(tax_table(st_sigPA_PeaSW))
write.table(seq_sigPA_PeaSW, file="Reduced Tax Tables/PresAbs_PEASW_COMP.csv", quote=FALSE,sep = ",", col.names=T)

#POSITIVE
sigPA_PeaSW_pos <- sigPA_PeaSW[sigPA_PeaSW$log2FoldChange>0,]
dim(sigPA_PeaSW_pos)
# 8 6
st_sigPA_PeaSW_pos <- subset_taxa(prune_taxa(rownames(sigPA_PeaSW_pos), physeq))
seq_sigPA_PeaSW_pos <- as.data.frame(tax_table(st_sigPA_PeaSW_pos))
write.table(seq_sigPA_PeaSW_pos, file="Reduced Tax Tables/PresAbs_PEASW_POS.csv", quote=FALSE,sep = ",", col.names=T)

#NEGATIVE
sigPA_PeaSW_neg <- sigPA_PeaSW[sigPA_PeaSW$log2FoldChange<0,]
dim(sigPA_PeaSW_neg)
# 29 6
st_sigPA_PeaSW_neg <- subset_taxa(prune_taxa(rownames(sigPA_PeaSW_neg), physeq))
seq_sigPA_PeaSW_neg <- as.data.frame(tax_table(st_sigPA_PeaSW_neg))
write.table(seq_sigPA_PeaSW_neg, file="Reduced Tax Tables/PresAbs_PEASW_NEG.csv", quote=FALSE,sep = ",", col.names=T)


