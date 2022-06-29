#



#Import reduced tax table
RFTM_comp <- read.csv("Reduced Tax Tables/PresAbs_RFTM_COMP.csv") 
#rownames(RFTM_comp) <- RFTM_comp[,1]
#RFTM_comp <- RFTM_comp[ -c(1) ]
view(RFTM_comp)

#Create and import S4 file for Sig of the resudec tax table 
write.table(sigPA_rftm, file="PEA_RFTM/S4_SigPA_RFTM_COMP.csv", quote=FALSE,sep = ",", col.names=TRUE)
sigPA_RFTM <- read.csv("PEA_RFTM/S4_SigPA_RFTM_COMP.csv") 
view(sigPA_RFTM)

#Removing all other columns besides Phylogeny and log2fold change value
RFTM_pa <- merge(RFTM_comp, sigPA_RFTM,  by = 'row.names', all = TRUE)
RFTM_pa <- RFTM_pa[ -c(2, 9, 11, 13:16) ]
view(RFTM_pa)






