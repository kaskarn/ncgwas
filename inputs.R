phenodir  <- "/nas02/depts/epi/CVDGeneNas/antoine/ECG_GWAS/WHI/phenotypes/" #phenotype file directory
pheno     <- paste0(phenodir,"ecg_whi_whites_fi.txt") #Full path of the .txt phenotype file
resdir    <- "/proj/epi/CVDGeneNas/antoine/dev_garnetmpi/test_results/" #where to store results
gpath     <- "/nas02/depts/epi/Genetic_Data_Center/whi_share/whi_1000g_fh_imp/ncdf-data/" #where the 1KG .nc files live
study     <- "GARNET" #the WHI study
outcome   <- "jt" #The outcome of interest
form      <- "~g+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+region+rr_d+age"

model     <- "linear" #"linear" or "glm"
#glmfam    <-  #Optional: family for GLM -- default is binomial.
#link      <-  #Optional: GLM link function -- default is family default (e.g. logit for binomial)
chr       <- 22 #Chromosomes to restrict to
idvar     <- "id"