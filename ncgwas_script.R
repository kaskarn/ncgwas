#!/usr/bin/env Rscript 


##################################################################################### 
              # Report bugs and issues to baldassa@email.unc.edu # 
##################################################################################### 

##################################################################################### 
#############################    DESCRIPTION     ####################################
#####################################################################################
#
# Runs GWAS on NCDF4-formatted genotype dosage files. Written with WHI analyses in mind, 
# although this should work with any data so long as it is in NCDF4 format.
# outputs have the following columns for each SNP:
#   index -- order of appearance of SNPs in file
#   snp -- snp name
#   coded -- coded allele
#   other -- referent allele
#   caf -- allele frequency
#   v -- dosage variance
#   n -- nonmissing observations for SNP
#   b -- SNP beta
#   se -- standard error on b
#   p -- p-value
#   conv (GLM only) -- 1/0 indicator of model convergence

##################################################################################### 
#############################       USAGE        ####################################
#####################################################################################
#
#   --source=SOURCE FILE
#   R file to source before options are parsed. Allows specifying inputs 
#   in an R scripts rather than, or in addition to the command line. Optional.
#   
#   -p PHENOTYPE FILE, --pheno=PHENOTYPE FILE
#   Path to phenotype file in csv, tab-delimited or space-delimited format. Required.
#   
#   -r RESULTS, --resdir=RESULTS
#   Results directory path. Default is folder "ncgwas_results" in working directory
#   
#   -g GENE DATA, --gpath=GENE DATA
#   Path to NCDF genetic data. Default is /nas02/depts/epi/Genetic_Data_Center/whi_share/whi_1000g_fh_imp/ncdf-data/
#     
#     -s STUDY, --study=STUDY
#   WHI study. Required.
#   
#   -o OUTCOME, --outcome=OUTCOME
#   Outcome variable to be used in models. Required.
#   
#   -f FORMULA, --form=FORMULA
#   Right-hand side of the model, starting with ~ in R-style. Required.
#   
#   -m MODEL, --model=MODEL
#   Model type: linear or GLM. Default is linear.
#   
#   --glmfam=GLM FAMILY
#   For GLMs: Distribution family of outcome (binomial, Gamma, gaussian,...) see ?family. 
#   Default is binomial if model is set to "glm"
#   
#   --link=GLM LINK
#   For GLMs: Link function. Default is NULL if model is set to "glm"
#   
#   -i ID VARIABLE, --idvar=ID VARIABLE
#   Name of ID variable in phenotype file. Default is leftmost variable in phenotype file
#   
#   -x MIN CAF, --mincaf=MIN CAF
#   Minimum allele frequency required for inclusion. Default is 0.01
#   
#   -c CHROMOSOME(S), --chr=CHROMOSOME(S)
#   Chromosomes to run the GWAS on, specified as an expression for an R vector, e.g. 1:22, c(1,3,22)
#   
#   --norun=NORUN
#   Stop before running the GWAS, to inspect log for debugging purposes
#   
#   -h, --help
#   Show this help message and exit
#   
#   
#
##################################################################################### 
#############################    INSTRUCTIONS    ####################################
##################################################################################### 
#
# [[ 0 ]] Install the development version of data.table from github, using
# devtool::install_github("Rdatatable/data.table", build_vignettes=FALSE)
#
#
#
# [[ 1 ]] Copy the file /nas02/apps/r-3.3.1/lib64/R/OFED-1.5.2/library/Rmpi/Rprofile into your
# working directory; rename it to ".Rprofile". This will ensure that R is launched with
# the proper MPI-enabled profile from the script's directory. If the file is somehow lost,
# it may also be downloaded from https://web.stanford.edu/group/farmshare/cgi-bin/wiki/index.php/Rmpi
#
#
#
# [[ 2 ]] Place this script file in your working directory. Make sure to have the latest version -- updates
# are pushed to https://github.com/kaskarn/ncgwas. 
#
#
#
# [[ 3 ]] Set input arguments (quotes are not needed, e.g. -o jt and -o "jt" will have the same effect) 
# 
## Inputs can be set by sourcing an R file with the --source options, by listing them individually, or 
# by a combination of these two ways. When an argument is provided both in the command line and in a 
# sourced file, the command line value prevails. This allows to quickly run models changing only a few 
# arguments at a time. An example input source file follows: 
#
# example_inputs.R: (Quotes are now necessary to designate character variables since this is sourced by R)
#
#   phenodir  <- "/nas02/depts/epi/CVDGeneNas/antoine/ECG_GWAS/WHI/phenotypes/"
#   pheno     <- paste0(phenodir,"ecg_whi_whites_fi.txt")
#   resdir    <- "/proj/epi/CVDGeneNas/antoine/dev_garnetmpi/test_results/"
#   gpath     <- "/nas02/depts/epi/Genetic_Data_Center/whi_share/whi_1000g_fh_imp/ncdf-data/"
#   study     <- "GARNET"
#   outcome   <- "qt" 
#   form      <- "~g+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+region+rr_d+age"
#   chr       <-  22
#   idvar     <- "id"
#
# end example_inputs.R
#
# Since the unspecified parameters have default values (see usage), this will run fine by calling
# ./ncgwas_script.R --source inputs.R.  This is useful to run different models of a base sourced file,
# only overriding a few arguments at the time in the command line, e.g. changing the outcome with -o new_outcome, 
# or the model with -m glm without needing to write a new inputs.R source file
#
#
#
# [[ 4 ]] Check everything looks right by running (in or out of bsub) ./ncgwas_script.R with all your intended arguments
# followed by --norun. This will provide information on the run without going through the actual computations.
#
#
#
# [[ 5 ]] run on the LSF platform, with mpirun: for example:
# 
# ex: bsub -n 30 -M 6 -o "test.txt" mpirun ./ncgwas_script.R --source inputs.R -m glm --glmfam Gamma --link log

###################################################################################### 
############################     START OF SCRIPT      ################################
###################################################################################### 

#Load libraries
library(Rcpp)
library(MASS)
library(Matrix)
library(Rmpi)
library(optparse)
library(data.table)
library(ncdf4)
library(RcppEigen)
library(speedglm)

#####################################################################################
################################## Parsing step  ####################################
#####################################################################################

#Default values not specified here to avoid overriding --source file with default values
option_list = list(
  make_option("--source", type = "character", 
              help = ".R file to source before options are parsed. Allows specifying inputs 
              in an R scripts rather than, or in addition to the command line. Optional.",
              metavar = "source file"),
  make_option(c("-p", "--pheno"), type = "character", 
              help = "Path to phenotype file in csv, tab-delimited or space-delimited format. Required.", metavar = "phenotype file"),
  make_option(c("-r", "--resdir"), type = "character", 
              help = "Results directory path. Default is folder \"ncgwas_results\" in working directory", metavar = "results"),
  make_option(c("-g", "--gpath"), type = "character",
              help = "Path to NCDF genetic data. Default is /nas02/depts/epi/Genetic_Data_Center/whi_share/whi_1000g_fh_imp/ncdf-data/", metavar = "Gene data"),
  make_option(c("-s", "--study"), type = "character", 
              help = "WHI study. Required.", metavar = "study"),	
  make_option(c("-o", "--outcome"), type = "character", 
              help = "Outcome variable to be used in models. Required.", metavar = "Outcome"),	
  make_option(c("-f", "--form"), type = "character", 
              help = "Right-hand side of the model, starting with ~ in R-style. Required.", metavar = "Formula"),	
  make_option(c("-m", "--model"), type = "character",
              help = "Model type: linear or GLM. Default is linear.", metavar = "Model"),
  make_option("--glmfam", type = "character",
              help = "For GLMs: Distribution family of outcome (binomial, Gamma, gaussian,...) see ?family. 
              Default is binomial if model is set to \"glm\"", metavar = "GLM Family"),
  make_option("--link", type = "character",
              help = "For GLMs: Link function. Default is NULL if model is set to \"glm\"", metavar = "GLM Link"),
  make_option(c("-i", "--idvar"), type = "character",
              help = "Name of ID variable in phenotype file. Default is leftmost variable in phenotype file", metavar = "ID VARIABLE"),
  make_option(c("-x", "--mincaf"), type = "double",
              help = "Minimum allele frequency required for inclusion. Default is 0 (keep all)", metavar = c("MIN CAF")),
  make_option(c("-c", "--chr"), type = "character",
              help = "Chromosomes to run the GWAS on, specified as an expression for an R vector, e.g. 1:22, c(1,3,22)",
              metavar = "chromosome(s)"),
  make_option(c("--norun"), type = "logical", action = "store_true", default = "FALSE",
              help = "Stop before running the GWAS, to inspect log for debugging purposes. norun must be specified from the command line",
              metavar = "norun")
)
cat("\n_____________________________________________________________________\n\n")
cat("_ __   ___ __ ___      ____ _ ___     ___  ___ _ __(_)_ __ | |_ \n")
cat("| '_ \\ / __/ _` \\ \\ /\\ / / _` / __|   / __|/ __| '__| | '_ \\| __|\n")
cat("| | | | (_| (_| |\\ V  V / (_| \\__ \\   \\__ \\ (__| |  | | |_) | |_ \n")
cat("|_| |_|\\___\\__, | \\_/\\_/\\__,_|___/___|___/\\___|_|  |_| .__/ \\__|\n")
cat("           |___/                 |_____|              |_|        \n")

cat("ncgwas_script.R\ncontact baldassa@email.unc.edu to complain about bugs or request feature\n\n\n")
opt_parser = OptionParser(usage = "%prog [options] file", option_list=option_list)
opt = parse_args(opt_parser)

if(!exists("opt")){
  cat("\nERROR:\tparse_args() returned an error\n\tdue to incorect arguments, see above error message. Check the program usage (--help)\n\n")
  if(mpi.comm.size() > 1) mpi.close.Rslaves()
  mpi.quit()
}

#Reads command-line options into like-named variables. Overrides --source for convenience
if(!is.null(opt$source)) source(opt$source)

#Print out the variables defined in the source script
cat("\n\nSource file arguments:\n")
for(i in option_list) {
  if(exists(i@dest) && !is.function(get(i@dest))) cat("\t", i@dest,":",get(i@dest),"\n")
}

cat("\nCommand line arguments:\n")
for(i in names(opt)) {
  cat("\t", i,":",opt[[i]])
  if(exists(i) && !is.function(get(i))) cat(" (OVERRIDEN)")
  cat("\n")
  assign(i, opt[[i]])
}
#Sets default if not specified in command line or --source file
if(!exists("gpath")) gpath <- "/nas02/depts/epi/Genetic_Data_Center/whi_share/whi_1000g_fh_imp/ncdf-data/"
if(!exists("resdir")) resdir <- "ncgwas_results"
if(!exists("model")) model <- "linear"
if(!exists("mincaf")) mincaf <- 1E-2
if(model == "glm") {if(!exists("glmfam")) glmfam <- "binomial"; if(!exists("link")) link <- NULL}


#Check nothing something important is missing
im <- c("outcome", "pheno")
if(sum(is.na(im_mis <- match(im, ls()))) > 0){
  cat("\nError: parameter(s): ", paste(im[is.na(im_mis)], collapse = ", "), "  must be provided\n")
  if(mpi.comm.size() > 1) mpi.close.Rslaves()
  mpi.quit()
}

#####################################################################################
#############################    End of parsing step  ###############################
#####################################################################################

#####################################################################################
#################################   Setup step   ####################################
#####################################################################################

#Turn formula input into formula object
form <- as.formula(form)

#Create directory/ies for results
dir.create(resdir,showWarnings = FALSE, recursive = TRUE)

#Load and pare down data
Epidata <- fread(pheno)
if(!exists("idvar")) idvar <- names(Epidata)[1]
invisible(Epidata[,g:=rnorm(nrow(Epidata))])
Epidata <- na.omit(Epidata[,c(idvar, outcome, all.vars(form)), with = F])
if(idvar != "Common_ID") setnames(Epidata, idvar, "Common_ID")
setkey(Epidata, Common_ID)

#order ids as in NC file
nc22 <- nc_open(paste0(gpath,study,'-chr',22,'-c.nc'))
ncids <- data.table(Common_ID = ncvar_get(nc22, "Common_ID"), 
                    ncid = seq_len(nc22$dim$Samples$len))
setkey(ncids,Common_ID)
dt_ana <- ncids[Epidata, nomatch = 0]
setkey(dt_ana, ncid)

#Record NCDF indices of participants with phenotype data
nckeep <- dt_ana$ncid

#Create X and y model matrices, for LMs
X <- as.matrix(dt_ana[,all.vars(form),with =FALSE][,int:=1])
y <- as.numeric(dt_ana[,get(outcome)])
nvar <- ncol(X)
gpos <- match("g", colnames(X))

#Finalize analytical datasets, for GLMs
dt_ana <- dt_ana[,c(outcome, all.vars(form)),with=FALSE]

#Number of workers, without the master thread
nworkers <- mpi.comm.size() - 1

#Fit functions: workhorse functions are #RcppEigen::fastLmPure for linear models 
#and speedglm::speedglm for generalized linear models
qfit_lm <- function(gnow){
  ind <- which(!is.na(gnow))
  X[,gpos] <- gnow
  tm <- fastLmPure(X[ind,],y[ind])
  list(tm$coefficients[gpos], tm$se[gpos])
}
qfit_glm <- function(gnow){
  ind <- which(!is.na(gnow))
  X[,gpos] <- gnow
  tm <- try(speedglm.wfit(y[ind],X[ind,],FALSE, 
                          family=do.call(glmfam,as.list(link)),
                          set.default=list(row.chunk=2000)), TRUE)
  if(class(tm) == "try-error") return(as.list(as.numeric(c(NA,NA,NA,0))))
  c(as.list(as.numeric(as.matrix(summary(tm)$coefficients[gpos,-3]))),1)
}

##### Print details of analyses to logfile #####
if(model == "linear"){
  cat("\nWill fit Linear models:\n")
}else if(model == "glm") {
  cat("\nWill fit Generalized Linear Models:\n")
  cat("\t Family:",glmfam,"\n")
  cat("\t Link:",link,"\n")
}
cat("\t Formula:", paste(form), "\n")
cat("\nDetails:\n")
cat("\t Chromosomes:", chr,"\n")
cat("\t Phenotype file N in NC files:", nrow(dt_ana),"\n")
cat("\t # of workers:", nworkers,"\n")
#fmem <- lapply(strsplit(system("free -m", intern = T), " "), function(i) i[i != ""])[[3]][4]

##### Error catching ####
if(nrow(dt_ana) < 10){
  cat("\nERROR: Only",nrow(dt_ana),"IDs found. Wrong ID variable specified, or wrong study IDs used.\n")
  if(mpi.comm.size() > 1) mpi.close.Rslaves()
  mpi.quit()
}
#Check model is linear or glm
if(!(model %in% c("linear", "glm"))){
  cat("\nERROR: model must be either linear, or glm.\n")
  if(mpi.comm.size() > 1) mpi.close.Rslaves()
  mpi.quit()
}
#Check MPI correctly setup
if(mpi.comm.size() <= 1){
  cat("\nERROR: too few MPI workers; did you call the script with mpirun, and the -n bsub option?\n\n")
  mpi.quit()
}

cat("\n\nALL LOOK OK\n\n")

###Stop here if debugging only ###
if(norun == TRUE){
  if(mpi.comm.size() > 1) mpi.close.Rslaves()
  mpi.quit()
}

#Split convenience function
splitup <- function(a, n) lapply(split(a[1]:a[2], cut(a[1]:a[2], n)), range)

#### MPI Setup ####

#Send ALL objects (because why not), and the needed libraries to worker threads
mpi.bcast.Robj2slave(all = TRUE) 
mpi.bcast.cmd({
  library(Rcpp); library(data.table)
  library(Matrix); library(MASS)
  library(RcppEigen); library(parallel)
  library(speedglm); library(ncdf4)
})



#Data for debugging

#####################################################################################
############################    End of setup step    ################################
#####################################################################################

#####################################################################################
#################################   LET'S ROLL   ####################################
#####################################################################################

#Start loop over chromosomes
for(i in chr){
  #misc: get #snps, make output file name, send chromosome # to workers...
  print(paste("Starting on chromosome",i,"at:",Sys.time()))
  mpi.bcast.Robj2slave(i)
  rname <- paste0(resdir,"Chr",i,"_",outcome,"_",study,"_results.csv")
  nsnp <- nc_open(paste0(gpath,study,'-chr',i,'-c.nc'))$dim$SNPs$len
  
  #splitup task into optimized # of chunks with even memory burden
  bits <- splitup(c(1,nsnp), ceiling(100/nworkers)*100)
  res <- mpi.parLapply(bits, function(k) {
    start <- k[1]; end <- k[2]; span <- k[2]-k[1]+1
    #Open nc file, get SNP names and create results dataset
    nc <- nc_open(paste0(gpath,study,'-chr',i,'-c.nc'))
    snp_names <- ncvar_get(nc, "SNP_Name",c(1,start), c(-1,span))
    res_part <- data.table(
      index = as.integer(start:end), 
      snp = snp_names, 
      coded = ncvar_get(nc,"Allele1_Reference", start, span), 
      other = ncvar_get(nc,"Allele2_Reference", start, span),
      caf = as.numeric(NA), b = as.numeric(NA), se = as.numeric(NA), 
      p = as.numeric(NA), j = seq_along(start:end))
    if(model == "glm") res_part[,conv := as.numeric(NA)]
    
    #Read dosages at relevant indices, and restrict to participants also in phenotype file
    p_aa <- ncvar_get(nc,"Prob_AA", start=c(start,1), count=c(span, -1))[,nckeep]
    p_ab <- ncvar_get(nc,"Prob_AB", start=c(start,1), count=c(span, -1))[,nckeep]
    dos <- p_aa*2 + p_ab
    
    #Add allele frequency, variance and nonomissing N
    res_part[,c("caf", "v", "n") := list(mean(dos[j,]/2, na.rm = TRUE), 
                                         var(dos[j,], na.rm = TRUE),
                                         sum(!is.nan(dos[j,]), na.rm = TRUE)), j]
    
    #Add regression results using qfit functions applied to every column of the dosage matrix,
    #wrapped with the data.table by= operator for speed. qfit_lm and qfit_glm are defined above.
    if(model == "linear"){ res_part[n > 0 & abs(1-caf) > mincaf & v > 0, c("b", "se") := qfit_lm(dos[j,]), j]
    }else res_part[n > 0 & abs(1-caf) > mincaf & v > 0, c("b", "se","p","conv") := qfit_glm(dos[j,]), j]
    
    #Return data.table copy to avoid memory leaks
    copy(res_part)
  })
  #For debugging
  if(!is.data.table(res[[1]])) print(res)
  warnings()
  
  #Combine worker outputs into single DT
  res <- rbindlist(res)
  
  #Get p-values for linear model (which don't already compute it)
  if(model == "linear") res[,p := 2*(1-pt(abs(b/se),n-nvar)) ]
  
  #Add chromosome column, and remove j-index column
  res[,chr := i]
  res[,j := NULL]
  
  #Write to file
  fwrite(res, rname, sep = ",")
  print(paste("Done with chromosome:",i,"at:",Sys.time()))
  print(paste("Size of chunked outputs: ", format(object.size(res), units = "MB")))
}

mpi.close.Rslaves()
mpi.quit()