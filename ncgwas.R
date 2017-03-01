#!/usr/bin/env Rscript 

# Please do not change this script without notifying Antoine Baldassari
# baldassa@email.unc.edu

###################################################################################### 
#######################        GET .Rprofile if needed         #######################
###################################################################################### 

# Throw error if Rprofile is missing, offer to install it if needed
if(!("package:Rmpi" %in% search()) && !file.exists(".Rprofile")){
  wd <- getwd()
  err <- 0
  rprofpath <- "/nas02/apps/r-3.3.1/lib64/R/OFED-1.5.2/library/Rmpi/Rprofile"
  cat("\nIncorrect R profile. Make sure the right .Rprofile is in your working directory.")
  cat("\nThe correct file is /nas02/apps/r-3.3.1/lib64/R/OFED-1.5.2/library/Rmpi/Rprofile")
  cat("\nand must be renamed to .Rprofile once copied into your working directory.\n\n")
  cat("Install Rprofile in working directory now? [Y/N]\n")
  doit <- readLines(con = "stdin", 1)
  if(doit == "Y"){
    if(!file.exists(rprofpath)){
      cat("\nFile was moved from",rprofpath,"\n")
      err <- 1
    }else if(file.exists(".Rprofile")){
      cat("\nFile named \".Rprofile\" already exists in your working directory")
      cat("\nAre you sure you want to overwrite it? [Y/N]\n")
      sure <- readLines(con = "stdin", 1)
      if(sure == "N") err <- 1
    }
    if(err == 0){
      cat("Copying Rprofile to working directory...\n")
      file.copy(rprofpath, tempdir())
      setwd(tempdir())
      file.rename("Rprofile", ".Rprofile")
      file.copy(".Rprofile", wd)
      file.remove(".Rprofile")
      setwd(wd)
      cat("Done.\n")
    }
  }else err <- 1
  
  if(err == 1) cat("\nContact antoine baldassa@email.unc.edu if you need assistance\n\n")
  if(err == 0) cat("\nncgwas_script should now work. Try running the script again.\n")
  q()
}


###################################################################################### 
############################        LIBRARIES         ################################
###################################################################################### 

#Load libraries
library(Rcpp)
library(MASS)
library(Matrix)
library(optparse)
library(data.table)
library(ncdf4)
library(RcppEigen)
library(speedglm)

#####################################################################################
##############################     Parsing step      ################################
#####################################################################################

#Default values not specified here to avoid overriding --source file with default values
option_list = list(
  make_option("--source", type = "character", help = "
		R file to source before options are parsed. Allows specifying inputs
		in an R scripts rather than, or in addition to the command line.",
              metavar = "file"),
  make_option(c("--pheno"), type = "character", help = "[Required]
		Path to phenotype file. Allowed types:
		  a. tab-, comma-, space- delimited files, with missing data marked as NA
		  b. SAS .sasb7dat files,
		  c. Stata .dat files", metavar = "file"),
  make_option(c("-s", "--study"), type = "character", help = "[Required]	
		WHI study; case-sensitive.", metavar = "study"),	
  make_option(c("-o", "--outcome"), type = "character",help = "[Required]	
		Outcome variable.", metavar = "varname"),			  
  make_option(c("-f", "--form"), type = "character", help = "[Required]	
		Right-hand side of the model written in R-formula style, starting as: ~ g + ... ", metavar = "Formula"),
  make_option(c("-m", "--model"), type = "character", help = "[Default: linear]	
		Model type: either linear or glm.", metavar = "model"),		
  make_option("--glmfam", type = "character",help = "[GLM option. Default: binomial]		
		For GLMs: Distribution family of outcome
		For a list of possible families and their default link function, see ?family (in R). ", metavar = "GLM family"),
  make_option("--glmlink", type = "character",help = "[GLM option. Default: R default for family]		
		For GLMs: Link function.", metavar = "GLM Link"),
  make_option(c("-r", "--resdir"), type = "character", help = "[Default: working directory]
		Path to results directory. Will create new folders along the path if needed.", metavar = "PATH"),
  make_option(c("-g", "--gpath"), type = "character", 
    help = "[Default: /nas02/depts/epi/Genetic_Data_Center/whi_share/whi_1000g_fh_imp/ncdf-data/]	
		Gene data directory", metavar = "PATH"),
  make_option(c("-i", "--idvar"), type = "character",help = "[Default: name of first column in phenotype file]		
		Name of ID variable in phenotype file.", metavar = "varname"),
  make_option(c("-x", "--mincaf"), type = "double",help = "[Default: 0]		
		Minimum minor allele frequency cutoff for inclusion. Should be between 0 and 1", metavar = "NUM"),
  make_option(c("-c", "--chr"), type = "character",help = "[Default: 1:22]		
		Chromosomes to run the GWAS on, specified as an expression for an R vector
		e.g. 1:22, c(1,3,22)", metavar = "chromosome(s)"),
  make_option(c("-a", "--annodir"), type = "character", 
    help = "[Default: /nas02/depts/epi/Genetic_Data_Center/whi_share/whi_1000g_fh_imp/snp-summary-files/]	
		Annotation data directory", metavar = "PATH"),
  make_option("--nomerge", type = "logical", action = "store_true", 
    help = "Separate output by chromosome"),
  make_option("--noanno", type = "logical", action = "store_true",
    help = "Don't merge with annotation file (has no effect for now)"),
  make_option(c("--norun"), type = "logical", action = "store_true", default = FALSE,
    help = "Stop after printing out summary of model. Useful for checking all arguments are as 
		intended before running the GWAS. 

		Important: --norun must be specified from the command line, and not from a sourced .R file")
)

opt_parser = OptionParser(usage = "%prog  [--pheno file] [--study name] [--outcome name] [--form formula]
                          [--source file] [--resdir path] [--gpath path] [--model linear | glm] 
                          [--idvar varname] [--mincaf number] [--chr R-expression] [--norun] [--help]",
                          option_list=option_list)
opt = parse_args(opt_parser)

cat("\n\n_ __   ___ __ ___      ____ _ ___     ___  ___ _ __(_)_ __ | |_ \n")
cat("| '_ \\ / __/ _` \\ \\ /\\ / / _` / __|   / __|/ __| '__| | '_ \\| __|\n")
cat("| | | | (_| (_| |\\ V  V / (_| \\__ \\   \\__ \\ (__| |  | | |_) | |_ \n")
cat("|_| |_|\\___\\__, | \\_/\\_/\\__,_|___/___|___/\\___|_|  |_| .__/ \\__|\n")
cat("           |___/                 |_____|              |_|        \n")
cat("\n_____________________________________________________________________\n")
cat("\nncgwas.R\ncontact baldassa@email.unc.edu to complain about bugs or request features\n")

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

#Then those defined in the command line
cat("\nCommand line arguments:\n")
for(i in names(opt)) {
  cat("\t", i,":",opt[[i]])
  if(exists(i) && !is.function(get(i))) cat("  (OVERRIDES variable in ",basename(opt$source),")",sep = "")
  cat("\n")
  assign(i, opt[[i]])
}

#Check nothing important is missing
im <- c("outcome", "pheno", "form", "study")
if(sum(im_mis <- !sapply(im, exists)) != 0){
  cat("\nERROR: missing arguments: ", paste(im[im_mis], collapse = ", "), "\n")
  if(mpi.comm.size() > 1) mpi.close.Rslaves()
  mpi.quit()
}

#Sets defaults if not specified in command line or --source file, and log them
cat("\nDefault arguments used:\n")
def <- list(
  resdir = "ncgwas_results", model = "linear",
  gpath = "/nas02/depts/epi/Genetic_Data_Center/whi_share/whi_1000g_fh_imp/ncdf-data/",
  nomerge = FALSE, noanno = FALSE, 
  annodir = "/nas02/depts/epi/Genetic_Data_Center/whi_share/whi_1000g_fh_imp/snp-summary-files/",
  glmfam = "binomial", glmlink = substitute(get(glmfam)()$link), mincaf = 0
)
if(exists("model") & model != "glm") def$glmfam <- def$glmlink <- NULL
for(i in sort(names(def))) {
  if(!exists(i)) {assign(i, eval(def[[i]])) ; cat("\t ",i,": ",get(i),"\n",sep="")}
}

err <- 0
#Check phenotype file exists
if(!file.exists(pheno)){
  cat("\nERROR: Phenotype file does not exist")
  err <- 1;
}

#Check gene folder exists and contains study
if(dir.exists(gpath)){
  if(length(grep(study,list.files(gpath))) == 0){
    cat("\nERROR: study not found in gene data.")
    err <- 1
  }
}else{
  cat("\nERROR: Gene directory does not exist")
  err <- 1;
}

#Check annotation directory exists
if(!dir.exists(annodir)){
  cat("\nERROR: Annotation files directory does not exist")
  err <- 1;
}

#Check model is linear or glm
if(!(model %in% c("linear", "glm"))){
  cat("\nERROR: model must be either linear, or glm.");
  err <- 1
}
#Quit if any error so far
if(err > 0){
  cat("\n\nExiting due to errors\n")
  if(mpi.comm.size() > 1) mpi.close.Rslaves()
  mpi.quit()  
}
#####################################################################################
#################################   Setup step   ####################################
#####################################################################################

#Create directory/ies for results
dir.create(resdir,showWarnings = FALSE, recursive = TRUE)

#Turn formula and chromosome inputs into right formats
if(is.character(form)) form <- as.formula(form)
if(is.character(chr)) chr <- eval(parse(text = chr))

#Load and pare down data
pext <- strsplit(pheno, "\\.")[[1]][2]
if(pext == "sas7bdat"){ 
  require(sas7bdat)
  Epidata <- as.data.table(read.sas7bdat(phen))
}else if(pext == "dat") {
  require(readstata13)
  Epidata <- as.data.table(read.dta13(phen))
}else Epidata <- fread(pheno)

#Set idvar to first column by default
if(!exists("idvar")) idvar <- names(Epidata)[1]
#Create dummy g variable
invisible(Epidata[,g:=rnorm(nrow(Epidata))])
#Only keep variables in formula, with full data
Epidata <- na.omit(Epidata[,c(idvar, outcome, all.vars(form)), with = F])
#Rename ID to "Common_ID" if not already the case, and sort by it.
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

#Define fit functions: workhorse functions are #RcppEigen::fastLmPure for linear models 
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
                          family=do.call(glmfam,as.list(glmlink)), #get(glmfam)(glmlink) maybe clearer
                          set.default=list(row.chunk=2000)), TRUE) 
  if(class(tm) == "try-error") return(as.list(as.numeric(c(NA,NA,NA,0))))
  c(as.list(as.numeric(as.matrix(summary(tm)$coefficients[gpos,-3]))),1)
}

#####################################################################################
##############################      Log   Details      ##############################
#####################################################################################

if(model == "linear"){
  cat("\nWill fit Linear models:\n")
}else if(model == "glm") {
  cat("\nWill fit Generalized Linear Models:\n")
  cat("\t Family:",glmfam,"\n")
  cat("\t Link:",get(glmfam)()$link,"\n")
}
cat("\t Formula:", paste(form), "\n")
cat("\nRun details:\n")
cat("\t Chromosomes:", chr,"\n")
cat("\t N Observations in both phenotype and genotype files:", nrow(dt_ana),"\n")
if(norun == FALSE) cat("\t # of workers:", nworkers,"\n")

#Sanity check
if(nrow(dt_ana) < 10){
  cat("\nERROR: Only",nrow(dt_ana),"IDs found. Wrong ID variable specified, or wrong study IDs used.\n")
  if(mpi.comm.size() > 1) mpi.close.Rslaves()
  mpi.quit()  
}

if(norun == TRUE){
  cat("\n\nLooks good! Suggested command to run in LSF is:")
  cat("\n\tbsub -n 25 -M 9 mpirun ./ncgwas_script [your options] \n\n")
  if(mpi.comm.size() > 1) mpi.close.Rslaves()
  mpi.quit()  
}

#Check MPI correctly setup
if(mpi.comm.size() <= 1){
  cat("\nERROR: too few MPI workers; did you call the script with mpirun, and the -n bsub option?\n\n")
  mpi.quit()
}

cat("\n\nALL LOOK OK\n\n")

#####################################################################################
##############################        MPI Setup        ##############################
#####################################################################################

#Split convenience function
splitup <- function(a, n) lapply(split(a[1]:a[2], cut(a[1]:a[2], n)), range)

#Send ALL objects (because why not), and the needed libraries to worker threads
mpi.bcast.Robj2slave(all = TRUE) 
mpi.bcast.cmd({
  libs <- c("Rcpp","data.table","Matrix","MASS","RcppEigen","parallel","speedglm",
            "ncdf4", "pryr")
  lapply(libs, function(i) suppressPackageStartupMessages(require(i, ch = TRUE)))
})

#####################################################################################
#############################        Workhorse       ################################
#####################################################################################

dofun <- function(k, ch){
  memsplit <- 4
  start <- k[1]; end <- k[2]; span <- end-start+1
  #Open nc file, get SNP names and create results dataset
  nc <- nc_open(paste0(gpath,study,'-chr',ch,'-c.nc'))
  pos <- as.integer(tstrsplit(ncvar_get(nc, "SNP_Name",c(1,start), c(-1,span)),":")[[2]])
  res_part <- data.table(
    index = as.integer(start:end), 
    pos = pos, 
    b = as.numeric(NA), se = as.numeric(NA), 
    p = as.numeric(NA), key = "index")
  if(model == "glm") res_part[,conv := as.numeric(NA)]
  #Read dosages at relevant indices, and restrict to participants also in phenotype file
  for(x in splitup(k,memsplit)){
    start <- x[1]; end <- x[2]; span <- end-start+1
    newind <- res_part[,between(index,start,end)]
    res_part[newind,j := seq_along(start:end)]
    dos <- ncvar_get(nc,"Prob_AB", start=c(start,1), count=c(span, -1))[,nckeep]
    dos <- dos + 2*ncvar_get(nc,"Prob_AA", start=c(start,1), count=c(span, -1))[,nckeep]
    
    #Add allele frequency, variance and nonomissing N
    res_part[newind,c("caf_sample", "v", "n") := list(mean(dos[j,]/2, na.rm = TRUE), 
                                         var(dos[j,], na.rm = TRUE),
                                         sum(!is.nan(dos[j,]), na.rm = TRUE)), j]
    
    #Add regression results using qfit functions applied to every column of the dosage matrix,
    #wrapped with the data.table by= operator for speed. qfit_lm and qfit_glm are defined above.
    if(model == "linear"){ 
      res_part[newind & n > 0 & abs(1-caf_sample) > mincaf & v > 0 , c("b", "se") := qfit_lm(dos[j,]), j]
    }else res_part[newind & n > 0 & abs(1-caf_sample) > mincaf & v > 0 , c("b", "se","p","conv") := qfit_glm(dos[j,]), j]
    
    #Return data.table copy to avoid memory leaks
    res_part[,j:=NULL]
  }
  mem <- mem_used()
  rm(dos);gc()
  list(k = k, res = copy(res_part), mem = mem)
}

#####################################################################################
#############################      RUN THE GWAS      ################################
#####################################################################################

if(!nomerge){
  if(identical(sort(chr), 1:22)){ rname <- paste0(resdir,"/allchr_",outcome,"_",study,"_results.csv")
  }else rname <- paste0(resdir,"/partial_run_",outcome,"_",study,"_results.csv")
}
    
cat("\nStarting main loop at:",format(Sys.time(),"%H:%M:%S"),"\n")
for(i in chr){
  #misc: get #snps, make output file name, send chromosome # to workers...
  mem <- 0
  #splitup task into optimized # of chunks
  nc <- nc_open(paste0(gpath,study,'-chr',i,'-c.nc'))
  nsnp <- nc$dim$SNPs$len
  bits <- splitup(c(1,nsnp), nworkers)
  
  tag <- floor(runif(1,1,1000))
  mpi.bcast.cmd(.mpi.worker.apply, n=nworkers, tag=tag)
  mpi.bcast.Robj(list(FUN=dofun,dot.arg=list(ch=i)),rank=0)
  mpi.scatter.Robj(c(list("master"),as.list(bits)),root=0)

  #Write previous chromosome output while new one computing
  if(i != chr[1]){
    snpSum <- biganno[snpSum]
    invisible(fwrite(snpSum, rname, sep = ",", append = (!nomerge & i != chr[2])))
  }
  if(nomerge) rname <- paste0(resdir,"/Chr",i,"_",outcome,"_",study,"_results.csv")
  
  #Read giant annotation file during first CHR run
  if(i == chr[1]){
    biganno <- fread("/proj/epi/CVDGeneNas/Giant/1KG/references/cptid_rsmid_map.1000G_ALL_p1v3.merged_mach_impute.v1.txt",
                   select=c("rsid", "chr", "pos"), 
                   colClasses = list("integer" = c("pos")), showProgress = FALSE)
    biganno <- biganno[chr != "X"]
    biganno[,chr := as.integer(chr)]
    setkey(biganno, chr, pos)
  }

  #Load annotations while new data computing
  load(paste0(annodir,study,"_snp_summary_chr",i,".Rdata"))
  setDT(snpSum)
  snpSum[,chr := as.integer(i)]
  snpSum[,ref_allele := factor(Baseline_Allele, labels = c("A", "C", "G", "T", "I", "D", "R"))]
  snpSum[,eff_allele := factor(Count_Allele, labels = c("A", "C", "G", "T", "I", "D", "R"))]
  
  #Grab output from worker threads
  for (j in seq_along(bits)){
     tmp <- mpi.recv.Robj(mpi.any.source(),tag)
     src <- mpi.get.sourcetag()[1]
     snpSum[tmp$k[1]:tmp$k[2],names(tmp$res) := tmp$res]
     if(tmp$mem > mem) mem <- tmp$mem
  }
  if(model == "linear") snpSum[,p := 2*(1-pt(abs(b/se),n-nvar))]
  setkey(snpSum, chr, pos)
  
  #Last output of results
  if(i == last(chr)){
    snpSum <- biganno[snpSum]
    fwrite(snpSum, rname, sep = ",", append = !nomerge & (i != first(chr)), showProgress = FALSE)
  }
  cat("\nDone with chromosome: ",i," at: ",format(Sys.time(),"%H:%M:%S"),
      ".\n\t Size of output: ", format(object.size(snpSum), units = "MB"),
      ".\n\t Max memory used: ", floor(mem/1E6), "MB\n", sep = "")
  warnings()
}

cat("\n\nAll done. Results may be found in ",resdir,"\n\n")
invisible(mpi.close.Rslaves())
mpi.quit()

#####################################################################################
#############################          DONE          ################################
#####################################################################################
