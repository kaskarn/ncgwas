#!/usr/bin/env Rscript 

#### contact baldassa@email.unc.edu for questions

#Load libraries
library(MASS)
library(Matrix)
library(Rmpi)
library(optparse)
library(data.table)
library(ncdf4)
library(RcppEigen)
library(speedglm)

###################################################
################# Parsing step  ###################
###################################################

#Default values not specified here to avoid overriding --source file with default values
option_list = list(
  make_option(c("--source"), type = "character", 
		help = "R file to source before options are parsed. Allows specifying inputs in an R scripts rather than. or in addition to the command line", 
		metavar = "source file"),
	make_option(c("-p", "--pheno"), type = "character", 
		help = "Phenotype file path", metavar = "phenotype file"),
	make_option(c("-r", "--resdir"), type = "character", 
		help = "Results directory path", metavar = "results"),
	make_option(c("-g", "--gpath"), type = "character",
	  help = "Path to NCDF genetic data", metavar = "Gene data"),
	make_option(c("-s", "--study"), type = "character", 
	  help = "WHI study", metavar = "study"),	
	make_option(c("-o", "--outcome"), type = "character", 
	  help = "Outcome variable to be used in models", metavar = "Outcome"),	
	make_option(c("-f", "--form"), type = "character", 
	  help = "Right-hand side of the model, starting with ~ in R-style", metavar = "Formula"),	
	make_option(c("-m", "--model"), type = "character",
	  help = "Model type: linear or GLM", metavar = "Model"),
	make_option(c("--family"), type = "character",
	  help = "For GLMs: Distribution family of outcome (binomial, Gamma, gaussian,...) see ?family", metavar = "GLM Family"),
	make_option(c("--link"), type = "character",
	  help = "For GLMs: Link function", metavar = "GLM Link"),
	make_option(c("-i", "--idvar"), type = "character",
	  help = "Name of ID variable in phenotype file")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#Reads command-line options into like-named variables. Overrides --source for convenience
if(!is.null(opt$source)) source(opt$source)
cat("\n\nCommand line arguments:\n")
for(i in names(opt)) {
  cat("\t", i,":",opt[[i]],"\n")
  assign(i, opt[[i]])
}
#Sets default if not specified in command line or --source file
if(!exists("gpath")) gpath <- "/nas02/depts/epi/Genetic_Data_Center/whi_share/whi_1000g_fh_imp/ncdf-data/"
if(!exists("resdir")) resdir <- "ncgwas_results"
if(!exists("model")) model <- "linear"
if(model == "glm") {if(!exists(family)) family <- "binomial"; if(!exists(link)) link <- NULL}

#Exit if something important is missing
im <- c("form", "outcome", "pheno")
if(sum(is.na(im_mis <- match(im, ls()))) > 0){
  print(paste("Parameter(s): ", paste(im[is.na(im_mis)], collapse = ", "), "  are missing"))
  mpi.close.Rslaves()
  mpi.quit()
}

###################################################
############### End of parsing step ###############
###################################################


###################################################
################   Setup step   ###################
###################################################

#Turn formula input into formula object
form <- as.formula(form)

#Create directory/ies for results
dir.create(resdir,showWarnings = FALSE, recursive = TRUE)

#Load and pare down data
Epidata <- fread(pheno)
invisible(Epidata[,g:=rnorm(nrow(Epidata))])
Epidata <- na.omit(Epidata[,c("id", outcome, all.vars(form)), with = F])
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
y <- dt_ana[,get(outcome)]
nvar <- ncol(X)

#Finalize analytical datasets, for GLMs
dt_ana <- dt_ana[,c(outcome, all.vars(form)),with=FALSE]

#Generate fit function from file inputs. Workhorse functions are
#RcppEigen::fastLmPure for linear models and speedglm::speedglm 
#for generalized linear models
qfit_setup <- function(gpos, model, family = family, link = link){
  if(model == "linear"){
    function(gnow){
      ind <- which(!is.na(gnow))
      X[,gpos] <- gnow
      tm <- fastLmPure(X[ind,],y[ind])
      list(tm$coefficients[gpos], tm$se[gpos])
    }
  }else if(model == "glm"){
    if (family == "") family <- "binomial"
    function(gnow){
      ind <- which(!is.na(gnow))
      dt_ana[,g:=gnow]
      tm <- speedglm(form,data=dt_ana[ind],
                     family=do.call(family,as.list(link)),
                     set.default=list(row.chunk=2000))
      as.list(summary(tm)[gpos,-3])
    }
  }else stop("Specify model as linear or logistic")
}
qfit <- qfit_setup(match("g", all.vars(form))+1, model)

#Split function
splitup <- function(a, n) lapply(split(a[1]:a[2], cut(a[1]:a[2], n)), range)

#Send objects and libraries to worker threads
mpi.bcast.Robj2slave(all = TRUE) 
mpi.bcast.cmd({
  library(data.table); library(ncdf4)
  library(Matrix); library(MASS)
  library(RcppEigen); library(parallel)
  library(speedglm)
})

#The number of workers, -1 because we omit the master thread
nworkers <- mpi.comm.size() - 1

#Data for debugging
cat("\nUniverse size:", mpi.universe.size())
cat("\nComm size:", mpi.comm.size(), "\n")

###################################################
###########    End of setup step    ###############
###################################################

###################################################
################   LET'S ROLL   ###################
###################################################

#Start loop over chromosomes
for(i in chr){
  
  #misc: get #snps, make output file name, send chromosome # to workers...
  rname <- paste0(resdir,"Chr",i,"_",outcome,"_",study,"_results.csv")
  mpi.bcast.Robj2slave(i)
  nc <- nc_open(paste0(gpath,study,'-chr',i,'-c.nc'))
  nsnp <- nc$dim$SNPs$len
  print(paste("Starting on chromosome",i,"at:",Sys.time()))
  
  #splitup task into 10 chunks to prevent memory overrun
  if(nworkers < 20){ parts <- splitup(c(1,nsnp), ceiling(20/nworkers))
  }else parts <- list(c(1,nsnp))
  for(p in parts){
    #Splitup task into indices for workers
    bits <- splitup(p,nworkers)
    res <- mpi.parLapply(bits, function(k) {
      #Open nc file, get SNP names and create results dataset
      nc <- nc_open(paste0(gpath,study,'-chr',i,'-c.nc'))
      snp_names <- ncvar_get(nc, "SNP_Name",c(1,k[1]), c(-1,k[2]-k[1]+1))
      res_part <- data.table(
        index = as.integer(k[1]:k[2]), 
        snp = snp_names, 
        coded = ncvar_get(nc,"Allele1_Reference", k[1], k[2]-k[1]+1), 
        other = ncvar_get(nc,"Allele2_Reference", k[1], k[2]-k[1]+1),
        caf = as.numeric(NA), b = as.numeric(NA),
        se = as.numeric(NA), p = as.numeric(NA), j = seq_along(k[1]:k[2]))
      
      #Read dosages at relevant indices, and restrict to participants also in phenotype file
      p_aa <- ncvar_get(nc,"Prob_AA", start=c(k[1],1), count=c(k[2]-k[1]+1, -1))[,nckeep]
      p_ab <- ncvar_get(nc,"Prob_AB", start=c(k[1],1), count=c(k[2]-k[1]+1, -1))[,nckeep]
      dos <- p_aa*2 + p_ab
      
      #Add allele frequency, variance and nonomissing N
      res_part[,c("caf", "v", "n") := list(mean(dos[j,]/2, na.rm = TRUE), 
                                           var(dos[j,], na.rm = TRUE),
                                           sum(!is.nan(dos[j,]), na.rm = TRUE)), j]
      
      #Add regression results using qfit over each column of the dosage matrix,
      #wrapped with the data.table by= operator. qfit is defined at the start of 
      #this file 
      
      if(model == "linear"){ res_part[n > 0 & v > 0, c("b", "se") := qfit(dos[j,]), j]
      }else res_part[n > 0 & v > 0, c("b", "se", "p") := qfit(dos[j,]), j]
      
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
    if(p[1] == 1){ fwrite(res, rname, sep = ",")
    }else fwrite(res, rname, sep = ",", append = TRUE)
  }
  print(paste("Done with chromosome:",i,"at:",Sys.time()))
  print(paste("Size of chunked outputs: ", format(object.size(res), units = "MB")))
}

mpi.close.Rslaves()
mpi.quit()
