library(OxyBS)
library(minfi) 
library(IlluminaHumanMethylation450kmanifest) 
library(IlluminaHumanMethylation450kanno.ilmn12.hg19) 
library(GEOquery)
#data(OxyBSSampleData)

################################
# 1. Get GEO data
################################
setwd("/storage/nmishra/CTE-CFE/")

# # # Find file names based on GEO "GSM" individual sample identifier
files <- list.files(getwd())[grep(".idat", list.files(getwd()))]
#idat_files <- list.files(file.path(destination), pattern="GSM")

pheno <- read.csv("SampleSheet.txt", header = TRUE, sep = "\t", check.names = FALSE)
# Double check that the IDAT files match from the phenotype file and directory of downloaded IDAT files
#length(intersect(pheno_base, base)) # 8 total array positions for 4 cerebellum samples

# Order the arrays such that samples match on BS and oxBS 
pheno <- pheno[naturalsort::naturalorder(pheno$title),]
pheno$pheno_base <- paste0(pheno$Sentrix_ID, "_", pheno$Sentrix_Position)
# Create 'flags' for sample treatment with BS/oxBS
BS_indices <- grep("-BS", pheno$title)
oxBS_indices <- grep("-oxBS", pheno$title)

# Name the array and slide position for BS treated samples
#arraysBS <- intersect(pheno$pheno_base[BS_indices], base)
arraysBS <- pheno$pheno_base[BS_indices]
# Name the array and slide position for oxBS treated samples
#arraysOxBS <- intersect(pheno$pheno_base[oxBS_indices], base)
arraysOxBS <- pheno$pheno_base[oxBS_indices]

dataDir <- "/storage/nmishra/CTE-CFE"
datListBS <- read.metharray(file.path(dataDir, arraysBS),verbose=TRUE)
datListOxBS <- read.metharray(file.path(dataDir, arraysOxBS),verbose=TRUE)

####### Detection of probes failed in detection P-value < 0.05 ######## 
####### In this analysis I am using all IDAT files in one run. 
####### It doesn't matter if we use all IDAT in one run or hydroxy and BS separately. We will get same detection P-value.
rgSet <- read.metharray.exp(dataDir) # Read all file. It will use SampleSheet.csv by default
detP <- detectionP(RGSet)
failed <- detP>0.05
failedProbes <-rownames(failed)[rowMeans(failed)>0.2]
#################################################################
# Create new preprocessFunNorm() function
preprocessFunnormRedGreen <- function (rgSet, nPCs = 2, sex = NULL, verbose = TRUE)
{
  #minfi:::.isRG(rgSet)
  rgSet <- updateObject(rgSet)
  if (verbose)
    cat("[preprocessFunnorm] Mapping to genome\n")
  gmSet <- mapToGenome(rgSet)
  subverbose <- max(as.integer(verbose) - 1L, 0)
  if (verbose)
    cat("[preprocessFunnorm] Quantile extraction\n")
  extractedData <- minfi:::.extractFromRGSet450k(rgSet)
  if (is.null(sex)) {
    gmSet <- addSex(gmSet, getSex(gmSet, cutoff = -3))
    sex <- rep(1L, length(gmSet$predictedSex))
    sex[gmSet$predictedSex == "F"] <- 2L
  }
  rm(rgSet)
  if (verbose)
    cat("[preprocessFunnorm] Normalization\n")
  CN <- getCN(gmSet)
  minfi:::.normalizeFunnorm450k(object = gmSet, extractedData = extractedData,
                                sex = sex, nPCs = nPCs, verbose = subverbose)
  
}

# Preprocess using functional normalization separately for both BS and oxBS arrays
rgBS <-  preprocessFunnormRedGreen(datListBS)
rgOxBS <-  preprocessFunnormRedGreen(datListOxBS)


####################################################
# 3. Define normalized signals and calculate Betas
####################################################

# Methylated signals from the BS and oxBS arrays
methBS <- assay(rgBS,"Meth")
methOxBS <- assay(rgOxBS,"Meth")
# Unmethylated signals from the BS and oxBS arrays
unmethBS <- assay(rgBS,"Unmeth")
unmethOxBS <- assay(rgOxBS,"Unmeth")

# Calculate Total Signals
signalBS <- methBS+unmethBS
signalOxBS <- methOxBS+unmethOxBS

# Calculate Beta Values
betaBS <- methBS/signalBS
betaOxBS <- methOxBS/signalOxBS

####################################################
# 4. Apply fitOxBS function to preprocessed values
####################################################

# Select the number of CpGs and Subjects to which the method will be applied 
nCpGs <- dim(unmethOxBS)[1]
nSpecimens <- dim(unmethOxBS)[2]
# 
# Create container for the OxyBS results
MethOxy <- array(NA,dim=c(nCpGs,nSpecimens,3))
dimnames(MethOxy) <- list(
  rownames(methBS)[1:nCpGs],
  colnames(methBS)[1:nSpecimens], c("C","5mC","5hmC"))

# Process results (one array at a time, slow)
#for(i in 1:nSpecimens){
#MethOxy[,i,] <-fitOxBS(betaBS[,i],betaOxBS[,i],signalBS[,i],signalOxBS[,i])
#}

# Process the results using parallelization and implementation of foreach function
# Import parallelization packages
library(foreach)
library(doParallel)
# 
# # Set-up parallel backend to use all but one available processors
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)
# 
# # Calculate the time of funciton
strt <- Sys.time()
print(strt)
# # Parallelized loop for 4 Cerebellum samples across all CpGs on the 450K array 
MethOxy <- foreach(i = 1:nSpecimens, .combine=cbind, .packages='OxyBS') %dopar% {
  fitOxBS(betaBS[,i],betaOxBS[,i],signalBS[,i],signalOxBS[,i])
}
# 
end <- Sys.time()
print(strt-end)
print(Sys.time()-strt) # Time difference of 31.42 mins on 16GB local machine
stopCluster(cl)
# 
# ####################################################
# # 5. Output - Examples and visualizations
# ####################################################
# 
# # Output is a 3D array for 5C, 5mC, and 5hmC levels
# # First specimen
MethOxy[,1,]
# # 5-hydroxymethylcytosine values for first few CpGs
head(MethOxy[,,3])
# 
# # Ranges for each cytosine modification
range(MethOxy[,,1]) # 5C
range(MethOxy[,,2]) # 5mC
range(MethOxy[,,3]) # 5hmC
# 
# # Check that results sum to one
# table(apply(MethOxy,1:2,sum))
# 
# # Some NaNs may be produced when using OxyBS
any(is.na(MethOxy)) # TRUE
# 
# # Compare OxyBS with naive approach
MethOxy0 <- MethOxy
MethOxy0[,,1] <- 1-betaBS
MethOxy0[,,2] <- betaOxBS
MethOxy0[,,3] <- betaBS-betaOxBS
# 

########
#last 6 character of colname
stringr::str_sub(colnames(betaBS), start=-6)
stringr::str_sub(colnames(betaBS), -17,-8)


tmpBS <- betaBS; tmpOxBS <- betaOxBS; tmpSignalBS <- signalBS; tmpSignalOxyBS <- signalOxBS
colnames(tmpBS) <- pheno$title[BS_indices]
colnames(tmpBS) <- gsub("-BS", "-BS-", colnames(tmpBS))
colnames(tmpOxBS) <- pheno$title[oxBS_indices]
colnames(tmpOxBS) <- gsub("-oxBS", "-BS-", colnames(tmpOxBS))
colnames(tmpSignalBS) <- pheno$title[BS_indices]
colnames(tmpSignalBS) <- gsub("-BS", "-BS-", colnames(tmpSignalBS))
colnames(tmpSignalOxyBS) <- pheno$title[oxBS_indices]
colnames(tmpSignalOxyBS) <- gsub("-oxBS", "-BS-", colnames(tmpSignalOxyBS))



temp <- ENmix::oxBS.MLE(tmpBS,tmpOxBS,tmpSignalBS,tmpSignalOxyBS)



result <- list("5C"=MethOxy[,seq(1, ncol(MethOxy), 3)], "5hC"=MethOxy[,seq(2, ncol(MethOxy), 3)], "5hmC"= MethOxy[,seq(3, ncol(MethOxy), 3)]) 
rownames(result$`5C`) <- rownames(betaBS); colnames(result$`5C`) <- colnames(betaBS)
rownames(result$`5hC`) <- rownames(betaBS); colnames(result$`5hC`) <- colnames(betaBS)
rownames(result$`5hmC`) <- rownames(betaBS); colnames(result$`5hmC`) <- colnames(betaBS)
