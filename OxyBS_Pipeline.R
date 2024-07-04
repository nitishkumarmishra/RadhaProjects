########### Here I will use OxyBS ###
########### But OxyBs is super slow. So I will use oxBS.MLE function from ENmix
### oxBS.MLE is very fast, its hmc value is almost same to OxyBS, but -ve value will be zero

library(OxyBS)
library(minfi) 
library(IlluminaHumanMethylation450kmanifest) 
library(IlluminaHumanMethylation450kanno.ilmn12.hg19) 
library(GEOquery)
data(OxyBSSampleData)

################################
# 1. Get GEO data
################################
setwd("C:/Users/nitish.mishra/Desktop")
# mainDir <- getwd()
# #subDir <- "OxyBS"
# 
# if (file.exists(subDir)){
#   setwd(file.path(mainDir))
# } else {
#   dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
#   setwd(file.path(mainDir, subDir))
#   
# }
destination <- "OxyBS"
dir.create(destination, recursive = T, showWarnings = F)
# Load series and platform data from GEO
getGEOSuppFiles("GSE63179", baseDir = "OxyBS") 

# Get Phenotype Data
gse <- getGEO("GSE63179")
pheno <- pData(gse[[1]])

# Define a directory 
destination <- "OxyBS/Data/"
dir.create(destination, recursive = T, showWarnings = F)

# Create a phenotype file for each GEO dataset for downstream analyses
write.table(pheno, "OxyBS/Data/PhenoData_GSE63179.csv", sep = ",", row.names = T, col.names = NA)

untar("OxyBS/GSE63179/GSE63179_RAW.tar", exdir = destination)

# gunzip the files in this newly untarred directory replacing all files
files <- list.files(destination)[grep(".idat.gz", list.files(destination))]
for(i in 1:length(files))
{
  gunzip(paste(destination, files[i], sep = ""), overwrite = T)
}

# Find file names based on GEO "GSM" individual sample identifier
idat_files <- list.files(file.path(destination), pattern="GSM")

# Define the basenames for each IDAT file
base <- unique(lapply(idat_files, function (x) {paste(unlist(strsplit(x, "_"))[1], unlist(strsplit(x, "_"))[2], unlist(strsplit(x, "_"))[3], sep = "_")}))

# Create phenotype basenames to match those listed with each IDAT file
pheno_base <- unique(lapply(as.character(pheno$supplementary_file), function (x) {paste(unlist(strsplit(unlist(strsplit(x, "/"))[7], "_"))[1], unlist(strsplit(x, "_"))[2], unlist(strsplit(x, "_"))[3], sep = "_")}))
pheno$pheno_base <- pheno_base

# Double check that the IDAT files match from the phenotype file and directory of downloaded IDAT files
length(intersect(pheno_base, base)) # 8 total array positions for 4 cerebellum samples

# Order the arrays such that samples match on BS and oxBS 
pheno_ordered <- pheno[order(pheno$title, decreasing=F), ]

# Create 'flags' for sample treatment with BS/oxBS
BS_indices <- grep("brain-BS", pheno_ordered$title)
oxBS_indices <- grep("oxBS", pheno_ordered$title)

# Name the array and slide position for BS treated samples
arraysBS <- intersect(pheno_ordered$pheno_base[BS_indices], base)

# Name the array and slide position for oxBS treated samples
arraysOxBS <- intersect(pheno_ordered$pheno_base[oxBS_indices], base)

dataDir <- "OxyBS/Data"
datListBS <- read.metharray(file.path(dataDir, arraysBS),verbose=TRUE)
datListOxBS <- read.metharray(file.path(dataDir, arraysOxBS),verbose=TRUE)

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


########### Here I am using ENmix::oxBS.MLE ########
########
#last 6 character of colname
#stringr::str_sub(colnames(betaBS), start=-6)
#stringr::str_sub(colnames(betaBS), -17,-8)


tmpBS <- betaBS; tmpOxBS <- betaOxBS; tmpSignalBS <- signalBS; tmpSignalOxyBS <- signalOxBS
colnames(tmpBS) <- pheno_ordered$title[grep("brain-BS", pheno_ordered$title)]
colnames(tmpOxBS) <- pheno_ordered$title[grep("brain-BS", pheno_ordered$title)]
colnames(tmpSignalBS) <- pheno_ordered$title[grep("brain-BS", pheno_ordered$title)]
colnames(tmpSignalOxyBS) <- pheno_ordered$title[grep("brain-BS", pheno_ordered$title)]

temp <- ENmix::oxBS.MLE(tmpBS,tmpOxBS,tmpSignalBS,tmpSignalOxyBS)

###############
save.image("OxyBD_pipeline.RData")
