# redirect console output to text file
setwd("/pub/mcfarlar/")
output <- file("LoadFilterOut.txt")
sink(file = output, append = TRUE, type = "output", split = TRUE)
sink(file = output, append = TRUE, type = "message")

# load required packages; put session info into text file
library(minfi)
library(ChAMP)
library(ChAMPdata)
sessionInfo()

# define location of IDAT files
baseDir <- "/data/users/mcfarlar/IDAT_Files"

# read IDAT files
targets <- read.metharray.sheet(baseDir)

myRaw <- read.metharray.exp(targets = targets)

myPvals <- detectionP(myRaw,type = "m+u")
myBetas <- getBeta(myRaw)
failed <- myPvals > 0.05
myBetas[failed] <- NA

pd <- pData(myRaw)

# save RGSet object & beta matrix
save(myRaw, file = "RawRGSet.RData")
write.csv(myBetas, file = "RawBetas.csv", row.names = TRUE)

# Noob normalization
myNoob <- preprocessNoob(myRaw)
myNoobBetas <- getBeta(myNoob)

# filtering
myFilter <- champ.filter(beta=myNoobBetas, pd = pd, detP = myPvals, autoimpute = FALSE, SampleCutoff = 0.01, ProbeCutoff = 0.2, detPcut = 0.05, filterBeads = FALSE, filterNoCG = FALSE, filterSNPs = TRUE, filterMultiHit = TRUE, filterXY = TRUE, fixOutlier = TRUE, arraytype = "EPIC")

# BMIQ normalization (TII probe normalization)
# slowest step

# first, clean up workspace to free up memory
rm(myNoob,myNoobBetas,myPvals,myRaw,myBetas)
gc()

myNormed <- champ.norm(beta = myFilter$beta, resultsDir = "./BMIQ", plotBMIQ = FALSE, arraytype="EPIC", cores = 16)
write.csv(myNormed, file = "NormalizedBetas.csv", row.names = TRUE)

# save R session, commands, and filtered/normalized beta values for later
save(file="LoadFilter.RData",list = ls(all.names = TRUE), envir = .GlobalEnv)
q(save = "no", status = 0)


#### code for later
#targets$HighLow <- 0
#targets$HighLow[targets$GP > quantile(targets$GP)[3]] <- 1
# QC plot
#mdsPlot(myRaw, numPositions = 1000, sampNames = c(1:400), sampGroups = pd$HighLow)
