###### Unified script to do multiple analyses. Adjust loop input to parallelize
# 1. Records significant probes associated with glyphosate and AMPA for each iteration and associated statistics
# 2. Records summary statistics for all probes included in the analysis for: p, q, delta-M, delta-B. Records quantiles for each: 2.5, 25, 50, 75, 97.5
# 3. Records significant regions associated with GP and AMPA for each iteration and associated statistics
# 4. Records significant probes associated with smoking status (current/former/never combined)

#### SETUP
# allow you to pass the iteration # from your bash script
args <- commandArgs(trailingOnly = TRUE)
parallel_iter = as.integer(args[1])
parallel_max = 20
parallel_width = 1000/parallel_max

# redirect console output to text file
output <- file(paste0("./GP/Subsampling/Subsampling_1000_", parallel_iter,".txt"))
sink(file = output, append = TRUE, type = "output", split = TRUE)
sink(file = output, append = TRUE, type = "message")

print(paste(Sys.time()," Beginning Setup..."))

# load previous workspace and required packages
load("./GP/ComBatOut.RData")
library(limma)
library(DMRcate)
library(flock)

#### DATA PREPARATION

# read in new sample sheet with updated variables
# setup sample sheet (substitute NA values w/ median or most frequent)
samplesRaw <- read.csv("SampleSheetGP.csv",header = TRUE, row.names = 1) 
samples <- samplesRaw
samples$Alcohol[is.na(samples$Alcohol)] <- 1
samples$Smoking[is.na(samples$Smoking)] <- 0
levels(samples$Race) <- c("Asian","Black","Hispanic","Other","Other","White")
samples$Race <- relevel(samples$Race, "White")
samples$EatOrganic[is.na(samples$EatOrganic)] <- 2
samples$HEI[is.na(samples$HEI)] <- 63.86

samples$Sample_Plate <- as.factor(samples$Sample_Plate)
samples$ColBatch <- paste(samples$Column_N,samples$Batch,sep = "_")
samples$ColBatch[samples$Batch == "Pilot"] <- "Pilot"

samples$logGP = log(samples$GP.Raw)
samples$logAMPA = log(samples$AMPA.Raw)



# ComBat-adjusted Ms: myAdjMs2 (all)
# ComBat-adjusted Bs: myAdjBs2 (all)
myAdjBs2 <- (2^(myAdjMs2))/((2^myAdjMs2) + 1)
myAdjBs2[myAdjBs2 = 1] <- 0.999999

# verify that order is preserved
print("rownames check:")
sum(rownames(samples) == colnames(myAdjMs2))

# read in sample members
members = read.csv("./GP/Subsampling/Subsample_Members_GP_1000.csv", header = FALSE, row.names = 1, stringsAsFactors = FALSE)

#### DATA ANALYSIS
print(paste(Sys.time()," Beginning Analysis..."))

arg = TRUE # for later
if (parallel_iter == 1){arg = TRUE} else {arg = FALSE}

for (i in (1 + parallel_width*(parallel_iter - 1)):(parallel_width*parallel_iter)){
  # trim dataset to 90% subsample for this iteration
  samples_i = samples[rownames(samples) %in% members[i,], ]
  
  Ms_i <- myAdjMs2[, colnames(myAdjMs2) %in% rownames(samples_i)]
  Bs_i <- myAdjBs2[, colnames(myAdjBs2) %in% rownames(samples_i)]
  
  # verify that order is preserved
  if (sum(rownames(samples_i) == colnames(Ms_i)) == 353) {
    print(paste("Now running iteration",i))
  } else {
    print(paste("Rownames check failed at iteration",i))
    break
  }
  
  # run main probe-level model
  myModel_i <- model.matrix(~ Race + BMI + Age + Batch + Column_N + Row + CD8T + CD4T + NK + Bcell + Mono + Neu + Alcohol + Smoking + EatOrganic + HEI + logGP + logAMPA + CR, data=samples_i)
  
  fit <- lmFit(Ms_i, myModel_i)
  fit <- eBayes(fit)
  
  betafit <- lmFit(Bs_i, myModel_i)
  betafit <- eBayes(betafit)
  
  # 1. Records significant probes associated with GP + AMPA for each iteration and associated statistics
  
  DMPs_i <- topTable(fit, coef="logGP", num=Inf, adjust.method = "BH", p.value = 0.05, confint = TRUE)
  DMPs_i$iter = i
  write.table(DMPs_i, file = "./GP/Subsampling/DMPs_GP_M_1000.csv", row.names = TRUE, append = TRUE, sep =",", col.names = arg)
  
  DMPs_AMPA_i <- topTable(fit, coef="logAMPA", num=Inf, adjust.method = "BH", p.value = 0.05, confint = TRUE)
  DMPs_AMPA_i$iter = i
  write.table(DMPs_AMPA_i, file = "./GP/Subsampling/DMPs_AMPA_M_1000.csv", row.names = TRUE, append = TRUE, sep =",", col.names = arg)
  
  # 2. Records summary statistics for all probes included in the analysis for: p, q, delta-M, delta-B. Compute quantiles later.
  # GP first
  allMDMPs_i <- topTable(fit, coef="logGP", num=Inf, adjust.method = "BH", confint = TRUE)
  allBetaDMPs_i <- topTable(betafit, coef="logGP", num=Inf, adjust.method = "BH", confint = TRUE)
  
  allMDMPs_i = allMDMPs_i[order(rownames(allMDMPs_i)), ]
  allBetaDMPs_i = allBetaDMPs_i[order(rownames(allBetaDMPs_i)), ]
  
  pvals = t(allMDMPs_i$P.Value)
  qvals = t(allMDMPs_i$adj.P.Val)
  deltaMs = t(allMDMPs_i$logFC)
  deltaBs = t(allBetaDMPs_i$logFC)
  
  colnames(pvals) = rownames(allMDMPs_i)
  colnames(qvals) = rownames(allMDMPs_i)
  colnames(deltaMs) = rownames(allMDMPs_i)
  colnames(deltaBs) = rownames(allBetaDMPs_i)
  
  # lock the files to avoid strange overwrites
  lock.p = flock::lock("./GP/Subsampling/GP_Pvals.csv")
  lock.q = flock::lock("./GP/Subsampling/GP_Qvals.csv")
  lock.M = flock::lock("./GP/Subsampling/GP_deltaMs.csv")
  lock.B = flock::lock("./GP/Subsampling/GP_deltaBs.csv")
  
  write.table(pvals, file = "./GP/Subsampling/GP_Pvals.csv", row.names = i, append = TRUE, sep =",", col.names = arg)
  write.table(qvals, file = "./GP/Subsampling/GP_Qvals.csv", row.names = i, append = TRUE, sep =",", col.names = arg)
  write.table(deltaMs, file = "./GP/Subsampling/GP_deltaMs.csv", row.names = i, append = TRUE, sep =",", col.names = arg)
  write.table(deltaBs, file = "./GP/Subsampling/GP_deltaBs.csv", row.names = i, append = TRUE, sep =",", col.names = arg)
  
  # unlock the files
  flock::unlock(lock.p)
  flock::unlock(lock.q)
  flock::unlock(lock.M)
  flock::unlock(lock.B)
  
  # then AMPA
  allMDMPs_i <- topTable(fit, coef="logAMPA", num=Inf, adjust.method = "BH", confint = TRUE)
  allBetaDMPs_i <- topTable(betafit, coef="logAMPA", num=Inf, adjust.method = "BH", confint = TRUE)
  
  allMDMPs_i = allMDMPs_i[order(rownames(allMDMPs_i)), ]
  allBetaDMPs_i = allBetaDMPs_i[order(rownames(allBetaDMPs_i)), ]
  
  pvals = t(allMDMPs_i$P.Value)
  qvals = t(allMDMPs_i$adj.P.Val)
  deltaMs = t(allMDMPs_i$logFC)
  deltaBs = t(allBetaDMPs_i$logFC)
  
  colnames(pvals) = rownames(allMDMPs_i)
  colnames(qvals) = rownames(allMDMPs_i)
  colnames(deltaMs) = rownames(allMDMPs_i)
  colnames(deltaBs) = rownames(allBetaDMPs_i)
  
  # lock the files to avoid strange overwrites
  lock.p = flock::lock("./GP/Subsampling/AMPA_Pvals.csv")
  lock.q = flock::lock("./GP/Subsampling/AMPA_Qvals.csv")
  lock.M = flock::lock("./GP/Subsampling/AMPA_deltaMs.csv")
  lock.B = flock::lock("./GP/Subsampling/AMPA_deltaBs.csv")
  
  write.table(pvals, file = "./GP/Subsampling/AMPA_Pvals.csv", row.names = i, append = TRUE, sep =",", col.names = arg)
  write.table(qvals, file = "./GP/Subsampling/AMPA_Qvals.csv", row.names = i, append = TRUE, sep =",", col.names = arg)
  write.table(deltaMs, file = "./GP/Subsampling/AMPA_deltaMs.csv", row.names = i, append = TRUE, sep =",", col.names = arg)
  write.table(deltaBs, file = "./GP/Subsampling/AMPA_deltaBs.csv", row.names = i, append = TRUE, sep =",", col.names = arg)
  
  # unlock the files
  flock::unlock(lock.p)
  flock::unlock(lock.q)
  flock::unlock(lock.M)
  flock::unlock(lock.B)
  
  # 3. Records significant regions associated with GP + AMPA for each iteration and associated statistics
  # GP
  myannotation <- cpg.annotate(datatype = "array", Ms_i, what = "M", arraytype = "EPIC", analysis.type = "differential", design = myModel_i, coef = "logGP")
  DMRs_i <- dmrcate(myannotation, lambda=1000, C=2)
  DMR_ranges_i <- extractRanges(DMRs_i, genome = "hg19")
  DMR_ranges_i$iter = i
  
  write.table(DMR_ranges_i, file = "./GP/Subsampling/DMRs_GP_M_1000.csv", row.names = TRUE, append = TRUE, sep =",", col.names = arg)
  
  #AMPA
  myannotationAMPA <- cpg.annotate(datatype = "array", Ms_i, what = "M", arraytype = "EPIC", analysis.type = "differential", design = myModel_i, coef = "logAMPA")
  DMRs_AMPA_i <- dmrcate(myannotationAMPA, lambda=1000, C=2)
  DMR_ranges_AMPA_i <- extractRanges(DMRs_AMPA_i, genome = "hg19")
  DMR_ranges_AMPA_i$iter = i
  
  write.table(DMR_ranges_AMPA_i, file = "./GP/Subsampling/DMRs_AMPA_M_1000.csv", row.names = TRUE, append = TRUE, sep =",", col.names = arg)
  
  # 4. Records significant probes associated with smoking status (current/former/never combined)
  
  SmokDMPs_i <- topTable(fit, coef="Smoking", num=Inf, adjust.method = "BH", p.value = 0.05, confint = TRUE)
  SmokDMPs_i$iter = i
  write.table(SmokDMPs_i, file = "./GP/Subsampling/DMPs_Smoking_GP_M_1000.csv", row.names = TRUE, append = TRUE, sep =",", col.names = arg)
  # note: will have to go back later and extract summary statistics for all smoking-associated probes. 
  
  #### code for 2nd round - extract summary statistics for smoking-associated probes
  #SmokDMPs_i <- topTable(fit, coef="Smoking", num=Inf, adjust.method = "BH", confint = TRUE)
  #SmokBetaDMPs_i <- topTable(betafit, coef="Smoking", num=Inf, adjust.method = "BH", confint = TRUE)
  
  #probes = read.csv("./Density393/Subsampling/SmokingProbesWithFrequencies.csv", header = TRUE)
  #probes = as.character(probes$name)
  
  #sigSmokDMPs_i = SmokDMPs_i[probes, ]
  #sigSmokBetaDMPs_i = SmokBetaDMPs_i[probes, ]
  
  #sigSmokDMPs_i = sigSmokDMPs_i[order(rownames(sigSmokDMPs_i)), ]
  #sigSmokBetaDMPs_i = sigSmokBetaDMPs_i[order(rownames(sigSmokBetaDMPs_i)), ]
  
  #Smokpvals = t(sigSmokDMPs_i$P.Value)
  #Smokqvals = t(sigSmokDMPs_i$adj.P.Val)
  #SmokdeltaMs = t(sigSmokDMPs_i$logFC)
  #SmokdeltaBs = t(sigSmokBetaDMPs_i$logFC)
  
  #colnames(Smokpvals) = rownames(sigSmokDMPs_i)
  #colnames(Smokqvals) = rownames(sigSmokDMPs_i)
  #colnames(SmokdeltaMs) = rownames(sigSmokDMPs_i)
  #colnames(SmokdeltaBs) = rownames(sigSmokBetaDMPs_i)
  
  #write.table(Smokpvals, file = "./Density393/Subsampling/Smoking_Pvals.csv", row.names = i, append = TRUE, sep =",", col.names = arg)
  #write.table(Smokqvals, file = "./Density393/Subsampling/Smoking_Qvals.csv", row.names = i, append = TRUE, sep =",", col.names = arg)
  #write.table(SmokdeltaMs, file = "./Density393/Subsampling/Smoking_deltaMs.csv", row.names = i, append = TRUE, sep =",", col.names = arg)
  #write.table(SmokdeltaBs, file = "./Density393/Subsampling/Smoking_deltaBs.csv", row.names = i, append = TRUE, sep =",", col.names = arg)
  
  arg = FALSE # write colnames only for first iteration
  
  write.table(i, file = "./GP/Subsampling/Counter.csv", row.names = TRUE, append = TRUE, sep =",", col.names = arg)
}

print(paste(Sys.time(), i, "Iterations Complete."))

#### WRAP UP
warnings()
sessionInfo()
q(save = "no", status = 0)