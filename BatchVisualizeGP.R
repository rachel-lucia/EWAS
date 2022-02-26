## Batch diagnostics

load("LoadFilter.RData")
load("RawMsPvalMasked.RData")

library(ChAMP)
library(limma)
library(RColorBrewer)

myNormedM <- log2(myNormed/(1-myNormed))


# limit to just 392 samples included in study
data <- read.csv("./SampleSheetGP.csv",header = TRUE, row.names = 1)

myNormedM <- myNormedM[,colnames(myNormedM) %in% rownames(data)]
myNormedB <- myNormed[,colnames(myNormed) %in% rownames(data)]

# visualize distance between groups for technical variables
mdsPlot(myNormedB, numPositions = 10000, sampGroups = data$Batch, sampNames = data$Batch, main = "Normalized Bs, labeled by Batch")
mdsPlot(myNormedB, numPositions = 10000, sampGroups = data$Sample_Plate, sampNames = data$Sample_Plate, main = "Normalized Bs, labeled by Plate")
mdsPlot(myNormedB, numPositions = 10000, sampGroups = data$Column, sampNames = data$Column, main = "Normalized Bs, labeled by Column")

# lifestyle variables/covariates
mdsPlot(myNormedB, numPositions = 10000, sampGroups = data$EatOrganic, sampNames = data$EatOrganic, main = "Normalized Bs, labeled by EatOrganic", pal = brewer.pal(3,"RdYlBu"))
mdsPlot(myNormedB, numPositions = 10000, sampGroups = data$Alcohol, sampNames = data$Alcohol, main = "Normalized Bs, labeled by Alcohol", pal = brewer.pal(4,"RdYlBu"))
mdsPlot(myNormedB, numPositions = 10000, sampGroups = data$Smoking, sampNames = data$Smoking, main = "Normalized Bs, labeled by Smoking", pal = brewer.pal(3,"RdYlBu"))
mdsPlot(myNormedB, numPositions = 10000, sampGroups = data$BMICat, sampNames = data$BMICat, main = "Normalized Bs, labeled by BMICat", pal = brewer.pal(4,"RdYlBu"))

data <- within(data, HEI.Quart <- as.integer(cut(HEI, quantile(HEI, probs=0:4/4, na.rm = TRUE), include.lowest=TRUE)))
mdsPlot(myNormedB, numPositions = 10000, sampGroups = data$HEI.Quart, sampNames = data$HEI.Quart, main = "Normalized Bs, labeled by HEI.Quart", pal = brewer.pal(4,"RdYlBu"))

levels(data$Race) <- c("As","Bl","Hi","Ot","Uk","Wh")
mdsPlot(myNormedB, numPositions = 10000, sampGroups = data$Race, sampNames = data$Race, main = "Normalized Bs, labeled by Race", pal = brewer.pal(6,"RdYlBu"))

# glyphosate and AMPA
mdsPlot(myNormedB, numPositions = 10000, sampGroups = data$GP.Tert, sampNames = data$GP.Tert, main = "Normalized Bs, labeled by GP.Tert", pal = brewer.pal(3,"RdYlBu"))
mdsPlot(myNormedB, numPositions = 10000, sampGroups = data$AMPA.Tert, sampNames = data$AMPA.Tert, main = "Normalized Bs, labeled by AMPA.Tert", pal = brewer.pal(3,"RdYlBu"))




# visualize many variables at once
data$logGP = log(data$GP)
data$logAMPA = log(data$AMPA)

keep <- c("Sentrix_ID","Batch","Sample_Plate","Column_N","Row",
          "Age","Race","Alcohol","Smoking","BMI",
          "EatOrganic","HEI","logGP","logAMPA")
SVDdata <- data[,keep]
SVDdata$Sample_Plate = as.factor(SVDdata$Sample_Plate)
SVDdata$Sentrix_ID = as.factor(SVDdata$Sentrix_ID)

champ.SVD(beta=myNormedB, pd=SVDdata, resultsDir = "./SVD/") # takes 20-30 minutes
