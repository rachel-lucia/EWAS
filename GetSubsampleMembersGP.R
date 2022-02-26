setwd("M:/DATA/Methylation Data/Glyphosate_Methylation_Data_Analysis/")

set.seed(793)

samples = read.csv("SampleSheetGP.csv", header = TRUE)
samples$Sample_Name = as.character(samples$Sample_Name)
samples$order = 1:392

for (i in 1:1000){
  # subsampling
  samples_i = samples[sample(1:392, 353, replace = FALSE), ]
  samples_i = samples_i[order(samples_i$order),] # sample also does not preserve order
  write.table(t(samples_i$Sample_Name), file = "Subsample_Members_GP_1000.csv", row.names = i, append = TRUE, sep =",", col.names = FALSE)
}