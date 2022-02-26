setwd("M:/DATA/Methylation Data/Glyphosate_Methylation_Data_Analysis/Results_CR_Covariate_Adjustment")
d = read.csv("GenomeWideMeansbyLocation.csv", header = TRUE, row.names = 1,stringsAsFactors = FALSE)

# plot by location within island
boxplot(d[,c("Overall","CpGI_None","N_Shelf","N_Shore","Island","S_Shore","S_Shelf")], ylim = c(0,1), ylab = "Genome-wide Mean Methylation")

# rollup CpG categories
d$Shore = ((d$N_Shore * 73460) + (d$S_Shore * 62656))/136116
d$IslandShore = ((d$Shore * 136116) + (d$Island * 143761))/279877
d$ShelfNone = ((d$N_Shelf * 26513) + (d$S_Shelf * 24585) + (d$CpGI_None * 418201))/469299
d$Shelf = ((d$N_Shelf * 26513) + (d$S_Shelf * 24585))/51098

d$Check = ((d$IslandShore * 279877) + (d$ShelfNone * 469299))/749176 # matches genome wide mean as expected

# same deal for encode - rollup
# counts for each
#1_Active_Promoter 10_Txn_Elongation       11_Weak_Txn      12_Repressed
#100618             50094             79469             62693
#13_Heterochrom/lo 14_Repetitive/CNV 15_Repetitive/CNV   2_Weak_Promoter
#274474               834               273             38571
#3_Poised_Promoter 4_Strong_Enhancer 5_Strong_Enhancer   6_Weak_Enhancer
#20484             26626             13369             27844
#7_Weak_Enhancer       8_Insulator  9_Txn_Transition
#26220             17563              9970

d$EncPromoter = ((d$X1_Active_Promoter * 100618) + (d$X2_Weak_Promoter * 38571) + (d$X3_Poised_Promoter * 20484))/159673
d$EncEnhancer = ((d$X4_Strong_Enhancer * 26626) + (d$X5_Strong_Enhancer * 13369) + (d$X6_Weak_Enhancer * 27844) + (d$X7_Weak_Enhancer * 26220))/94059
d$EncTranscribed = ((d$X9_Txn_Transition * 9970) + (d$X10_Txn_Elongation * 50094) + (d$X11_Weak_Txn * 79469))/139533
d$EncCNV = ((d$X14_Repetitive.CNV * 834) + (d$X15_Repetitive.CNV * 273))/ 1107

d$EncStrongEnhancer = ((d$X4_Strong_Enhancer * 26626) + (d$X5_Strong_Enhancer * 13369))/39995
d$EncWeakEnhancer = ((d$X6_Weak_Enhancer * 27844) + (d$X7_Weak_Enhancer * 26220))/54064

d$EncInactive = ((d$EncCNV * 1107) + (d$X13_Heterochrom.lo * 274474))/ 275581


# there are 74 NA probes not accounted for here
d$Check = ((d$EncPromoter * 159673) + (d$EncEnhancer * 94059) + (d$EncTranscribed * 139533) + (d$X8_Insulator * 17563) + (d$X12_Repressed * 62693) + (d$X13_Heterochrom.lo * 274474) + (d$EncCNV * 1107) )/749102 
View(d[,c("Overall","Check")])
# adds up to 749102 - should be 749176

table = aggregate(d[,30:65], by = list(d$GP.Tert), mean)
table.ampa = aggregate(d[,30:65], by = list(d$AMPA.Tert), mean)

# univariate plots by genomic context
plot(x = table$Group.1, y = table$Overall)
plot(x = table$Group.1, y = table$Shore)
plot(x = table$Group.1, y = table$Island)
plot(x = table$Group.1, y = table$IslandShore)
plot(x = table$Group.1, y = table$ShelfNone)
plot(x = table$Group.1, y = table$EncPromoter)
plot(x = table$Group.1, y = table$EncEnhancer)
plot(x = table$Group.1, y = table$EncTranscribed)
plot(x = table$Group.1, y = table$EncInactive)

# prepare data: all probes
d[d$Race == "Black", "Race"] = "Other"
d[d$Race == "Unknown", "Race"] = "Other"



means = aggregate(d[d$Race == race,"Overall"], by = list(d[d$Race == race,"GP.Tert"]), FUN = mean)
sds = aggregate(d[d$Race == race,"Overall"], by = list(d[d$Race == race,"GP.Tert"]), FUN = sd)
ses = sds$x/sqrt(table(d[d$Race == race,"GP.Tert"]))
means$ci.l = means$x - 1.96 * ses
means$ci.u = means$x + 1.96 * ses
colnames(means)[1:2] = c("GP.Tert","mean")

# plot: all probes
par(mar = c(5.1, 4.1, 3.1, 2.1)) # bottom left top right, default (5.1, 4.1, 4.1, 2.1)

plot(x = means$GP.Tert, y = means$mean, type = "l",ylim = c(0.999*min(means[,2:4]), 1.001*max(means[,2:4])), xaxt = "n", 
     lwd = 2, xlab = "", ylab = "", main = race)
axis(side = 1, at = c(1,2,3), labels = c("Tertile 1","Tertile 2","Tertile 3"), padj = 1, mgp = c(3,0.5,0))
title(xlab = "Glyphosate Tertile", ylab = expression("Mean methylation ("*beta*") value"))
title(ylab = "All probes", line = 2.25)
points(x = means$GP.Tert, y = means$mean, pch = 19, col = "black")
lines(x = means$GP.Tert, y = means$ci.l, lty = "dotted")
lines(x = means$GP.Tert, y = means$ci.u, lty = "dotted")

# dotplot
library("RColorBrewer")
library("beeswarm")
cols = brewer.pal(3,"Blues")
cols = cols[3:1]

beeswarm(d$Overall ~ d$GP.Tert, pch = 21, cex = 1, method = "swarm", priority = "random", col = "black", bg = cols, at = c(1,2,3), spacing = 0.6,corral = "random", corralWidth = 0.8, axes = FALSE, xlab = "", ylab = "")
axis(side = 1, at = c(1,2,3), labels = c("Tertile 1","Tertile 2","Tertile 3"), padj = 1, mgp = c(3,0.5,0))
axis(side = 2, at = c(seq(from = 0.61,to = 0.63, by = 0.005)), labels = c(seq(from = 0.61,to = 0.63, by = 0.005)), las = 1, mgp = c(3,0.5,0))
title(xlab = "Glyphosate Tertile", ylab = expression("Mean methylation ("*beta*") value"))
title(ylab = "All probes", line = 2.25)

# boxplot
boxplot(d$Overall ~ d$GP.Tert, pch = 19, boxwex = 0.5, medcol = "white", medlwd = 1, whisklty = 1, staplelty = 0, outcex =0.4, col = cols, at = c(1,2,3), xaxt = "n", xlab = "", ylab = "")
axis(side = 1, at = c(1,2,3), labels = c("Tertile 1","Tertile 2", "Tertile 3"),mgp = c(3,0,0), padj = 1, tck = -0.015)
title(ylab = expression("Methylation ("*beta*") value"), xlab = "Glyphosate Tertile")
title(ylab = "All probes", line = 2.25)


#### promoters
# dotplot
par(mar = c(5.1, 5.1, 3.1, 2.1)) # bottom left top right, default (5.1, 4.1, 4.1, 2.1)
beeswarm(d$EncPromoter ~ d$GP.Tert, pch = 21, cex = 1, method = "swarm", priority = "random", col = "black", bg = cols, at = c(1,2,3), spacing = 0.6,corral = "random", corralWidth = 0.8, axes = FALSE, xlab = "", ylab = "")
axis(side = 1, at = c(1,2,3), labels = c("Tertile 1","Tertile 2","Tertile 3"), padj = 1, mgp = c(3,0.5,0))
axis(side = 2, at = c(seq(from = 0.13,to = 0.15, by = 0.005)), labels = c(seq(from = 0.13,to = 0.15, by = 0.005)), las = 1, mgp = c(3,0.5,0))
title(xlab = "Glyphosate Tertile", line = 3.5, ylab = expression("Mean methylation ("*beta*") value"))
title(ylab = "Probes in promoters", line = 2.5)

# boxplot
boxplot(d$EncPromoter ~ d$GP.Tert, pch = 19, boxwex = 0.5, medcol = "white", medlwd = 1, whisklty = 1, staplelty = 0, outcex =0.4, col = cols, at = c(1,2,3), xaxt = "n", xlab = "", ylab = "")
axis(side = 1, at = c(1,2,3), labels = c("Tertile 1","Tertile 2", "Tertile 3"),mgp = c(3,0,0), padj = 1, tck = -0.015)
title(ylab = expression("Methylation ("*beta*") value"), xlab = "Glyphosate Tertile")
title(ylab = "Probes in promoters", line = 2.25)



#### islands
# dotplot
par(mar = c(5.1, 5.1, 3.1, 2.1)) # bottom left top right, default (5.1, 4.1, 4.1, 2.1)
beeswarm(d$Island ~ d$GP.Tert, pch = 21, cex = 1, method = "swarm", priority = "random", col = "black", bg = cols, at = c(1,2,3), spacing = 0.6,corral = "random", corralWidth = 0.8, axes = FALSE, xlab = "", ylab = "")
axis(side = 1, at = c(1,2,3), labels = c("Tertile 1","Tertile 2","Tertile 3"), padj = 1, mgp = c(3,0.5,0))
axis(side = 2, at = c(seq(from = 0.17,to = 0.19, by = 0.005)), labels = c(seq(from = 0.17,to = 0.19, by = 0.005)), las = 1, mgp = c(3,0.5,0))
title(xlab = "Glyphosate Tertile", line = 3.5, ylab = expression("Mean methylation ("*beta*") value"))
title(ylab = "Probes in CpG Islands", line = 2.5)

# boxplot
boxplot(d$Island ~ d$GP.Tert, pch = 19, boxwex = 0.5, medcol = "white", medlwd = 1, whisklty = 1, staplelty = 0, outcex =0.4, col = cols, at = c(1,2,3), xaxt = "n", xlab = "", ylab = "")
axis(side = 1, at = c(1,2,3), labels = c("Tertile 1","Tertile 2", "Tertile 3"),mgp = c(3,0,0), padj = 1, tck = -0.015)
title(ylab = expression("Methylation ("*beta*") value"), xlab = "Glyphosate Tertile")
title(ylab = "Probes in CpG Islands", line = 2.25)





model = lm(Overall ~ GP.Tert, data = d)
summary(model)

model = lm(EncPromoter ~ GP.Tert, data = d)
summary(model)

model = lm(Island ~ GP.Tert, data = d)
summary(model)



## all significant

model = lm(Overall ~ log(GP.ngmL), data = d)
summary(model)

model = lm(EncPromoter ~ log(GP.ngmL), data = d)
summary(model)

model = lm(Island ~ log(GP.ngmL), data = d)
summary(model)

model = lm(EncEnhancer ~ log(GP.ngmL), data = d)
summary(model)

model = lm(EncTranscribed ~ log(GP.ngmL), data = d)
summary(model)

model = lm(EncInactive ~ log(GP.ngmL), data = d)
summary(model)

model = lm(CpGI_None ~ log(GP.ngmL), data = d)
summary(model)


### NS but look kinda goofy w/ all the duplicate values
model = lm(Overall ~ log(AMPA.ngmL), data = d)
summary(model)

model = lm(EncPromoter ~ log(AMPA.ngmL), data = d)
summary(model)

model = lm(Island ~ log(AMPA.ngmL), data = d)
summary(model)

model = lm(EncEnhancer ~ log(AMPA.ngmL), data = d)
summary(model)

model = lm(EncTranscribed ~ log(AMPA.ngmL), data = d)
summary(model)

model = lm(EncInactive ~ log(AMPA.ngmL), data = d)
summary(model)

model = lm(CpGI_None ~ log(AMPA.ngmL), data = d)
summary(model)



# what about adjusting for stuff?

## function to record results of linear model for specified variable
recordModel <- function(mod, var){
  #store = as.data.frame(cbind(mod[["call"]][2],t(mod[["coefficients"]][var,c(1,2,4)]), mod$r.squared, mod$adj.r.squared))
  store = as.data.frame(matrix(nrow = 1,ncol = 6))
  
  store[1,1] = as.character(mod[["call"]][2])
  store[1,2:4] = t(mod[["coefficients"]][var,c(1,2,4)])
  store[1,5] = mod$r.squared
  store[1,6] = mod$adj.r.squared
  
  colnames(store) = c("Call",paste0(var,".Coef"),paste0(var,".SE"), paste0(var,".P"),"R2","AdjR2")
  rownames(store) = NULL
  return(store)
}

d$Race = relevel(as.factor(d$Race),"White")

## linear models: ng/mL
model1 = lm(Overall ~ log(GP.ngmL), data = d) # AIC  -3404.409
write.table(recordModel(summary(model1),"log(GP.ngmL)"), file = "GenomeWideModelsResults.csv", row.names = FALSE, sep = ",", append = TRUE)

model1 = lm(Overall ~ log(AMPA.ngmL), data = d) # AIC  -3399.111
write.table(recordModel(summary(model1),"log(AMPA.ngmL)"), file = "GenomeWideModelsResults.csv", row.names = FALSE, sep = ",", append = TRUE)

adj.min = lm(Overall ~ log(GP.ngmL) + log(AMPA.ngmL) + Age + Race + CR, data = d) # AIC -3412.723 keep
adj.tech = lm(Overall ~ log(GP.ngmL) + log(AMPA.ngmL) + Age + Race + CR + Batch + Column_N + Row, data = d) # AIC -3405.03 drop

adj.cells = lm(Overall ~ log(GP.ngmL) + log(AMPA.ngmL) + Age + Race + CR + CD8T + CD4T + NK + Bcell + Mono + Neu, data = d) # AIC -3579.871 keep
adj.bmi = lm(Overall ~ log(GP.ngmL) + log(AMPA.ngmL) + Age + Race + CR + CD8T + CD4T + NK + Bcell + Mono + Neu + BMI, data = d) # AIC -3587.233 keep
adj.behav = lm(Overall ~ log(GP.ngmL) + log(AMPA.ngmL) + Age + Race + CR + CD8T + CD4T + NK + Bcell + Mono + Neu + BMI + Smoking + Alcohol, data = d) # AIC -3574.605 drop

adj.organic = lm(Overall ~ log(GP.ngmL) + log(AMPA.ngmL) + Age + Race + CR + CD8T + CD4T + NK + Bcell + Mono + Neu + BMI + EatOrganic, data = d) # AIC -3579.68 drop
adj.diet = lm(Overall ~ log(GP.ngmL) + log(AMPA.ngmL) + Age + Race + CR + CD8T + CD4T + NK + Bcell + Mono + Neu + BMI + HEI, data = d) # AIC -3401.23 drop
adj.physical = lm(Overall ~ log(GP.ngmL) + log(AMPA.ngmL) + Age + Race + CR + CD8T + CD4T + NK + Bcell + Mono + Neu + BMI + PhysAct, data = d) # AIC -3433.46 drop

# adj.bmi is best model by forward selection
write.table(recordModel(summary(adj.bmi),"log(GP.ngmL)"), file = "GenomeWideModelsResults.csv", row.names = FALSE, sep = ",", append = TRUE)
write.table(recordModel(summary(adj.bmi),"log(AMPA.ngmL)"), file = "GenomeWideModelsResults.csv", row.names = FALSE, sep = ",", append = TRUE)

# run a last one with all dropped variables
adj.all = lm(Overall ~ log(GP.ngmL) + log(AMPA.ngmL) + Age + Race + CR + CD8T + CD4T + NK + Bcell + Mono + Neu + BMI + Smoking + Alcohol + EatOrganic + HEI + PhysAct, data = d) # AIC -3259.279

write.table(recordModel(summary(adj.all),"log(GP.ngmL)"), file = "GenomeWideModelsResults.csv", row.names = FALSE, sep = ",", append = TRUE)
write.table(recordModel(summary(adj.all),"log(AMPA.ngmL)"), file = "GenomeWideModelsResults.csv", row.names = FALSE, sep = ",", append = TRUE)



# stratify by race
adj.bmi.race = lm(Overall ~ log(GP.ngmL) + log(AMPA.ngmL) + Age + CR + CD8T + CD4T + NK + Bcell + Mono + Neu + BMI, data = d[d$Race == "Other",])
write.table(recordModel(summary(adj.bmi.race),"log(GP.ngmL)"), file = "GenomeWideModelsResults.csv", row.names = FALSE, sep = ",", append = TRUE)
write.table(recordModel(summary(adj.bmi.race),"log(AMPA.ngmL)"), file = "GenomeWideModelsResults.csv", row.names = FALSE, sep = ",", append = TRUE)




Promoters = lm(EncPromoter ~ log(GP.ngmL), data = d)
NonIslands = lm(CpGI_None ~ log(GP.ngmL), data = d)
Island = lm(Island ~ log(GP.ngmL), data = d)
Genebody = lm(Gene_Body ~ log(GP.ngmL), data = d)
Intergenic = lm(Intergenic ~ log(GP.ngmL), data = d)
IslandShore = lm(IslandShore ~ log(GP.ngmL), data = d)
ShelfNone = lm(ShelfNone ~ log(GP.ngmL), data = d)
Enhancer = lm(EncEnhancer ~ log(GP.ngmL), data = d)
Transcribed = lm(EncTranscribed ~ log(GP.ngmL), data = d)
Inactive = lm(EncInactive ~ log(GP.ngmL), data = d)

write.table(recordModel(summary(Promoters),"log(GP.ngmL)"), file = "GenomeWideModelsResults.csv", row.names = FALSE, sep = ",", append = TRUE)
write.table(recordModel(summary(NonIslands),"log(GP.ngmL)"), file = "GenomeWideModelsResults.csv", row.names = FALSE, sep = ",", append = TRUE)
write.table(recordModel(summary(Island),"log(GP.ngmL)"), file = "GenomeWideModelsResults.csv", row.names = FALSE, sep = ",", append = TRUE)
write.table(recordModel(summary(Genebody),"log(GP.ngmL)"), file = "GenomeWideModelsResults.csv", row.names = FALSE, sep = ",", append = TRUE)
write.table(recordModel(summary(Intergenic),"log(GP.ngmL)"), file = "GenomeWideModelsResults.csv", row.names = FALSE, sep = ",", append = TRUE)
write.table(recordModel(summary(IslandShore),"log(GP.ngmL)"), file = "GenomeWideModelsResults.csv", row.names = FALSE, sep = ",", append = TRUE)
write.table(recordModel(summary(ShelfNone),"log(GP.ngmL)"), file = "GenomeWideModelsResults.csv", row.names = FALSE, sep = ",", append = TRUE)
write.table(recordModel(summary(Enhancer),"log(GP.ngmL)"), file = "GenomeWideModelsResults.csv", row.names = FALSE, sep = ",", append = TRUE)
write.table(recordModel(summary(Transcribed),"log(GP.ngmL)"), file = "GenomeWideModelsResults.csv", row.names = FALSE, sep = ",", append = TRUE)
write.table(recordModel(summary(Inactive),"log(GP.ngmL)"), file = "GenomeWideModelsResults.csv", row.names = FALSE, sep = ",", append = TRUE)

## linear models: ug/G
model1 = lm(Overall ~ log(GP.ugG), data = d) # AIC -3402.4
write.table(recordModel(summary(model1),"log(GP.ugG)"), file = "GenomeWideModelsResults.csv", row.names = FALSE, sep = ",", append = TRUE)

model1 = lm(Overall ~ log(AMPA.ugG), data = d) # AIC  -3398.34
write.table(recordModel(summary(model1),"log(AMPA.ugG)"), file = "GenomeWideModelsResults.csv", row.names = FALSE, sep = ",", append = TRUE)

adj.min = lm(Overall ~ log(GP.ugG) + log(AMPA.ugG) + Age + Race, data = d) # AIC -3412.085 keep
adj.tech = lm(Overall ~ log(GP.ugG) + log(AMPA.ugG) + Age + Race + Batch + Column_N + Row, data = d) # AIC -3404.375 drop

adj.cells = lm(Overall ~ log(GP.ugG) + log(AMPA.ugG) + Age + Race + CD8T + CD4T + NK + Bcell + Mono + Neu, data = d) # AIC -3579.742 keep
adj.bmi = lm(Overall ~ log(GP.ugG) + log(AMPA.ugG) + Age + Race + CD8T + CD4T + NK + Bcell + Mono + Neu + BMI, data = d) # AIC -3587.112 keep
adj.behav = lm(Overall ~ log(GP.ugG) + log(AMPA.ugG) + Age + Race + CD8T + CD4T + NK + Bcell + Mono + Neu + BMI + Smoking + Alcohol, data = d) # AIC -3574.568 drop

adj.organic = lm(Overall ~ log(GP.ugG) + log(AMPA.ugG) + Age + Race + CD8T + CD4T + NK + Bcell + Mono + Neu + BMI + EatOrganic, data = d) # AIC -3579.565 drop
adj.diet = lm(Overall ~ log(GP.ugG) + log(AMPA.ugG) + Age + Race + CD8T + CD4T + NK + Bcell + Mono + Neu + BMI + HEI, data = d) # AIC -3400.573 drop
adj.physical = lm(Overall ~ log(GP.ugG) + log(AMPA.ugG) + Age + Race + CD8T + CD4T + NK + Bcell + Mono + Neu + BMI + PhysAct, data = d) # AIC -3433.418 drop

# adj.bmi is best model by forward selection
write.table(recordModel(summary(adj.bmi),"log(GP.ugG)"), file = "GenomeWideModelsResults.csv", row.names = FALSE, sep = ",", append = TRUE)
write.table(recordModel(summary(adj.bmi),"log(AMPA.ugG)"), file = "GenomeWideModelsResults.csv", row.names = FALSE, sep = ",", append = TRUE)

# run a last one with all dropped variables
adj.all = lm(Overall ~ log(GP.ugG) + log(AMPA.ugG) + Age + Race + CD8T + CD4T + NK + Bcell + Mono + Neu + BMI + Smoking + Alcohol + EatOrganic + HEI + PhysAct, data = d) # AIC -3258.74

write.table(recordModel(summary(adj.all),"log(GP.ugG)"), file = "GenomeWideModelsResults.csv", row.names = FALSE, sep = ",", append = TRUE)
write.table(recordModel(summary(adj.all),"log(AMPA.ugG)"), file = "GenomeWideModelsResults.csv", row.names = FALSE, sep = ",", append = TRUE)



