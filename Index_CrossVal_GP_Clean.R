##### Data setup

setwd("M:/DATA/Methylation Data/Glyphosate_Methylation_Data_Analysis/Results_Validate")
library(glmnet)
library(pROC)

freqs = read.csv("TopHitsValidate.csv", header = TRUE, stringsAsFactors = FALSE)
keep = freqs[freqs$which == "GP" & freqs$freq >= 900, "name"]

keep450 = c("cg18722557","cg26787244","cg20540608","cg20531550","cg25597976",
            "cg26483235","cg13170005","cg01430385","cg05730283","cg00355690",
            "cg02519806","cg11062848","cg04915788","cg13310154")

keeptop5 = c("cg26483235","cg06993862","cg26787244","cg05730283","cg02519806")

betas = read.csv("NormalizedBetas_OfInterest_Validate.csv", row.names = 1, header = TRUE)
betas = betas[rownames(betas) %in% keep,]

tbetas = t(betas)
samples = read.csv("SampleSheetGP_Validate.csv", row.names = 1, header = TRUE)
samples$match = paste("X",rownames(samples), sep = "")

covars = c("Validation","GP.Raw","AMPA.Raw","GP.Raw.Tert","AMPA.Raw.Tert","GP.Above.Median")
data = merge(x = samples[,c("match",covars)], y = tbetas, by.x = "match", by.y = "row.names")
rownames(data) = data$match
data = data[, 2:ncol(data)]

data$logGP = log(data$GP.Raw)
data$logAMPA = log(data$AMPA.Raw)


##### LASSO regression

# reserve validation data for later
developData = data[data$Validation == 0,]
validateData = data[data$Validation == 1,]

# assign folds for 10-fold cross-validation
##### setup
set.seed(142)

storePredictions = as.data.frame(matrix(nrow = 0, ncol = 3))
colnames(storePredictions) = c("Sample","logGP","Prediction")

storeFitStats = as.data.frame(matrix(nrow = 0, ncol = 3))
colnames(storeFitStats) = c("RSQ","RMSE","lambda")

# repeat 10-fold cross-val 100 times to find optimal lambda
for (j in 1:100){
  developData$rand = runif(332)
  developData = developData[order(developData$rand),]
  developData$fold = c(rep(1,34),rep(2,34),rep(3,33),rep(4,33),rep(5,33),rep(6,33),rep(7,33),rep(8,33),rep(9,33),rep(10,33))
  
  for (i in 1:10){
    trainData = developData[developData$fold != i,c("logGP",keep)]
    testData = developData[developData$fold == i,c("logGP",keep)]
    
    y_train = trainData$logGP
    x_train = model.matrix(logGP ~ ., data=trainData)[,-1]
    
    fit = glmnet(x_train, y_train, alpha = 1, family = "gaussian", standardize = FALSE) # alpha = 1 means lasso, alpha = 0 means ridge
    
    cv_fit = cv.glmnet(x_train, y_train, alpha = 1, family = "gaussian", standardize = FALSE) # nfolds default = 10
    opt_lambda = cv_fit$lambda.1se # 1se (most parsimonious model), min (most optimized model)
    
    # get results from test set
    y_test = testData$logGP
    x_test = model.matrix(logGP ~ ., data=testData)[,-1]
    
    y_preds = as.numeric(predict(fit, s = opt_lambda, newx = x_test))
    #plot(y_preds ~ y_test, xlab = "log(Glyphosate)", ylab = "Prediction", main = paste("Iteration",i), pch = 19)
    
    predictions = as.data.frame(matrix(nrow = nrow(testData), ncol = 3))
    colnames(predictions) = c("Sample","logGP","Prediction")
    
    predictions$Sample = rownames(testData)
    predictions$logGP = y_test
    predictions$Prediction = y_preds
    
    storePredictions = rbind(storePredictions,predictions)
    
    # store the rsq, rmse, lambda value
    sst = sum((y_test - mean(y_test))^2)
    sse = sum((y_preds - y_test)^2)
    rsq = 1 - sse/sst
    rmse = sqrt(sse/length(y_preds))
    
    stats = c(rsq,rmse,opt_lambda)
    storeFitStats = rbind(storeFitStats,stats)
  }
  
  print(paste0("Done with ",j," cross-validations"))
  j = j + 1
  
}

colnames(storePredictions) = c("Sample","logGP","Prediction")
colnames(storeFitStats) = c("RSQ","RMSE","lambda")

# do some overall plots/stats
plot(storePredictions$Prediction ~ storePredictions$logGP, xlab = "log(Glyphosate)", ylab = "Prediction", main = "Overall", pch = 19)
mean(storeFitStats$RSQ) #0.2826 lasso; 0.1062 ridge; w/ min lambda 0.3276 lasso; 0.1901 450k only; 0.1827304 top 5 only
mean(storeFitStats$RMSE) #0.8888 lasso; 1.0002 ridge; w/ min lambda 0.8566 lasso; 0.9473 450k only; 0.952989 top 5 only
mean(storeFitStats$lambda) # 0.0005097722 lasso; w/ just 450k: 0.0006884053; w/ just top 5 0.0006021781
median(storeFitStats$lambda) # 0.0005086675 lasso; w/ just 450k: 0.0006484963; w/ just top 5 0.0006015248

final_lambda = mean(storeFitStats$lambda)



#### combine models into one model
y_develop = developData$logGP
x_develop = model.matrix(logGP ~ ., data=developData[,c("logGP",keep)])[,-1]
y_validate = validateData$logGP
x_validate = model.matrix(logGP ~ ., data=validateData[,c("logGP",keep)])[,-1]

final_fit = glmnet(x_develop, y_develop, alpha = 1, family = "gaussian", standardize = FALSE) # alpha = 1 means lasso, alpha = 0 means ridge

y_preds_validate = as.numeric(predict(final_fit, s = final_lambda, newx = x_validate))
validateData$prediction = y_preds_validate

cor.test(validateData$logGP, validateData$prediction)

y_preds_develop = as.numeric(predict(final_fit, s = final_lambda, newx = x_develop))
developData$prediction = y_preds_develop

cor.test(developData$logGP, developData$prediction)

coeffs  = predict(final_fit, type = 'coefficients', s = final_lambda)
coeffs = as.data.frame(coeffs[, 1])
colnames(coeffs) = "beta"
coeffs$probe = rownames(coeffs)

# test without outlier
cor.test(validateData[rownames(validateData) != "X203072530165_R04C01","logGP"], validateData[rownames(validateData) != "X203072530165_R04C01","prediction"])

########### few notes about model building

# performance notes from 100 iterations group
# freq > 70 corr = 0.22, p = 0.09
# freq > 80 corr = 0.25, p = 0.05

# using M-values instead of beta: generally comparable but slightly worse
# correlation WITH outlier is actually worse (0.2577, p = 0.0468) because she's even more dramatically out
# recovers when she's dropped but not fully (0.4019, p = 0.0016)

# would ridge be better than lasso with this small # probes
# nope it's far worse
# corr = 0.35 p = 0.006 (lasso), corr 0.2819 p = 0.0291 (ridge)

# 1se vs min for lambda selection
# min is better in bag but worse out of bag
# min lambda corr = 0.2617, p = 0.0434 (lasso)

# so probs best choice is the original. lasso, 1se lambda

# AMPA prediction very poor

# predict creatinine-adjusted values? marginally better prediction actually.

### test in non-batch-corrected data
unadj.betas = read.csv("NormalizedNoComBatBetas_OfInterest_Validate.csv", row.names = 1, header = TRUE)

unadj.betas = unadj.betas[rownames(unadj.betas) %in% keep,]
tunadj.betas = t(unadj.betas)

unadjData = merge(x = samples[,c("match",covars)], y = tunadj.betas, by.x = "match", by.y = "row.names")
rownames(unadjData) = unadjData$match
unadjData = unadjData[, 2:ncol(unadjData)]

unadjData$logGP = log(unadjData$GP.Raw)
unadjData$logAMPA = log(unadjData$AMPA.Raw)  

# then apply fitted model
x_unadj = model.matrix(logGP ~ ., data=unadjData[,c("logGP",keep)])[,-1]
y_preds_unadj = as.numeric(predict(final_fit, s = final_lambda, newx = x_unadj))

unadjData$prediction = y_preds_unadj

plot(x = unadjData$logGP, y = unadjData$prediction, col = "black", pch = 19)
cor.test(unadjData$logGP, unadjData$prediction)
# cor = 0.435 p < 0.0001 in ALL data
# cor = 0.252 P = 0.0523 in VALIDATION data
# developing and then testing model on non-batch-adjusted data performs OK (technically significant correlation but barely, does OK in tertile space)
# just testing on non-batch-adjusted data: does OK but not quite significant unless outlier excluded. tertiles all right esp if outlier excluded.

### extra validation w/ saliva samples.
# get saliva data
salivaSamples = samples[samples$Validation == 2,]
salivaSamples$match = paste("X",rownames(salivaSamples), sep = "")

salivaBetas = read.csv("NormalizedNoComBatBetas_OfInterest_Saliva.csv", row.names = 1, header = TRUE)
salivaBetas = salivaBetas[rownames(salivaBetas) %in% keep,]
tsalivabetas = t(salivaBetas)

salivaData = merge(x = salivaSamples[,c("match",covars)], y = tsalivabetas, by.x = "match", by.y = "row.names")
rownames(salivaData) = salivaData$match
salivaData = salivaData[, 2:ncol(salivaData)]

salivaData$logGP = log(salivaData$GP.Raw)
salivaData$logAMPA = log(salivaData$AMPA.Raw)  

# then apply fitted model
x_saliva = model.matrix(logGP ~ ., data=salivaData[,c("logGP",keep)])[,-1]
y_preds_saliva = as.numeric(predict(final_fit, s = final_lambda, newx = x_saliva))

salivaData$prediction = y_preds_saliva

plot(x = salivaData$logGP, y = salivaData$prediction, col = "red", pch = 19)
cor.test(salivaData$logGP, salivaData$prediction)
# NS

### extra validation with 8 excluded samples
# get data
excluded = samples[samples$Validation == 3,]
excluded_names = paste("X",rownames(excluded), sep = "")
t.exc.betas = t(unadj.betas[rownames(unadj.betas) %in% keep,colnames(unadj.betas) %in% excluded_names])

excludedData = merge(x = excluded[,c("match",covars)], y = t.exc.betas, by.x = "match", by.y = "row.names")
rownames(excludedData) = excludedData$match
excludedData = excludedData[, 2:ncol(excludedData)]

excludedData$logGP = log(excludedData$GP.Raw)
excludedData$logAMPA = log(excludedData$AMPA.Raw)

# then apply fitted model
x_excluded = model.matrix(logGP ~ ., data=excludedData[,c("logGP",keep)])[,-1]
y_preds_excluded = as.numeric(predict(final_fit, s = final_lambda, newx = x_excluded))

excludedData$prediction = y_preds_excluded

plot(x = excludedData$logGP, y = excludedData$prediction, col = "blue", pch = 19)
cor.test(excludedData$logGP, excludedData$prediction)
# NS


#### tertile prediction

anova = aov(prediction ~ GP.Raw.Tert, data = validateData)
summary(anova) # p = 0.01045

ROC = roc(response = validateData[validateData$GP.Raw.Tert == 1 | validateData$GP.Raw.Tert == 3,"GP.Raw.Tert"], predictor = validateData[validateData$GP.Raw.Tert == 1 | validateData$GP.Raw.Tert == 3,"prediction"], plot = FALSE, auc = TRUE, ci = TRUE)
ROC$auc
ROC$ci
# 1 vs. 3, 0.7525 95% CI: 0.5958-0.9092 (DeLong) with outlier (lasso)
# 1 vs. 2, 0.665 95% CI: 0.4942-0.8358 (DeLong) with outlier (lasso)
# 2 vs. 3, 0.590 95% CI: 0.4075-0.7725 (DeLong) with outlier (lasso)

ROC =  roc(response = validateData[,"GP.Above.Median"], predictor = validateData[,"prediction"], plot = FALSE, auc = TRUE, ci = TRUE)
ROC$auc
ROC$ci
# 0.6274 95% CI: 0.4841-0.7707 (DeLong)
# not great, Bob


#### generate figure: panel A- train, panel B - validation, panel C - boxplot, panel D - AUC

tiff("./IndexPlotHighRes.tiff", units = "mm", width = 180, height = 180, res = 600)

par(mfrow = c(2,2))
par(mar = c(4.1, 4.1, 1.1, 1.1)) # bottom left top right, default (5.1, 4.1, 4.1, 2.1)

line2 = lm(prediction ~ logGP,data = developData)
plot(x = developData$logGP, y = developData$prediction, 
     xlab = "", ylab = "Predicted Glyphosate, ng/mL", main = "", cex = 0.8, 
     pch = 19, xlim = c(-4.5, 0.5), ylim = c(-4.5, 0.5), yaxt = "n", xaxt = "n")
title(xlab = "Actual Glyphosate, ng/mL", line = 2.25)
axis(side = 2, at = log(c(0.01, 0.05, 0.15, 0.5, 1, 1.5)), labels = c(0.01, 0.05, 0.15, 0.5, 1, 1.5), las = 1, cex.axis = 0.8)
axis(side = 1,  at = log(c(0.01, 0.05, 0.15, 0.5, 1, 1.5)), labels = c(0.01, 0.05, 0.15, 0.5, 1, 1.5), cex.axis = 0.8)
abline(line2, lty = "dotted")
legend("bottomleft", legend = "Training Set (n = 332)", cex = 0.8, bty = "n", pch = 19, col = "black")
text(x = -4.5, y = 0.25, labels = "Pearson r = 0.63\np<0.0001", pos = 4, cex = 0.8)
text(x = -5.75, y = 0.6, "A", cex = 1.7, pos = 4, offset = 0, xpd = TRUE)

line = lm(prediction ~ logGP,data = validateData)
plot(x = validateData$logGP, y = validateData$prediction, 
     xlab = "", ylab = "Predicted Glyphosate, ng/mL", main = "", cex = 0.8, 
     pch = 19, col = "darkorchid3",
     xlim = c(-4.5, 0.5), ylim = c(-4.5, 0.5), yaxt = "n", xaxt = "n")
title(xlab = "Actual Glyphosate, ng/mL", line = 2.25)
axis(side = 2, at = log(c(0.01, 0.05, 0.15, 0.5, 1, 1.5)), labels = c(0.01, 0.05, 0.15, 0.5, 1, 1.5), las = 1, cex.axis = 0.8)
axis(side = 1,  at = log(c(0.01, 0.05, 0.15, 0.5, 1, 1.5)), labels = c(0.01, 0.05, 0.15, 0.5, 1, 1.5), cex.axis = 0.8)
abline(line, col = "darkorchid3", lty = "dotted")
legend("bottomleft", legend = "Validation Set (n = 60)", cex = 0.8, bty = "n", pch = 19, col = "darkorchid3")
text(x = -4.5, y = 0.25, labels = "Pearson r = 0.35\np=0.006", pos = 4, cex = 0.8)
text(x = -5.75, y = 0.6, "B", cex = 1.7, pos = 4, offset = 0, xpd = TRUE)

boxplot(prediction ~ GP.Raw.Tert, data = validateData,
        xlab = "", ylab = "Methylation Index", main = "", cex.axis = 0.8, outpch = 19, outcex = 0.6, outcol = "darkorchid3",
        boxwex = 0.4, whisklty = "solid", medcol = "white", medlwd = 2, col = "darkorchid3", las = 1)
title(xlab = "Glyphosate Tertile", line = 2.25)
text(x = -0.3, y = -0.975, "C", cex = 1.7, pos = 4, offset = 0, xpd = TRUE)

ROC = roc(response = validateData[validateData$GP.Raw.Tert == 1 | validateData$GP.Raw.Tert == 3,"GP.Raw.Tert"],
          predictor = validateData[validateData$GP.Raw.Tert == 1 | validateData$GP.Raw.Tert == 3,"prediction"],
          plot = TRUE, las = 1,  
          main = "", col = "darkorchid3", mar = c(4.1, 4.1, 1.1, 1.1), cex.axis = 0.8)
text(x = 1, y = 0.95, labels = "AUC = 0.75\n95% CI 0.60-0.91", pos = 4, cex = 0.8)
text(x = 1.25, y = 1.025, "D", cex = 1.7, pos = 4, offset = 0, xpd = TRUE)

dev.off()