load("./GP/ComBatOut.RData")


# get genome-wide means for ComBat corrected data, overall and by genomic location
# prepare data
manifest = read.csv("Manifest_B4_Included.csv", header = TRUE, row.names = 1)
droppedProbes = read.csv("ProbestoRemove.csv",row.names = 1, header = TRUE)
ENCODE = read.csv("Manifest_Encode_Annotated.csv", header = TRUE)

manifest_merge = merge(x = manifest, y = ENCODE, by.x = "row.names", by.y = "name")
rownames(manifest_merge) = manifest_merge$Name

Ms <- myAdjMs2[!rownames(myAdjMs2) %in% rownames(droppedProbes), ]

annoMs = merge(x = Ms, y = manifest_merge, by = "row.names")

rownames(annoMs) = annoMs$Row.names

Ms2 <- annoMs[,2:393]
Bs <- (2^(Ms2))/((2^Ms2) + 1)



# calculate averages
means <- colMeans(Bs)
betas <- as.data.frame(Bs)

means2 <- aggregate(betas, list(annoMs$Relation_to_UCSC_CpG_Island), mean)
levels(means2$Group.1)[levels(means2$Group.1) == ""] = "CpGI_None"

# ENCODE Stuff
means7 <- aggregate(betas, list(annoMs$type_blood), mean)

# Gene Bodies
means5 <- aggregate(betas, list(annoMs$RML_Body), mean)
means5$Group.1 = "Gene_Body"

# Intergenic
means6 <- aggregate(betas, list(annoMs$RML_Intergenic), mean)
means6$Group.1 = "Intergenic"


# bind together write to files
meansLocation <- rbind(means2, means5, means6, means7)

write.csv(means, file = "./GP/GenomeWideMeans.csv", row.names= TRUE)
write.csv(meansLocation, file = "./GP/GenomeWideMeansbyLocation.csv", row.names = TRUE)