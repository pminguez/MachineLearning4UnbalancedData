rm(list=ls())

## PCA with genes
library(pca3d)

data.dir <- "/home/pablo/genetica/NeumoRandomForest/data/"
#input.file <- "data_curated_sorted.moredata_age.txt"
input.file <- "input.knn.matrix.moredata_age.txt"

variables.to.remove <- c(9) ## Position of variables to remove from the analysis
setwd(data.dir)

## First we input missing data using KNN
#input.matrix <- read.table(input.file, sep="\t", header=T, na.strings = "NULL", 
#                           stringsAsFactors = F, colClasses = "numeric")[,-variables.to.remove]

input.matrix <- read.table(input.file, sep="\t", na.strings = "NULL", header=T)[,-variables.to.remove]

deads <- 1:31
alives <- 32:1966

columns <- c(deads, alives)

#pca.matrix <- t(as.matrix(sapply(expression.matrix[,1:length(columns)], as.numeric)))
pca.matrix <- as.matrix(sapply(input.matrix, as.numeric))

classes.names <- c(rep("dead", times=length(deads)), rep("alive",times=length(alives)))
classes <- c(rep("dead", times=length(deads)), rep("alive",times=length(alives)))
colors <- c(rep("red", times=length(deads)), rep("darkgreen",times=length(alives)))


#colnames(pca.matrix) <- rownames(pca.matrix)

groups <- as.factor(classes)
pca <- prcomp(pca.matrix, cor=TRUE, scores=TRUE)
#pca <- prcomp(pca.matrix)

pca3d(pca, components = c(1,2,3), col = colors, title = NULL, new = TRUE, radius=0.3,
      axes.color = "black", bg = "white", group = groups, shape="sphere",
      show.shadows = FALSE,show.plane = F,
      show.ellipses = F, ellipse.ci = 0.5)


pca2d(pca, components = c(1,2), col = colors, title = NULL, new = TRUE, radius=0.7,
      axes.color = "black", bg = "white", group = groups, shape="sphere",
      show.shadows = FALSE,show.plane = FALSE, fancy=FALSE,
      show.ellipses = T, ellipse.ci = 0.9)
dev.off()

dist(rbind(feature.centroids[feature.centroids$Type == "barolo",]
           [-1],feature.centroids[feature.centroids$Type == "grignolino",]
           [-1]), method = "euclidean")


