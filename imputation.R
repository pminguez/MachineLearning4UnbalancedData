rm(list=ls())

library(irlba)
library(DMwR)
library(pca3d)
library(rgl)

setwd("/home/pablo/genetica/NeuralNetNeumo/data/")

input.file <- "data_curated.txt"
#input.file <- "input.knn.matrix.txt"
input.matrix <- read.table(input.file, sep="\t", header=T)


## Impute values on NAs using KNN
input.knn.matrix <- knnImputation(input.matrix, k=10)[c(1:8,10:18)] 


## Write matrices KNN
input.knn.matrix.file <- "input.knn.matrix.txt"

write.table(t(input.knn.matrix), file=input.knn.matrix.file,
            sep="\t", quote=F)








## Read classes
classes.file <- "classes.txt"
classes.matrix <- read.table(classes.file, sep="\t")

classes.names <- as.vector(classes.matrix[,1])
classes <- as.vector(classes.matrix[,2])
colors <- as.vector(classes.matrix[,4])

# pca.dineof <- prcomp(whole.dineof.matrix, cor=TRUE, scores=TRUE)
pca.knn <- prcomp(input.knn.matrix, cor=TRUE, scores=TRUE)

groups <- as.factor(classes)

## Print PCA dineof
# pca3d(pca.dineof, components = 1:3, col = colors, title = NULL, new = TRUE, radius=1.4,
#       axes.color = "black", bg = "white", group = groups, shape="sphere",
#       show.shadows = FALSE,show.plane = FALSE,
#       show.ellipses = F, ellipse.ci = 0.65)

# Print PCA knn
pca3d(pca.knn, components = 1:3, col = colors, title = NULL, new = TRUE, radius=1.4,
      axes.color = "black", bg = "white", group = groups, shape="sphere",
      show.shadows = FALSE,show.plane = FALSE,
      show.ellipses = F, ellipse.ci = 0.65)


