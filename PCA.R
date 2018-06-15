rm(list=ls())

## PCA with genes
library(pca3d)

ec <- c(1,5,9,13,17,21,25,29,33,37)
tp <- c(3,7,11,15,19,23,27,31,35,39)
hc <- c(2,6,10,14,18,22,26,30,34,38)
tx <- c(4,8,12,16,20,24,28,32,36,40)

columns <- c(ec,hc,tx)
data.dir <- "/home/pablo/genetica/VIH/data/expression_matrix/"

cels <- "TCD8"
#tag <- "lessAnnotation."
tag <- ""


#expression.matrix <- read.table(paste(data.dir, cels, ".normalized_dataset.PCA.txt", sep=""), 
#                               sep="\t", quote="", header=TRUE, row.names=NULL)[,columns+1]

expression.matrix <- read.table(paste(data.dir, cels, ".hgnc.noChrY.", tag, "merged.PCA.txt", sep=""), 
                                sep="\t", quote="", header=TRUE, row.names=1)[,columns]

pca.matrix <- t(as.matrix(sapply(expression.matrix[,1:length(columns)], as.numeric)))

# classes.names <- c(
# "TRM_226726_EC","TRM_226726_HC","TRM_226726_TP","TRM_226726_TX",
# "TRM_226727_EC","TRM_226727_HC","TRM_226727_TP","TRM_226727_TX",
# "TRM_226728_EC","TRM_226728_HC","TRM_226728_TP","TRM_226728_TX",
# "TRM_226729_EC","TRM_226729_HC","TRM_226729_TP","TRM_226729_TX",
# "TRM_226730_EC","TRM_226730_HC","TRM_226730_TP","TRM_226730_TX",
# "TRM_226731_EC","TRM_226731_HC","TRM_226731_TP","TRM_226731_TX",
# "TRM_226732_EC","TRM_226732_HC","TRM_226732_TP","TRM_226732_TX",
# "TRM_226733_EC","TRM_226733_HC","TRM_226733_TP","TRM_226733_TX",
# "TRM_226734_EC","TRM_226734_HC","TRM_226734_TP","TRM_226734_TX",
# "TRM_226735_EC","TRM_226735_HC","TRM_226735_TP","TRM_226735_TX")[columns]

classes.names <- c(
"TCD8_22303_EC","TCD8_22303_HC","TCD8_22303_TP","TCD8_22303_TX",
"TCD8_22306_EC","TCD8_22306_HC","TCD8_22306_TP","TCD8_22306_TX",
"TCD8_22307_EC","TCD8_22307_HC","TCD8_22307_TP","TCD8_22307_TX",
"TCD8_22317_EC","TCD8_22317_HC","TCD8_22317_TP","TCD8_22317_TX",
"TCD8_22345_EC","TCD8_22345_HC","TCD8_22345_TP","TCD8_22345_TX",
"TCD8_26736_EC","TCD8_26736_HC","TCD8_26736_TP","TCD8_26736_TX",
"TCD8_26737_EC","TCD8_26737_HC","TCD8_26737_TP","TCD8_26737_TX",
"TCD8_26738_EC","TCD8_26738_HC","TCD8_26738_TP","TCD8_26738_TX",
"TCD8_26739_EC","TCD8_26739_HC","TCD8_26739_TP","TCD8_26739_TX",
"TCD8_26740_EC","TCD8_26740_HC","TCD8_26740_TP","TCD8_26740_TX")[columns]





classes <- c("EC", "HC", "TP", "TX", 
             "EC", "HC", "TP", "TX",
             "EC", "HC", "TP", "TX",
             "EC", "HC", "TP", "TX", 
             "EC", "HC", "TP", "TX",
             "EC", "HC", "TP", "TX", 
             "EC", "HC", "TP", "TX", 
             "EC", "HC", "TP", "TX",
             "EC", "HC", "TP", "TX", 
             "EC", "HC", "TP", "TX")[columns]

rownames(pca.matrix) <- classes.names
colnames(pca.matrix) <- rownames(expression.matrix)

colors <- c("limegreen", "skyblue3", "indianred1", "red",
            "limegreen", "skyblue3", "indianred1", "red", 
            "limegreen", "skyblue3", "indianred1", "red", 
            "limegreen", "skyblue3", "indianred1", "red", 
            "limegreen", "skyblue3", "indianred1", "red",
            "limegreen", "skyblue3", "indianred1", "red", 
            "limegreen", "skyblue3", "indianred1", "red", 
            "limegreen", "skyblue3", "indianred1", "red", 
            "limegreen", "skyblue3", "indianred1", "red", 
            "limegreen", "skyblue3", "indianred1", "red")[columns]

groups <- as.factor(classes)
pca <- prcomp(pca.matrix, cor=TRUE, scores=TRUE)
#pca <- prcomp(pca.matrix)

pca3d(pca, components = c(1,2,3), col = colors, title = NULL, new = TRUE, radius=1,
      axes.color = "black", bg = "white", group = groups, shape="sphere",
      show.shadows = FALSE,show.plane = F,
      show.ellipses = F, ellipse.ci = 0.5)


pca2d(pca, components = c(1,2), col = colors, title = NULL, new = TRUE, radius=1.4,
      axes.color = "black", bg = "white", group = groups, shape="sphere",
      show.shadows = FALSE,show.plane = FALSE, fancy=FALSE,
      show.ellipses = T, ellipse.ci = 0.5)
#dev.off()

dist(rbind(feature.centroids[feature.centroids$Type == "barolo",]
           [-1],feature.centroids[feature.centroids$Type == "grignolino",]
           [-1]), method = "euclidean")


