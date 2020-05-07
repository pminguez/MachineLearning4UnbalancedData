rm(list=ls())

library(DMwR)
library(GGally)
library(ggplot2)
library(pROC)
library(caret)
library(pca3d)
library(RColorBrewer)
library(pheatmap)
library(factoextra)
library(reshape2)
library(tidyverse)
library(pls)
library(ggplot2)

## PARAMETERS TO CUSTOMIZE
#resampling <- "no" ## Say yes if the data is unbalanced

k <- 10 ## k-cross fold validation
times.the.size <- 1 #number of times that should multiply the cases to make the train/test datasets
subsampling.times <- 10
seed.sum <- 1
lucky.number <- 101
percentage.test <- 0.8
black.and.white <- "no"

## Takes values from command line
args <- commandArgs(trailingOnly = TRUE)

## To run on RStudio change to real values behind ##
variable.type <- args[1] ## drugs or comobidities or clinical1
tag <- args[2] ## age_range: "all","31-50","51-70","71-100","31-100"
extra.tag <- args[3] ## UEvsREST (meaning: UCI and EXITUS vs REST); EXITUSvsREST; UCIvsEXITUS
resampling <- args[4]
data.dir <- args[5]
results.dir <- args[6]
input.file <- args[7]

matrix.frame <- read.table(input.file, sep="\t", header=T, na.strings = "NULL", 
                           stringsAsFactors = F, colClasses = "numeric")

## Subsampling j times
number.to.keep <- sum(matrix.frame$OUTPUT==1)*times.the.size
if(resampling=="no"){
  number.to.keep <- sum(matrix.frame$OUTPUT==0)
  subsampling.times <- 1
}

## Vectors that save each performance to calculate mean afterwards
auc.vect <- rep(NA,k*subsampling.times) ## store the error in this vector (AUC)
accuracy.vect <- rep(NA, k*subsampling.times) ## store the accuracy in this vector
sensitivity.vector <- rep(NA, k*subsampling.times)
specificity.vector <- rep(NA, k*subsampling.times)

## Importance
importance.frame <- data.frame(Parameter=character(0),value=numeric(0))

count <- 0

## Draw the first part of the roc curve
if(black.and.white == "yes"){
  my.colors <- rep("grey", 10)
}else{
  my.colors <- c("red", "deepskyblue", "yellow", "green4", "hotpink", "orange", "blueviolet", "navy", "sienna", "darksalmon")
}

## an object to save the AUC plot
png(filename=paste(results.dir, "auc.", variable.type, ".", extra.tag, ".", tag, ".results.png", sep=""), width=500, height=480)

# We repeat the analysis if we go for subsampling
for (j in 1:subsampling.times){
  matrix.frame <- read.table(input.file, sep="\t", header=T)
  
  first.value <- sum(matrix.frame$OUTPUT==1)+1
  last.value <- dim(matrix.frame)[1]
  rows <- c(1:sum(matrix.frame$OUTPUT==1), sample(first.value:last.value, number.to.keep, replace=F))
  matrix.frame <- matrix.frame[rows,]
  
  ## We record vectors to draw a ROC curve per each resampling
  classes.local.vector <- vector(mode="logical")
  predictions.local.vector <- vector(mode="logical")
  
  ## 10-fold cross validation loop
  for(i in 1:k){
    set.seed(i + seed.sum + (j*100))
    count <- count + 1
    
    index <- sample(1:nrow(matrix.frame),round(percentage.test*nrow(matrix.frame)))
    
    ## creating training and test sets
    train.nn = matrix.frame[index, ]
    test.nn = matrix.frame[-index, ]
    
    ## To avoid 
    while(length(unique(test.nn$OUTPUT)) == 1 || length(unique(train.nn$OUTPUT)) == 1){
      set.seed(i + seed.sum + (j*100) + lucky.number)
      index <- sample(1:nrow(matrix.frame),round(0.9*nrow(matrix.frame)))
      
      ## creating training and test sets
      train.nn = matrix.frame[index, ]
      test.nn = matrix.frame[-index, ]
      
      lucky.number <- lucky.number + 1
    }
    
    ## Convert response variable to factor (MORTALITY) to avoid NN does regression
    train.nn$OUTPUT <- as.character(train.nn$OUTPUT)
    train.nn$OUTPUT <- as.factor(train.nn$OUTPUT)
    test.nn$OUTPUT <- as.character(test.nn$OUTPUT)
    test.nn$OUTPUT <- as.factor(test.nn$OUTPUT)
    
    ## Get the model and prediciton
    crtl <- trainControl(method="none",classProbs = T, summaryFunction=twoClassSummary)
    my.grid <- expand.grid(.decay = seq(from = 0.1, to = 0.5, by = 0.1), .size = c(1,2,3,4,5))
    nn.caret <- train(OUTPUT~., data=train.nn, method="nnet", trainControl=crtl, tuneGrid=my.grid)
    
    #plotnet(nn.caret)
    predict_test.nn.prob <- predict(nn.caret, test.nn, type = "prob")
    predict_test.nn.class <- predict(nn.caret, test.nn)
    
    t <- table(predictions=predict_test.nn.class, real=test.nn$OUTPUT) ## predictions vs real
    t
    
    results <- data.frame(actual = test.nn$OUTPUT, prediction = predict_test.nn.prob$"0")
    
    ## Accuracy
    accuracy.vect[count] <- sum(diag(t))/sum(t) ## percentage accuracy
    
    ## ROC element
    my.roc <- roc(predictor=predict_test.nn.prob$"0", response=results$actual)
    ## AUC calculated using all classes, other methods only allow two classes
    auc.vect[count] <- my.roc$auc[1]
    
    ## Calculate sensitivity & specificity
    specificity.vector[count] <- mean(mean(my.roc$specificities))
    sensitivity.vector[count] <- mean(mean(my.roc$sensitivities))
    
    classes.local.vector <- c(classes.local.vector, as.numeric(levels(results$actual))[results$actual])
    predictions.local.vector <- c(predictions.local.vector, results$prediction)
    
    ## Importance
    for(h in 1:length(colnames(matrix.frame)[1:length(colnames(matrix.frame))-1])){
      measurement <- colnames(matrix.frame)[h]
      importance.frame <- rbind(importance.frame, data.frame(Parameter=measurement, 
                                                             value=as.vector(varImp(nn.caret)$importance[h,1])))
    }
  }
  
  ## Generate the ROC for a re-sampled matrix
  intermediate.roc <- roc(factor(classes.local.vector), predictions.local.vector, ci=T)
  
  if(j==1){
    plot.roc(intermediate.roc, col=my.colors[j], add=F, lwd=3, cex=1.2)
  }else{plot.roc(intermediate.roc, col=my.colors[j], add=T, lwd=3, cex=1.2)}
}

while (!is.null(dev.list())){dev.off()}

## Now deal with other paramenters (as before adding the ROC curve)
print(paste("Average AUC:", mean(auc.vect)))
print(paste("Average Accuracy:", mean(accuracy.vect)))
print(paste("Average Sensitivity:", mean(sensitivity.vector)))
print(paste("Average Specificity:", mean(specificity.vector)))
