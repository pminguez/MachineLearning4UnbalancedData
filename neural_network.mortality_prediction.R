#######################################################################################
## Pablo Minguez 10/09/2018                                                          ##
## Bioinformatics Unit - Department of Genetics & Genomics                           ##
## IIS-Fundacion Jimenez Diaz                                                        ##
## Code to perform neural network to data with unbalanced classes                    ##
## MORTALITY parameter has to be the last column and should be named MORTALITY       ##
## Everybody can download and modify the script and use it under own responsability  ##
## ################################################################################# ##

rm(list=ls())

library(DMwR)
library(neuralnet)
library(nnet)
library(boot)
library(GGally)
library(ggplot2)
library(plyr)
library(readr)
library(pROC)

## PARAMETERS TO CUSTOMIZE
resampling <- "yes" ## Say yes if the data is unbalanced
k <- 10 ## k-cross fold validation
times.the.size <- 1 #number of times that should multiply the cases with mortality to make the train/test datasets
subsampling.times <- 10
seed.sum <- 100
lucky.number <- 101
hidden.layers <- c(2)

data.dir <- "/home/pablo/genetica/NeumoRandomForest/data/neural_network/"
setwd(data.dir)
results.dir <- "/home/pablo/genetica/NeumoRandomForest/results/"
variables.to.remove <- c(9) ## Position of variables to remove from the analysis

## Name of files
input.file <- "data_curated_sorted_moredata_age_v2.csv"
input.knn.matrix.file <- "input.knn.matrix_moredata_age.txt"

## First we input missing data using KNN
input.matrix <- read.table(input.file, sep="\t", header=T, na.strings = "NULL", 
                                        stringsAsFactors = F, colClasses = "numeric")[,-variables.to.remove] 

## Impute values on NAs using KNN  and remove any variable in vector variables.to.remove
## Should do a study of % of missing values per each paramenter)
input.knn.matrix <- knnImputation(as.matrix(input.matrix), k=6)
write.table(input.knn.matrix, file=input.knn.matrix.file,sep="\t", quote=F)
matrix.file <- paste(data.dir, input.knn.matrix.file, sep="")
matrix.frame <- read.table(matrix.file, sep="\t", header=T)

## Plot data description
## ggpairs(matrix.frame)

## Subsampling j times
number.to.keep <- sum(matrix.frame$MORTALITY==1)*times.the.size
if(resampling=="no"){
  number.to.keep <- sum(matrix.frame$MORTALITY==0)
  subsampling.times <- 1
}

## Vectors that save each performance to calculate mean afterwards
auc.vect <- rep(NA,k*subsampling.times) ## store the error in this vector (AUC)
accuracy.vect <- rep(NA, k*subsampling.times) ## store the accuracy in this vector
sensitivity.vector <- rep(NA, k*subsampling.times)
specificity.vector <- rep(NA, k*subsampling.times)
MSE.nn.vect <- rep(NA, k*subsampling.times)
RMSE.nn.vect <- rep(NA, k*subsampling.times)

## This is to build the ROC curve after the 10 fold cross-validation
classes.global.vector <- vector(mode="logical")
predictions.global.vector <- vector(mode="logical")

count <- 0

# We repeat the analysis if we go for subsampling
for (j in 1:subsampling.times){
  
  matrix.frame <- read.table(matrix.file, sep="\t", header=T)
  
  ## Make formula to feed the neural network
  all.vars <- colnames(matrix.frame)
  predictor.vars <- all.vars[!all.vars%in%"MORTALITY"]
  form <- as.formula(paste("MORTALITY ~", paste(all.vars[!all.vars %in% "MORTALITY"], collapse = "+")))
  
  first.value <- sum(matrix.frame$MORTALITY==1)+1
  last.value <- dim(matrix.frame)[1]
  rows <- c(1:sum(matrix.frame$MORTALITY==1), sample(first.value:last.value, number.to.keep, replace=F))
  matrix.frame <- matrix.frame[rows,]
  
  ## 10-fold cross validation loop
  for(i in 1:k){
    set.seed(i + seed.sum)
    count <- count + 1
    
    ## Normalize the matrix using min-max scale method
    maxs <- apply(matrix.frame, 2, max) 
    mins <- apply(matrix.frame, 2, min)
    
    scaled.data <- as.data.frame(scale(matrix.frame, center = mins, scale = maxs - mins))
    index <- sample(1:nrow(scaled.data),round(0.7*nrow(scaled.data)))
    
    ## creating training and test sets
    train.nn = scaled.data[index, ]
    test.nn = scaled.data[-index, ]
    test <- matrix.frame[-index, ]
    
    ## linear.output variable is set to FALSE,
    ## given the impact of the independent variables on the dependent variable (dividend) is assumed to be non-linear
    # nn = neuralnet(formula=form, train.nn, hidden = hidden.layers, linear.output = F, threshold=0.01)
    nn = neuralnet(formula=form, train.nn, hidden = hidden.layers, linear.output = T)
    predict_test.nn <- compute(nn, test.nn[,1:dim(test.nn)[2] - 1]) # we remove MORTALITY column for the prediction
    predict_test.nn <- (predict_test.nn$net.result * (max(matrix.frame$MORTALITY) - 
                                                        min(matrix.frame$MORTALITY))) + min(matrix.frame$MORTALITY)
    test.r <- (test.nn$MORTALITY)*(max(matrix.frame$MORTALITY)-min(matrix.frame$MORTALITY))+min(matrix.frame$MORTALITY)
    
    MSE.nn <- sum((test.r - predict_test.nn)^2)/nrow(test.nn)
    MSE.nn.vect[count] <- MSE.nn
    
    RMSE.nn = (sum((test$MORTALITY - predict_test.nn)^2) / nrow(test)) ^ 0.5
    RMSE.nn.vect[count] <- RMSE.nn
    
    results <- data.frame(actual = test.nn$MORTALITY, prediction = predict_test.nn)
    roundedresults<-sapply(results,round,digits=0)
    roundedresultsdf=data.frame(roundedresults)
    attach(roundedresultsdf)
    t <- table(actual,prediction)
    
    ## plot(test.nn$MORTALITY, predict_test.nn,col='red',main='Real vs predicted NN',pch=18,cex=0.7)
    ## abline(0,1,lwd=2)
    
    ## Accuracy
    accuracy.vect[count] <- sum(diag(t))/sum(t) ## percentage accuracy

    ## ROC element
    my.roc <- roc(actual, results$prediction, ci=T)
    ## AUC calculated using all classes, other methods only allow two classes
    auc.vect[count] <- my.roc$auc[1]
    
    ## Calculate sensitivity & specificity
    specificity.vector[count] <- mean(mean(my.roc$specificities))
    sensitivity.vector[count] <- mean(mean(my.roc$sensitivities))
    
    # feed the vectors for ROC curve
    classes.global.vector <- c(classes.global.vector, actual)
    predictions.global.vector <- c(predictions.global.vector, results$prediction)
  }
}
  
result.roc <- roc(factor(classes.global.vector), predictions.global.vector, ci=T)
plot.roc(result.roc, col="grey")

## Now deal with other paramenters (as before adding the ROC curve)
print(paste("Average AUC:", mean(auc.vect)))
print(paste("Average Accuracy:", mean(accuracy.vect)))
print(paste("Average Sensitivity:", mean(sensitivity.vector)))
print(paste("Average Specificity:", mean(specificity.vector)))
print(paste("Average MSE:", mean(MSE.nn.vect)))
print(paste("Average RMSE:", mean(RMSE.nn.vect)))


#import 'gar.fun' from Github
require(devtools)
source_url('https://gist.githubusercontent.com/fawda123/6206737/raw/d6f365c283a8cae23fb20892dc223bc5764d50c7/gar_fun.r')

predictor.vars <- all.vars[!all.vars%in%"MORTALITY"]
cols<-colorRampPalette(c('lightgreen','lightblue'))(length(predictor.vars))
gar.fun(predictor.vars, nn)

par(mar=c(3,4,1,1),family='serif')
gar.fun(predictor.vars,nn)













