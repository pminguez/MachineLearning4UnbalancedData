#######################################################################################
## Pablo Minguez 14/05/2018                                                          ##
## Bioinformatics Unit - Department of Genetics & Genomics                           ##
## IIS-Fundacion Jimenez Diaz                                                        ##
## Code to perform random forest to data with unbalanced classes                     ##
## MORTALITY parameter has to be the last column and should be named MORTALITY       ##
## Everybody can download and modify the script and use it under own responsability  ##
## ################################################################################# ##

rm(list=ls())

library(randomForest)
library(rfUtilities)
library(pROC)
library(HandTill2001)
library(ggplot2)
library(irlba)
library(DMwR)

## PARAMETERS TO CUSTOMIZE
resampling <- "no"
k <- 10 ## k-cross fold validation

data.dir <- "/home/pablo/genetica/NeuralNetNeumo/data/"
results.dir <- "/home/pablo/genetica/NeuralNetNeumo/results/"
variables.to.remove <- c(8) ## Position of variables to remove from the analysis

## Name of files
input.file <- "data_curated_sorted.moredata.txt"
input.knn.matrix.file <- "input.knn.matrix.moredata.txt"
roc.curve.file <- "rf_curvasroc.png"

setwd(data.dir)

## Transpose the file to count NAs
nas.matrix <- read.table(nas.file, sep="\t", header=T)
write.table(t(nas.matrix), file=nas.t.file,sep="\t", quote=F)

## First we input missing data using KNN
input.matrix <- read.table(input.file, sep="\t", header=T, na.strings = "NULL", 
                                        stringsAsFactors = F, colClasses = "numeric")

## Sort according to MORTALITY (1s first) ## In case dinput.natrix is not sorted (1s in the first columns and 0s in the last ones)
#input.matrix <- input.matrix[order(input.matrix$MORTALITY, decreasing = TRUE),]

## Impute values on NAs using KNN

input.knn.matrix <- knnImputation(as.matrix(input.matrix), k=10)[,-variables.to.remove] 
## Variable INR is removed (41% of NAs)

## Write matrices KNN
write.table(input.knn.matrix, file=input.knn.matrix.file,sep="\t", quote=F)

matrix.file <- paste(data.dir, input.knn.matrix.file, sep="")

## Vectors that save each performance to calculate mean afterwards
auc.vect <- rep(NA,k) ## store the error in this vector (AUC)
accuracy.vect <- rep(NA, k) ## store the accuracy in this vector
sensitivity.global.vector <- rep(NA, k)
specificity.global.vector <- rep(NA, k)

## This is to build the ROC curve after the 10 fold cross-validation
## 70 is the number of rows of cv.test * 10
classes.global.vector <- rep(NA, 7*k)
probabilities.global.vector <- rep(NA, 7*k)

importance.frame <- data.frame(Parameter=character(0),mda=numeric(0),mdg=numeric(0))

seed.sum <- 100
lucky.number <- 100

for(i in 1:k){
  set.seed(i + seed.sum)
  
  matrix.frame <- read.table(matrix.file, sep="\t", header=T)
  
  ## This is the case we want to introduce the same number of samples in both classes
  if(resampling=="yes"){
    number.to.keep <- sum(matrix.frame$MORTALITY==1)
    rows <- c(1:sum(matrix.frame$MORTALITY==1), sample(sum(matrix.frame$MORTALITY==1)+1:dim(matrix.frame)[1], number.to.keep, replace=F))
    matrix.frame <- matrix.frame[rows,]
    rownames(matrix.frame) <- 1:dim(matrix.frame)[1]
  }
  
  # Define train and test sets
  samp <- sample(nrow(matrix.frame), 0.9 * nrow(matrix.frame))
  cv.train <- matrix.frame[samp, ]
  cv.test <- matrix.frame[-samp, ]
  
  ## This in principle avoids that we have all classes in the cv.test data.frame
  while(length(unique(cv.test$MORTALITY)) != length(unique(matrix.frame$MORTALITY))){
    set.seed(i + seed.sum + lucky.number)
    samp <- sample(nrow(matrix.frame), 0.9 * nrow(matrix.frame))
    cv.train <- matrix.frame[samp, ]
    cv.test <- matrix.frame[-samp, ]
    lucky.number <- lucky.number + 1
  }
  
  ## Convert response variable to factor (MORTALITY) to avoid RF does regression
  cv.train$MORTALITY <- as.character(cv.train$MORTALITY)
  cv.train$MORTALITY <- as.factor(cv.train$MORTALITY)
  cv.test$MORTALITY <- as.character(cv.test$MORTALITY)
  cv.test$MORTALITY <- as.factor(cv.test$MORTALITY)
  
  
  ## Find best mtry number
  bestmtry <- tuneRF(cv.train, cv.train$MORTALITY, ntreeTry = 1000, 
                     stepFactor = 1.2, improve= 0.01,
                     trace=F, plot=F)
  

  ## Fit model ...
  fit <- randomForest(MORTALITY~., data=cv.train, importance=TRUE, mtry=bestmtry, ntree=1000)
  
  prediction.class <- predict(fit, newdata=cv.test[,-dim(cv.test)[2]], type="class")
  prediction.prob <- predict(fit, newdata=cv.test[,-dim(cv.test)[2]], type="prob")
  
  ## Importance
  importance(fit)
  #varImpPlot(fit)
  
  for(j in 1:length(colnames(matrix.frame)[1:length(colnames(matrix.frame))-1])){
    measurement <- colnames(matrix.frame)[j]
    importance.frame <- rbind(importance.frame, data.frame(Parameter=measurement, 
                                                           mda=as.vector(importance(fit)[j,dim(importance(fit))[2]-1]), 
                                                           mdg=as.vector(importance(fit)[j,dim(importance(fit))[2]])))
  }
  
  ## Table predictions vs real
  t <- table(predictions=prediction.class, real=cv.test$MORTALITY) ## predictions vs real
  t
  
  ## AUC calculated using all classes, other methods only allow two classes
  auc.vect[i] <- auc(multcap(response = cv.test$MORTALITY,predicted = as.matrix(prediction.prob)))

  ## Accuracy
  accuracy.vect[i] <- sum(diag(t))/sum(t) ## percentage accuracy
  
  ## Calculate sensitivity & specificity
  sensitivity.vector <- rep(NA, dim(t)[1])
  specificity.vector <- rep(NA, dim(t)[1])
  
  for (z in 1:dim(t)[1]){
    tps <- t[z,z] ## True Positives for a class is when the prediction is equal to the real in that class only
    tns <- sum(t[-z,-z]) ## True Negatives are those that do not belong to the class and are not assigned to that class
    fps <- sum(t[z,-z]) 
    fns <- sum(t[-z,z])
    
    sens <- tps/(tps+fns)
    sensitivity.vector[z] <- sens
    
    spec <- tns/(tns+fps)
    specificity.vector[z] <- spec
  }
  
  specificity.global.vector[i] <- mean(specificity.vector)
  sensitivity.global.vector[i] <- mean(sensitivity.vector)
  
  # feed the vectors for ROC curve
  for(new.i in 1:length(prediction.prob[,1])){
    modified.i <- ((i - 1) * length(prediction.prob[,1])) + new.i
    probabilities.global.vector[modified.i] <- prediction.prob[new.i, 2]
    classes.global.vector[modified.i] <- as.vector(cv.test$MORTALITY[new.i])
  }
  
  predictions <- as.data.frame(prediction.prob)
  predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
  predictions$observed <- cv.test$MORTALITY
}

## ROC curve with pROC package
## Prints both curves, 2 and all classes
png(paste(results.dir, roc.curve.file, sep=""), height=7, width=7, units="in", res=300)
# Draw ROC curve . 2 curves (one normal, the other multiclass)
result.roc <- roc(factor(classes.global.vector), probabilities.global.vector) 

## Plot both
plot.roc(result.roc, col="red")
dev.off()

#### Data for two classes ####

## Now deal with other paramenters (as before adding the ROC curve)
print(paste("Average AUC:", mean(auc.vect)))
print(paste("Average Accuracy:", mean(accuracy.vect)))
print(paste("Average Sensitivity:", mean(sensitivity.global.vector)))
print(paste("Average Specificity:", mean(specificity.global.vector)))

## Plot Importance parameters
levels(importance.frame$Parameter)

## This reorders the legend
importance.reorder <- reorder(importance.frame$Parameter, importance.frame$mda, FUN=mean)
levels.reordered <- rev(attributes(importance.reorder)$levels)

importance.frame$Parameter <- factor(importance.frame$Parameter, levels = levels.reordered)

p<-ggplot(importance.frame, aes(x=reorder(Parameter, mda, FUN=mean), y=mda, color=Parameter)) +
  geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_y_continuous("Mean Decrease Accuracy", limits=c(min(importance.frame$mda), 
                                                        max(importance.frame$mda))) +
  scale_x_discrete("Parameters")
p


importance.reorder <- reorder(importance.frame$Parameter, importance.frame$mdg, FUN=mean)
levels.reordered <- rev(attributes(importance.reorder)$levels)

importance.frame$Parameter <- factor(importance.frame$Parameter, levels = levels.reordered)


p<-ggplot(importance.frame, aes(x=reorder(Parameter, mdg, FUN=mean), y=mdg, color=Parameter)) +
  geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_y_continuous("Mean Decrease Gini", limits=c(min(importance.frame$mdg), 
                                                    max(importance.frame$mdg))) +
  scale_x_discrete("Parameters")
p

