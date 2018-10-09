#######################################################################################
## Pablo Minguez 24/09/2018                                                          ##
## Bioinformatics Unit - Department of Genetics & Genomics                           ##
## IIS-Fundacion Jimenez Diaz                                                        ##
## Code to perform neural network to data with unbalanced classes                    ##
## MORTALITY parameter has to be the last column and should be named MORTALITY       ##
## Everybody can download and modify the script and use it under own responsability  ##
## ################################################################################# ##

rm(list=ls())

library(DMwR)
library(GGally)
library(ggplot2)
library(pROC)
library(caret)

## PARAMETERS TO CUSTOMIZE
resampling <- "yes" ## Say yes if the data is unbalanced
k <- 10 ## k-cross fold validation
times.the.size <- 1 #number of times that should multiply the cases with mortality to make the train/test datasets
subsampling.times <- 10
seed.sum <- 1
lucky.number <- 101
percentage.test <- 0.9
black.and.white <- "yes"

data.dir <- "/home/pablo/genetica/NeumoMachineLearning/data/neural_network/"
setwd(data.dir)
results.dir <- "/home/pablo/genetica/NeumoMachineLearning/results/"
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

## Importance
importance.frame <- data.frame(Parameter=character(0),value=numeric(0))

count <- 0

## Draw the first part of the roc curve
if(black.and.white == "yes"){
  my.colors <- rep("grey", 10)
}else{
  my.colors <- c("red", "deepskyblue", "yellow", "green4", "hotpink", "orange", "blueviolet", "navy", "sienna", "darksalmon")
}

# We repeat the analysis if we go for subsampling
for (j in 1:subsampling.times){
  
  matrix.frame <- read.table(matrix.file, sep="\t", header=T)
  
  first.value <- sum(matrix.frame$MORTALITY==1)+1
  last.value <- dim(matrix.frame)[1]
  rows <- c(1:sum(matrix.frame$MORTALITY==1), sample(first.value:last.value, number.to.keep, replace=F))
  matrix.frame <- matrix.frame[rows,]
  
  ## We record vectors to draw a ROC curve per each resampling
  classes.local.vector <- vector(mode="logical")
  predictions.local.vector <- vector(mode="logical")
  
  ## 10-fold cross validation loop
  for(i in 1:k){
    set.seed(i + seed.sum + (j*100))
    count <- count + 1
    
    ## Normalize the matrix using min-max scale method
    maxs <- apply(matrix.frame, 2, max) 
    mins <- apply(matrix.frame, 2, min)
    
    scaled.data <- as.data.frame(scale(matrix.frame, center = mins, scale = maxs - mins))
    index <- sample(1:nrow(scaled.data),round(percentage.test*nrow(scaled.data)))
    
    ## creating training and test sets
    train.nn = scaled.data[index, ]
    test.nn = scaled.data[-index, ]
    test <- matrix.frame[-index, ]
    
    ## To avoid 
    while(length(unique(test.nn$MORTALITY)) == 1){
      set.seed(i + seed.sum + (j*100) + lucky.number)
      index <- sample(1:nrow(scaled.data),round(0.9*nrow(scaled.data)))
      
      ## creating training and test sets
      train.nn = scaled.data[index, ]
      test.nn = scaled.data[-index, ]
      test <- matrix.frame[-index, ]
      
      lucky.number <- lucky.number + 1
    }
    
    ## Convert response variable to factor (MORTALITY) to avoid NN does regression
    train.nn$MORTALITY <- as.character(train.nn$MORTALITY)
    train.nn$MORTALITY <- as.factor(train.nn$MORTALITY)
    test.nn$MORTALITY <- as.character(test.nn$MORTALITY)
    test.nn$MORTALITY <- as.factor(test.nn$MORTALITY)
    
    ## Get the model and prediciton
    crtl <- trainControl(method="none",classProbs = T, summaryFunction=twoClassSummary)
    my.grid <- expand.grid(.decay = seq(from = 0.1, to = 0.5, by = 0.1), .size = c(1,2,3,4,5))
    nn.caret <- train(MORTALITY~., data=train.nn, method="nnet", trainControl=crtl, tuneGrid=my.grid)
    
    #plotnet(nn.caret)
    predict_test.nn.prob <- predict(nn.caret, test.nn, type = "prob")
    predict_test.nn.class <- predict(nn.caret, test.nn)
    
    t <- table(predictions=predict_test.nn.class, real=test.nn$MORTALITY) ## predictions vs real
    t
    
    results <- data.frame(actual = test.nn$MORTALITY, prediction = predict_test.nn.prob$"0")

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
    
    ## Importance .. creo que estará bien ...
    for(h in 1:length(colnames(matrix.frame)[1:length(colnames(matrix.frame))-1])){
      measurement <- colnames(matrix.frame)[h]
      importance.frame <- rbind(importance.frame, data.frame(Parameter=measurement, 
                                                             value=as.vector(varImp(nn.caret)$importance[h,1])))
    }
  }
  
  ## Generate the ROC for a re-sampled matrix
  intermediate.roc <- roc(factor(classes.local.vector), predictions.local.vector, ci=T)
  
  if(j==1){
    plot.roc(intermediate.roc, col=my.colors[j], add=F)
  }else{plot.roc(intermediate.roc, col=my.colors[j], add=T)}
}

## Now deal with other paramenters (as before adding the ROC curve)
print(paste("Average AUC:", mean(auc.vect)))
print(paste("Average Accuracy:", mean(accuracy.vect)))
print(paste("Average Sensitivity:", mean(sensitivity.vector)))
print(paste("Average Specificity:", mean(specificity.vector)))

## Plot Importance parameters
levels(importance.frame$Parameter)

## This reorders the legend
importance.reorder <- reorder(importance.frame$Parameter, importance.frame$value, FUN=mean)
levels.reordered <- rev(attributes(importance.reorder)$levels)

importance.frame$Parameter <- factor(importance.frame$Parameter, levels = levels.reordered)

## And print the parameters
p<-ggplot(importance.frame, aes(x=reorder(Parameter, value, FUN=mean), y=value)) +
  geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_y_continuous("Importance", limits=c(min(importance.frame$value), 
                                                        max(importance.frame$value))) +
  scale_x_discrete("Parameters")
p
