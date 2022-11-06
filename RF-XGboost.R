library(xgboost)
library(caret)
library(pROC)
library(caTools)
library(data.table)
library(readr)
library(ggplot2)

#Get args; needs: input file name, name of "case" label from metadata column, # threads, # rounds
args = (commandArgs(TRUE))

### check there are input and output files from command line
if (length(args) < 4) {
  stop("Please supply proper arguments: file name, name of 'case' label (from metadata column), # threads, # rounds.", call.=FALSE)
}

### get arguments from command line
infile = args[1]
metadata_label <- args[2]
numThreads <- args[3]
numRounds <- args[4]
#metadata_label <- "TN_Breast_Cancer"
#numThreads <- 8
#numRounds <- 2

setwd("/zgrouphome/fslg_PickettLabGroup/RandomForest/")
#setwd("/Users/bep8/Library/CloudStorage/Box-Box/Work/Proposals/2022/TN_BreastCancer/RF/All")
#infile <- "Naomi_TNBC_original_count_matrix_human.tsv"

#load in file
mydata = as.data.frame(read.table(infile, sep = " "), stringsAsFactors = FALSE)
print("Done reading in data")

#move phenotype column to first
mydata <- mydata[,c(ncol(mydata),1:(ncol(mydata)-1))]
colnames(mydata)[1] <- "Type"
#remove rows where "Type" column is NA
mydata.narf <- as.data.frame(subset(mydata, Type!="NA"), stringsAsFactors = FALSE)
mytemp <- as.data.frame(mydata.narf$Type,drop=FALSE)

#apply z-score normalization across each column of the df
mydata.narf <- as.data.frame(scale(mydata.narf[2:ncol(mydata.narf)]))
mydata.narf <- cbind(mytemp,mydata.narf)
colnames(mydata.narf)[1] <- "Type"
rm(mytemp)

#remove columns with either NA or NaN anywhere in the column
mydata.narf <- mydata.narf[ , colSums(is.na(mydata.narf)) == 0]
print("Done formatting data")

#make the test & train tables
print("Preparing forest analysis")
set.seed(111)
test_size = floor(0.2 * nrow(mydata.narf))
samp = sample(nrow(mydata.narf), test_size,replace = FALSE)
y_train = mydata.narf[-samp,1]
x_train = mydata.narf[-samp,-c(1)]
y_test = mydata.narf[samp,1]
x_test = mydata.narf[samp,-c(1)]

x_train1 <- as.matrix(x_train)
x_test1 <- as.matrix(x_test)

#convert the y's to binary values
#y_train1 <- ifelse(y_train==levels(y_train)[2],1,0)
#y_test1 <- ifelse(y_test==levels(y_test)[2],1,0)
y_train1 <- ifelse(y_train==metadata_label,1,0)
y_test1 <- ifelse(y_test==metadata_label,1,0)
#metadata_label

#convert to xgboost matrices
ex_train = xgb.DMatrix(x_train1, label = y_train1)
ex_test = xgb.DMatrix(x_test1, label = y_test1)

#train model
print("Running Forest Analysis...")
print("Please look for the round number when the AUC value stops decreasing, then re-run using that number for the 'round number' parameter (to prevent over-fitting)")
bstDMatrix <- xgboost(#data = dtrain, 
    data = ex_train,#$data,
    #label = y_train1,#$label,
    ##booster = "gbtree", #default = gbtree
    ##max.depth = 4, #default = 6
    num_parallel_tree = 10000,
    subsample = 0.5,
    #colsample_bytree = 0.5,
    ##eta = 0.3, #range 0-1. Low value is more robust to overfitting; default = 0
    ##gamma = 0, #minimum loss reduction. Larger number is more conservative; default = 0
    ##max_depth = 6, #Range 1 to infinity. Maximum depth of tree.
    nthread = numThreads, #8
    nrounds = numRounds, #2
    objective = "binary:logistic",#"reg:linear",
    eval_metric = "auc",
    set.seed(111),
    #prediction = TRUE,
    verbose = 1
  )
print("Finished Forest Analysis")

#feature importance
print("Calculating Metrics")
importance <- xgb.importance(feature_names = bstDMatrix[["feature_names"]], 
                             model = bstDMatrix)
#gain: improvement in accuracy brought by a feature
#cover: relative quantity of observations concerned by a feature
#frequency: simpler way to measure gain by counting # times feature is used in trees (shouldn't be used)
#head(importance)
write.table(importance, file = paste0(infile, "-Importance_test.tsv"), append = FALSE, row.names=FALSE, na="",col.names=TRUE, sep="\t")

#produce importance figure
pdf("Importance.pdf")
xgb.plot.importance(importance_matrix = importance,
                    top_n = 20,
                    plot = TRUE,
                    )
dev.off()

#generate & save confusion matrix
pred <- predict(bstDMatrix,ex_train)
pred <-  as.numeric(pred > 0.5)
cm <- confusionMatrix(factor(pred),factor(y_train1))

#write results to file
#first the confusion matrix
cf_outfile <- paste0(infile, "-Confusion_Matrix_Summary.txt")
write.table("Prediction\tReference", file = cf_outfile, row.names = FALSE, col.names = FALSE)
write.table(cm[["table"]], file = cf_outfile, append = TRUE, row.names=TRUE, col.names=TRUE, sep="\t")

#then save the confusion matrix statistics to the same file
Label <- as.vector(c("Accuracy","95% CI","No Information Rate","P-Value [Acc > NIR]","Kappa","Mcnemar's Test P-Value","Sensitivity","Specificity","Pos Pred Value","Neg Pred Value","Precision","Recall","F1","Prevalence","Detection Rate","Detection Prevalence","Balanced Accuracy","'Positive' Class"))
Value <- as.vector(c(cm[["overall"]][["Accuracy"]],paste0(cm[["overall"]][["AccuracyLower"]]," - ",cm[["overall"]][["AccuracyUpper"]]),cm[["overall"]][["AccuracyNull"]],cm[["overall"]][["AccuracyPValue"]],cm[["overall"]][["Kappa"]],cm[["overall"]][["McnemarPValue"]],cm[["byClass"]][["Sensitivity"]],cm[["byClass"]][["Specificity"]],cm[["byClass"]][["Pos Pred Value"]],cm[["byClass"]][["Neg Pred Value"]],cm[["byClass"]][["Precision"]],cm[["byClass"]][["Recall"]],cm[["byClass"]][["F1"]],cm[["byClass"]][["Prevalence"]],cm[["byClass"]][["Detection Rate"]],cm[["byClass"]][["Detection Prevalence"]],cm[["byClass"]][["Balanced Accuracy"]],cm[["positive"]]))
cf_df <- as.data.frame(cbind(Label,Value))
write.table(cf_df,file = cf_outfile,append=TRUE,col.names = TRUE,row.names = FALSE, sep = "\t")

#Generate ROC for trained model
pdf("ROC.pdf")
roc_test <- roc(y_train1, 
                pred, algorithm = 6,#auto 
                plot=TRUE, 
                print.auc=TRUE,
                #ci=TRUE,
                #ci.alpha=0.9,
                grid=TRUE,
                #smooth=TRUE,
                max.auc.polygon=TRUE
                #transpose = FALSE,
                #auc.polygon=TRUE
                )
dev.off()
print("Analysis Complete")

