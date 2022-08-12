library(ggplot2)
library(corrplot)
library(reshape2)
library(ggthemes)
library(dplyr)
library(randomForest)
library(ROCR)
library(readxl)
library(tidyverse)

#options(expression = 1000000)

#load in file
setwd("~/Library/CloudStorage/Box-Box/Work/Lab/SARS-CoV2/SARS2-severity")
#setwd("/fslhome/bep8/fsl_groups/fslg_PickettLabGroup/RandomForest")
#load(file = "TCGAgeneMaster_RFinput.rda")
#setwd("~/Desktop/RF")
#infile <- "VAccine H+H Final March 25 numbers (David Edits)-bep.txt"
infile <- "count_matrix_human.excell.xlsx"
mydata <- read_excel(path = infile)

##concatenate row1 w/ colnames
mydata1 <- colnames(mydata)
mydata1 <- c("SRR_ID",mydata1)
mydata1 <- gsub("-","_",mydata1)
mydata1 <- as.character(mydata1)
colnames(mydata) <- mydata1
#mydata <- mydata[-c(1,2),]
#mydata <- mydata[-1,]
rm(mydata1)
mydata <- mydata[,c(ncol(mydata),1:(ncol(mydata)-1))] #move phenotype column to first

colnames(mydata)[1] <- "Type"
mydata.narf <- as.data.frame(subset(mydata, Type!="NA"))
mydata.narf$"44256" <- NULL #number in colname
mydata.narf$"44257" <- NULL #number in colname
#mydata.narf[,1366] <- NULL #breaks rf
#mydata.narf$"BMP15" <- NULL #NaN
mydata.narf$SRR_ID <- NULL
#mydata.narf1 <- t(na.omit(t(mydata.narf)))
#mydata.narf1 <- as.data.frame(mydata.narf1)

mydata.narf$Type <- gsub("low","l",mydata.narf$Type)
mydata.narf$Type <- gsub("high","h",mydata.narf$Type)
mytemp <- as.vector(mydata.narf$Type)

mydata.narf <- as.data.frame(scale(mydata.narf[2:ncol(mydata.narf)]))
mydata.narf <- cbind(mytemp,mydata.narf)
colnames(mydata.narf)[1] <- "Type"

#remove columns with either NA or NaN anywhere in the column
mydata.narf <- mydata.narf[ , colSums(is.na(mydata.narf)) == 0]
mydata.narf <- mydata.narf[ , colSums(is.nan(mydata.narf)) == 0]

#colnames(mydata) <- gsub("[?|.| |-|/|'|(|)|’|,|+]","_",colnames(mydata))#remove bad characters
##sapply(mydata.na, class)
#mydata.narf <- na.roughfix(mydata.na) #impute values using column medians, where original values were "NA"

#mydata.narf2 <- mydata.narf
#mydata.narf <- mydata.narf2

#mydata.narf <- mydata.narf[,c(1,1350:1367)]

set.seed(111)
test_size = floor(0.3 * nrow(mydata.narf))
samp = sample(nrow(mydata.narf), test_size,replace = FALSE)
y_train = mydata.narf[-samp,1]
x_train = mydata.narf[-samp,-c(1)]
y_test = mydata.narf[samp,1]
x_test = mydata.narf[samp,-c(1)]

#convert labels to categorical
#y_train = as.factor(y_train)
#y_test = as.factorfactor(y_test)

#Create training set and testing set
train = cbind(y_train,x_train)
test = cbind(y_test,x_test)

colnames(train)[1] = "Type"
colnames(test)[1] = "Type"

train <- as.data.frame(na.omit(train))
test <- as.data.frame(na.omit(test))

numTry <- sqrt(ncol(mydata.narf))

#mydata.narf$Type <- as.character(mydata.narf$Type)
model_1 = randomForest(Type~., data = train, ntree = 10000, importance = TRUE, do.trace = 10, proximity = TRUE, mtry = numTry)#10000
save(model_1, file = "RF-model.rda")
print(model_1)
pred_1 = predict(model_1, x_test)
summary_table <- table(y_test, pred_1)
print("Summary Table: ")
print(summary_table)
accuracy_m1 = mean(y_test == pred_1)
paste0("Accuracy: ",accuracy_m1)
varImpPlot(model_1)
importance = importance(model_1)
varImportance = data.frame(Features = row.names(importance),
                           Mean_Decrease_Gini =round(importance[, "MeanDecreaseGini"],2))
write.table(varImportance, file = paste0("varImportance_",numTry,"-Importance_test.tsv"), append = FALSE, row.names=TRUE, na="",col.names=TRUE, sep="\t")
#varImportance = data.frame(Variables = row.names(importance),
#                           Importance =round(importance[, "MeanDecreaseAccuracy"],2))
rankImportance=varImportance%>%mutate(Rank=paste("#",dense_rank(desc(Mean_Decrease_Gini))))

ggplot(rankImportance,aes(x=reorder(Features,Mean_Decrease_Gini),
                          y=Mean_Decrease_Gini,fill=Mean_Decrease_Gini))+ 
  geom_bar(stat="Identity") + 
  geom_text(aes(x = Features, y = 0.5, label = Rank),
            hjust=0, vjust=0.55, size = 4, colour = "white") +
  labs(x = "Features") +
  coord_flip() + 
  theme_classic()
outfile_name <- paste0("Gini_Bars-",numTry,"test.pdf")
ggsave(outfile_name)

#Generate AUC data and ROC curve
predictions=as.vector(model_1$votes[,2])
pred=prediction(predictions,train$Type)

perf_AUC=performance(pred,"auc") #Calculate the AUC value
AUC=perf_AUC@y.values[[1]]
paste0("AUC: ", AUC)

perf_ROC=performance(pred,"tpr","fpr") #plot the actual ROC curve
plot(perf_ROC, main="ROC plot")
text(0.5,0.5,paste("AUC = ",format(AUC, digits=5, scientific=FALSE)))
