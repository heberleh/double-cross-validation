library(rpart)
library(MASS)
library(class)
library(e1071)
library(dismo)

library("parallel")
#library("foreach")
#library("doParallel")


require(foreach)
#require(doParallel)
require(doSNOW)

data = vector()




db <- read.table("./dataset/spectral_counts_no_zeros_input.txt", header=TRUE,sep="\t")

#diminuir a dimensionalidade apenas para testes rápidos do código
#db2 <- t(db[1:500,])
db2 <-t(db)

x <- as.matrix(db2[3:nrow(db2),2:(ncol(db2))])
class(x) <- "numeric"
#colnames(x)<-1:ncol(x)
colnames(x)<-db2[2,2:ncol(db2)]
#yvector<-as.vector(unlist(db2[3:ncol(db2),1]))
y <- as.factor(db2[3:nrow(db2),1])

#x = log(x)
#scale samples with mean zero and standard deviation one
#for(i in 1:nrow(x))x[i,] = (x[i,]-mean(x[i,]))/sd(x[i,])
#scale features with mean zero and standard deviation one
#for(i in 1:ncol(x))x[,i] = (x[,i]-mean(x[,i]))/sd(x[,i])
#x = 2*atan(x/2)


model <- svm(x, y,scale=TRUE)
pred <- predict(model, x)
table(pred, y)
summary(model)












