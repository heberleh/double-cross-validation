#Installing required packages
list.of.packages <- c("MASS", "class", "e1071", "dismo", "caret", "rpart")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library("rpart")
library("MASS")
library("class")
library("e1071")
library("dismo")
library("caret")


# amostras de cada classe geradas aleatoriamente:
# 1 a 6 -> 3
# 7 a 12 -> 11
# 13 a 18 -> 13

# N selecionado como N do double cross
N = as.integer(read.csv("./results/svm-rfe/double_selected_N_svmrfe.txt", header=FALSE)[1,1])

# Esclha entre "tunar" ou não os parâmetros do SVM para cada valor de N, para cada ranking.
TUNE = FALSE
tuned = NULL # iniciando variável de tune


# leitura de treino
db <- read.table("./dataset/current/train.txt", header=TRUE,sep="\t")
db_test <- read.table("./dataset/current/independent_test.txt", header=TRUE,sep="\t")

# transposta, colunas serão proteínas
db2 <-t(db)
db2_test <-t(db_test)

# descarta header
matrix <- as.matrix(db2[3:nrow(db2),2:(ncol(db2))])
matrix_test <- as.matrix(db2_test[3:nrow(db2_test),2:(ncol(db2_test))])

class(matrix) <- 'numeric'
class(matrix_test) <- 'numeric'

aux <- matrix
nonzeros <-  (colSums(abs(aux)>0.000000))

# save constant proteins from the train set
zeros_index = which(as.logical(nonzeros == FALSE))
zeros_names = db[zeros_index,2]



#nomes das colunas (proteínas) e classes das amostras
colnames(matrix)<-db2[2,2:ncol(db2)]
colnames(matrix_test)<-db2_test[2,2:ncol(db2_test)]


#aux <- matrix
#nonzeros <-  (colSums(abs(aux))>0.000)

# save constant proteins from the train set
#zeros_index = which(as.logical(nonzeros == FALSE))+1  #+1 because of 1 line of header
#zeros_names = db[zeros_index,2]

# remove proteins wich were constant-zero in the used train set of double-cross
#matrix <- matrix[,nonzeros]

#scale *features* with mean zero and standard deviation one
#for(i in 1:ncol(matrix))matrix[,i] = (matrix[,i]-mean(matrix[,i]))/sd(matrix[,i])



y <- as.factor(db2[3:nrow(db2),1])
ytest <- as.factor(db2_test[3:nrow(db2_test),1])


# define train and test sets
x <- matrix
xtest <- matrix_test


# definindo os valores de classes
#r = 5
#y <- as.factor(c(rep('normal',r),rep('carcinoma',r),rep('melanoma',r)))
#r=1
#ytest <- as.factor(c(rep('normal',r),rep('carcinoma',r),rep('melanoma',r)))


# indexes of ranking list from double cross validation
rankdb <- read.csv("./results/svm-rfe/rank_index_svmrfe.csv", header=FALSE)
ranking_index = as.matrix(rankdb)


#remove zero-expressed proteins
for (protein_zero in zeros_index){
    ranking_index = ranking_index[ranking_index!=protein_zero]
}



predN = c()
cat("Ranking\n")
rank_accuracy <- {}
range = 2:length(ranking_index)

for(nfeatures in range){
  cat("\nRanking  ")
  cat(nfeatures)
  cat("\n\n")
  # escolhe os melhores parâmetros para o svm do conjunto atual (considerando as proteínas selecionadas atualmente neste loop)
  if (TUNE){
    tuned <- tune.svm(x=x[,ranking_index[1:nfeatures]], y=y, gamma = 10^(-6:-1), cost = 10^(1:2))
  }else{
    best = list(parameters=list(2))
    tuned = list(best=best)
    tuned$best.parameters[1] = 10^(-4)
    tuned$best.parameters[2] = 10
  }
  # constroi o modelo svm
  svmModel = svm(x=x[,ranking_index[1:nfeatures]], y=y, cost = tuned$best.parameters[2], gamma=tuned$best.parameters[1], cachesize=100,  scale=T, type="C-classification", kernel="linear")
  # classifica o conjunto de teste baseando-se no modelo criado  
  pred <- predict(svmModel, xtest[,ranking_index[1:nfeatures]])
  # se a predição corresponde à predição do N escolhido pelo Double, armazene para executar o caret
  if (nfeatures == N){
    predN <- pred
  }
 exp_pred = factor(ytest, levels=levels(pred))

  # verifica quantas amostras foram classificadas corretamente
  accuracy <- (length(which(as.logical(pred == exp_pred)))/length(ytest))
  rank_accuracy<-rbind(rank_accuracy,accuracy)
}
cat("\n\n")




cat("gniknaR\n")
unlike_rank_accuracy <- {}
max = length(ranking_index)
range = 1:(max-1)
for(i in range){
  cat("\ngniknaR  ")
  cat(i)
  cat("\n\n")
  start = max-i
  # escolhe os melhores parâmetros para o svm do conjunto atual (considerando as proteínas selecionadas atualmente neste loop)
  if(TUNE){
    tuned <- tune.svm(x=x[,ranking_index[start:max]], y=y, gamma = 10^(-6:-1), cost = 10^(1:2))
  }else{
    best = list(parameters=list(2))
    tuned = list(best=best)
    tuned$best.parameters[1] = 10^(-4)
    tuned$best.parameters[2] = 10
  }
  # constroi o modelo svm
  svmModel = svm(x=x[,ranking_index[start:max]], y=y, cost = tuned$best.parameters[2], gamma=tuned$best.parameters[1], cachesize=100,  scale=T, type="C-classification", kernel="linear")
  # classifica o conjunto de teste baseando-se no modelo criado
  pred <- predict(svmModel, xtest[,ranking_index[start:max]])
  exp_pred = factor(ytest, levels=levels(pred))

  # verifica quantas amostras foram classificadas corretamente
  accuracy <- (length(which(as.logical(pred == exp_pred)))/length(ytest))
  unlike_rank_accuracy<-rbind(unlike_rank_accuracy,accuracy)   
}



cat("\n\n")
cat("Random\n")
## generate a random ordering
set.seed(1) ## make reproducible here, but not if generating many random samples
random_rank <- sample(ncol(x))

#remove zero-expressed proteins
for (protein_zero in zeros_index){
    random_rank = random_rank[random_rank!=protein_zero]
}


random_rank_accuracy <- {}
range = 2:length(ranking_index)

for(nfeatures in range){
  cat("\nRandom  ")
  cat(nfeatures)
  cat("\n\n")
  # escolhe os melhores parâmetros para o svm do conjunto atual (considerando as proteínas selecionadas atualmente neste loop)
  if(TUNE){
    tuned <- tune.svm(x=x[,random_rank[1:nfeatures]], y=y, gamma = 10^(-6:-1), cost = 10^(1:2))
  }else{
    best = list(parameters=list(2))
    tuned = list(best=best)
    tuned$best.parameters[1] = 10^(-4)
    tuned$best.parameters[2] = 10
  }
  # constroi o modelo svm
  svmModel = svm(x=x[,random_rank[1:nfeatures]], y=y, cost = tuned$best.parameters[2], gamma=tuned$best.parameters[1], cachesize=100,  scale=T, type="C-classification", kernel="linear")
  # classifica o conjunto de teste baseando-se no modelo criado
  pred <- predict(svmModel, xtest[,random_rank[1:nfeatures]])
  exp_pred = factor(ytest, levels=levels(pred))

  # verifica quantas amostras foram classificadas corretamente
  accuracy <- (length(which(as.logical(pred == exp_pred)))/length(ytest))
  random_rank_accuracy<-rbind(random_rank_accuracy,accuracy)  
}

n_range= 2:length(ranking_index)
pdf("./results/svm-rfe/independent/ranking_svmrfe_indepentend_test_validation_N_values.pdf")
plot(n_range, rank_accuracy)
dev.off()

pdf("./results/svm-rfe/independent/gniknar_svmrfe_indepentend_test_validation_N_values.pdf")
plot(n_range, unlike_rank_accuracy)
dev.off()

pdf("./results/svm-rfe/independent/random_ranking_svmrfe_indepentend_test_validation_N_values.pdf")
plot(n_range, random_rank_accuracy)
dev.off()

# salva a matriz contendo os valores de N=1..max e suas respectivas taxas de acertos de todos os rankings
results = cbind(n_range,rank_accuracy,unlike_rank_accuracy,random_rank_accuracy)
colnames(results) <- c("N","Rank","knaR","random rank")
write.matrix(results, file = "./results/svm-rfe/independent/rankings_scores_svmrfe_independent_test.csv", sep = ",")


ytest = factor(ytest, levels=levels(predN))

sink("./results/svm-rfe/independent/caret_svmrfe_independent_test.txt")
confusionMatrix(predN, ytest, positive=NULL)
sink()



