library(rpart)
library(MASS)
library(class)
library(e1071)
library(dismo)
library(caret)
library(pamr)

# amostras de cada classe geradas aleatoriamente:
# 1 a 6 -> 3
# 7 a 12 -> 11
# 13 a 18 -> 13

# N selecionado como N do double cross
N = as.integer(read.csv("./results/nsc/double_selected_N_nsc.txt", header=FALSE)[1,1])

# Esclha entre "tunar" ou não os parâmetros do SVM para cada valor de N, para cada ranking.
TUNE = TRUE
tuned = NULL # iniciando variável de tune


# leitura de treino
db <- read.table("./dataset/current/train.txt", header=TRUE,sep="\t")

# leitura do teste independente
db_test <- read.table("./dataset/current/independent_test.txt", header=TRUE,sep="\t")


# transposta, colunas serão proteínas
db2 <-t(db)
db2_test <- t(db_test)

# descarta header
matrix <- as.matrix(db2[3:nrow(db2),2:(ncol(db2))])
matrix_test <- as.matrix(db2_test[3:nrow(db2_test),2:(ncol(db2_test))])

class(matrix) <- 'numeric'
class(matrix_test) <- 'numeric'


#aux <- matrix[-test_index,]
#nonzeros <-  (colSums(abs(aux)>0.000))

# save constant proteins from the train set
#zeros_index = which(as.logical(nonzeros == FALSE))+1  #+1 because of 1 line of header
#zeros_names = db[zeros_index,2]

# remove proteins wich were constant-zero in the used train set of double-cross
#matrix <- matrix[,nonzeros]


#scale *features* with mean zero and standard deviation one
#for(i in 1:ncol(matrix))matrix[,i] = (matrix[,i]-mean(matrix[,i]))/sd(matrix[,i])

#nomes das colunas (proteínas) e classes das amostras
colnames(matrix)<-db2[2,2:ncol(db2)]
colnames(matrix_test)<-db2_test[2,2:ncol(db2_test)]



y <- as.factor(db2[3:nrow(db2),1])
ytest <- as.factor(db2_test[3:nrow(db2_test),1])



# define train and test sets
x <- matrix
xtest <- matrix_test


# indexes of ranking list from double cross validation
rankdb <- read.csv("./results/nsc/rank_index_nsc.csv", header=FALSE)
ranking_index = as.matrix(rankdb) 


# criando variável data para o NSC
traindata = list()
testdata = list()
traindata$x = t(x)
traindata$y = y
testdata$x = t(xtest)
testdata$y = ytest

thresholds_number = 90

predN = c()
cat("Ranking\n")
rank_accuracy <- {}
range = 2:nrow(traindata$x)
#range = 1550:nrow(traindata$x)
#range = 2:30 #FOR DEBUG
for(nfeatures in range){
  cat("\nRanking  ")
  cat(nfeatures)
  cat("\n\n")
  model = NULL
  localtraindata <- traindata
  localtraindata$x <- localtraindata$x[ranking_index[1:nfeatures],]
  if (TUNE){
    #treinamento
    tuned <- pamr.train(localtraindata, scale.sd=TRUE)
    scales <- pamr.adaptthresh(tuned)
    model <- pamr.train(localtraindata, threshold.scale=scales, n.threshold=thresholds_number, scale.sd=TRUE)
  }else{
    model <- pamr.train(localtraindata, n.threshold=thresholds_number, scale.sd=TRUE)
  }
  
  # classifica o conjunto de teste baseando-se no modelo criado  
  pred = pamr.predict(model, testdata$x[ranking_index[1:nfeatures],], type = "class", threshold=0.000000000000000000000000000000000000000000000000000000000000)

  # se a predição corresponde à predição do N escolhido pelo Double, armazene para executar o caret
  if (nfeatures == N){
    predN <- pred
  }
  # verifica quantas amostras foram classificadas corretamente
  accuracy <- (length(which(as.logical(pred == testdata$y)))/length(testdata$y))
  cat("rank accuracy: ")
  cat(accuracy)
  cat("\n")
  rank_accuracy<-rbind(rank_accuracy, accuracy)
}
cat("\n\n")



cat("gniknaR\n")
unlike_rank_accuracy <- {}
max = ncol(x) #x = t(traindata$x)
#max=1567 #FOR DEBUG

aux <- as.matrix(traindata$x)
class(aux) <- "numeric"
zeros <- (rowSums(abs(aux)) == 0.0000000)
zeros_index <- which(zeros==TRUE) #128 zeros-rows
accuracy <- 0.0
#range <-start:max-1  #range de todas as últimas linhas = zeros... + 1 linha não zerada. O NSC não roda com apenas 1 variável.
# 1568:1696
aux_range = 1:(length(zeros_index))
for (i in aux_range){
   unlike_rank_accuracy<-rbind(unlike_rank_accuracy,accuracy)
}
#define o novo max. corresponde ao index maximo das linhas de não zeros.-1
# max = 1569
# start = max - i
#  start: max
# 1569 - 1568 : 1569  -> 1:1569. Ok, na posição 1569 tem linha nao zerada.
# 1569-1-  >          1568:1569. Ok, tem pelo menos 2 variáveis não zeradas.
max <-  max - length(zeros_index)
range = 1:(max-1)
#range = 1499:(max-1)
for(i in range){
  cat("\ngniknaR  ")
  cat(i)
  cat("\n\n")
  start = max-i
  model = NULL
  localtraindata <- traindata
  localtraindata$x <- localtraindata$x[ranking_index[start:max],]
  aux <- as.matrix(localtraindata$x) 
 
  if (TUNE){
    #treinamento
    tuned <- pamr.train(localtraindata, scale.sd=TRUE)
    scales <- pamr.adaptthresh(tuned)
    model <- pamr.train(localtraindata, threshold.scale=scales,n.threshold=thresholds_number, scale.sd=TRUE)
  }else{
    model <- pamr.train(localtraindata, n.threshold=thresholds_number, scale.sd=TRUE)
  }
  
  # classifica o conjunto de teste baseando-se no modelo criado
  pred = pamr.predict(model, testdata$x[ranking_index[start:max],], type = "class", threshold=0.000000000000000000000000000000000000000000000000000000000000)
  
  # verifica quantas amostras foram classificadas corretamente
  accuracy <- (length(which(as.logical(pred == testdata$y)))/length(testdata$y))
  unlike_rank_accuracy<-rbind(unlike_rank_accuracy,accuracy)
}





cat("\n\n")
cat("Random\n")
## generate a random ordering
set.seed(1) ## make reproducible here, but not if generating many random samples
random_rank <- sample(ncol(x))   #x = t(traindata$x)
random_rank_accuracy <- {}
range =  2:nrow(traindata$x)
#range =  1500:nrow(traindata$x)
#range= 2:30 FOR DEBUG
for(nfeatures in range){
  cat("\nRandom  ")
  cat(nfeatures)
  cat("\n\n")
  model = NULL
  localtraindata <- traindata
  localtraindata$x <- localtraindata$x[random_rank[1:nfeatures],]
  if (TUNE){
    #treinamento
    tuned <- pamr.train(localtraindata, scale.sd=TRUE)
    scales <- pamr.adaptthresh(tuned)
    model <- pamr.train(localtraindata, threshold.scale=scales,n.threshold=thresholds_number, scale.sd=TRUE)
  }else{
    model <- pamr.train(localtraindata, n.threshold=thresholds_number, scale.sd=TRUE)
  }
  
  # classifica o conjunto de teste baseando-se no modelo criado  
  pred = pamr.predict(model, testdata$x[random_rank[1:nfeatures],], type = "class", threshold=0.000000000000000000000000000000000000000000000000000000000000)
  
  # se a predição corresponde à predição do N escolhido pelo Double, armazene para executar o caret
  if (nfeatures == N){
    predN <- pred
  }
  # verifica quantas amostras foram classificadas corretamente
  accuracy <- (length(which(as.logical(pred == testdata$y)))/length(testdata$y))
  random_rank_accuracy<-rbind(random_rank_accuracy, accuracy)
}

n_range= 2:nrow(traindata$x)


# salva a matriz contendo os valores de N=1..max e suas respectivas taxas de acertos de todos os rankings
results = cbind(n_range,rank_accuracy,unlike_rank_accuracy,random_rank_accuracy)
colnames(results) <- c("N","Rank","knaR","random rank")
write.matrix(results, file = "./results/nsc/independent/rankings_scores_nsc_independent_test.csv", sep = ",")

#n_range=2:30 FOR DEBUG
pdf("./results/nsc/independent/ranking_nsc_indepentend_test_validation_N_values.pdf")
plot(n_range, rank_accuracy)
dev.off()

pdf("./results/nsc/independent/gniknar_nsc_indepentend_test_validation_N_values.pdf")
plot(n_range, unlike_rank_accuracy)
dev.off()

pdf("./results/nsc/independent/random_ranking_nsc_indepentend_test_validation_N_values.pdf")
plot(n_range, random_rank_accuracy)
dev.off()

sink("./results/nsc/independent/caret_nsc_independent_test.txt")
confusionMatrix(predN, ytest, positive=NULL)
sink()



