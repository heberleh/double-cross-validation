library(rpart)
library(MASS)
library(class)
library(e1071)
library(dismo)
library(caret)


# amostras de cada classe geradas aleatoriamente:
# 1 a 6 -> 3
# 7 a 12 -> 11
# 13 a 18 -> 13

# N selecionado como N do double cross
N = as.integer(read.csv("./double_selected_N_svmrfe.txt", header=FALSE)[1,1])

# Esclha entre "tunar" ou não os parâmetros do SVM para cada valor de N, para cada ranking.
TUNE = FALSE
tuned = NULL # iniciando variável de tune


# leitura de treino e teste
db <- read.table("./spectral_counts_no_zeros_input.txt", header=TRUE,sep="\t")

# transposta, colunas serão proteínas
db2 <-t(db)

# descarta header
matrix <- as.matrix(db2[3:nrow(db2),2:(ncol(db2))])
class(matrix) <- 'numeric'


#nomes das colunas (proteínas) e classes das amostras
colnames(matrix)<-db2[2,2:ncol(db2)]
y <- as.factor(db2[3:nrow(db2),1])

test_index = c(3,11,13)

aux <- matrix[-test_index,]
nonzeros <-  (colSums(abs(aux))>0.000)

# save constant proteins from the train set
#zeros_index = which(as.logical(nonzeros == FALSE))+1  #+1 because of 1 line of header
#zeros_names = db[zeros_index,2]

# remove proteins wich were constant-zero in the used train set of double-cross
matrix <- matrix[,nonzeros]

#scale *features* with mean zero and standard deviation one
#for(i in 1:ncol(matrix))matrix[,i] = (matrix[,i]-mean(matrix[,i]))/sd(matrix[,i])

# define train and test sets
x <- matrix[-test_index,]
xtest <- matrix[test_index,]
ytest <- y[test_index]
y <- y[-test_index]


# definindo os valores de classes
#r = 5
#y <- as.factor(c(rep('normal',r),rep('carcinoma',r),rep('melanoma',r)))
#r=1
#ytest <- as.factor(c(rep('normal',r),rep('carcinoma',r),rep('melanoma',r)))


# indexes of ranking list from double cross validation
rankdb <- read.csv("./rank_index_svmrfe.csv", header=FALSE)
ranking_index = as.matrix(rankdb)

predN = c()
cat("Ranking\n")
rank_accuracy <- {}
range = 2:ncol(x)
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
  # verifica quantas amostras foram classificadas corretamente
  accuracy <- (length(which(as.logical(pred == ytest)))/length(ytest))
  rank_accuracy<-rbind(rank_accuracy,accuracy)
}
cat("\n\n")

cat("gniknaR\n")
unlike_rank_accuracy <- {}
max = ncol(x)
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
  # verifica quantas amostras foram classificadas corretamente
  accuracy <- (length(which(as.logical(pred == ytest)))/length(ytest))
  unlike_rank_accuracy<-rbind(unlike_rank_accuracy,accuracy)   
}

cat("\n\n")
cat("Random\n")
## generate a random ordering
set.seed(1) ## make reproducible here, but not if generating many random samples
random_rank <- sample(ncol(x))
random_rank_accuracy <- {}
range = 2:ncol(x)
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
  # verifica quantas amostras foram classificadas corretamente
  accuracy <- (length(which(as.logical(pred == ytest)))/length(ytest))
  random_rank_accuracy<-rbind(random_rank_accuracy,accuracy)  
}

n_range= 2:ncol(x)
pdf("plot_ranking_svmrfe.pdf")
plot(n_range, rank_accuracy)
dev.off()

pdf("plot_gniknar_svmrfe.pdf")
plot(n_range, unlike_rank_accuracy)
dev.off()

pdf("plot_random_ranking_svmrfe.pdf")
plot(n_range, random_rank_accuracy)
dev.off()

# salva a matriz contendo os valores de N=1..max e suas respectivas taxas de acertos de todos os rankings
results = cbind(n_range,rank_accuracy,unlike_rank_accuracy,random_rank_accuracy)
colnames(results) <- c("N","Rank","knaR","random rank")
write.matrix(results, file = "rankings_scores_svmrfe.csv", sep = ",")


sink("caret_svmrfe_independent_test.txt")
confusionMatrix(predN, ytest, positive=NULL)
sink()



