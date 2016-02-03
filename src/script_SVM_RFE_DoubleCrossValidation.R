
#Installing required packages
list.of.packages <- c("MASS", "class", "e1071", "dismo", "caret", "rpart", "parallel")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


library("rpart")
library("MASS")
library("class")
library("e1071")
library("dismo")
library("caret")
library("parallel")
require(foreach)
require(doSNOW)




# LOOK FOR THE SPECIFIC "TEXT" TO EDIT INPUT OR PARAMETERS:
# #=========================== SETUP YOUR INPUT AND PARAMETERS ==============================
# SERACH FOR THE ABOVE STRING







# FUNCTIONS, DO NOT EDIT, UNLESS YOU WANT TO CHANGE THE SCRIPT
data = vector()

################################################
# Feature Ranking with SVM-RFE
################################################
svmrfeFeatureRanking = function(x,y){
    n = ncol(x)
    
    survivingFeaturesIndexes = seq(1:n)
    featureRankedList = vector(length=n)
    rankedFeatureIndex = n
    
    while(length(survivingFeaturesIndexes)>0){
        #train the support vector machine
        svmModel = svm(x[, survivingFeaturesIndexes], y, cost = 10, cachesize=500,  scale=F, type="C-classification", kernel="linear" )
        
        #compute the weight vector
        w = t(svmModel$coefs)%*%svmModel$SV
        
        #compute ranking criteria
        rankingCriteria = w * w
        
        #rank the features
        ranking = sort(rankingCriteria, index.return = TRUE)$ix
        
        #update feature ranked list
        featureRankedList[rankedFeatureIndex] = survivingFeaturesIndexes[ranking[1]]
        rankedFeatureIndex = rankedFeatureIndex - 1
        
        #eliminate the feature with smallest ranking criterion
        (survivingFeaturesIndexes = survivingFeaturesIndexes[-ranking[1]])

    }
    
    return (featureRankedList)
}


################################################
# Feature Ranking with Average Multiclass SVM-RFE
################################################

auxm <- list()
svmrfeFeatureRankingForMulticlass = function(x,y,tuned){
    n = ncol(x)
    
    survivingFeaturesIndexes = seq(1:n)
    featureRankedList = vector(length=n)
    featureRankedListWeights = list()
    rankedFeatureIndex = n
    rankingCriteriaFull = vector()
    while(length(survivingFeaturesIndexes)>0){
        #train the support vector machine
        svmModel = svm(x[, survivingFeaturesIndexes], y, tuned$best.parameters[2], gamma = tuned$best.parameters[1], scale=T, type="C-classification", kernel="linear" )
        
        #compute the weight vector
        multiclassWeights = svm.weights(svmModel)
		
		    #calcular o erro baseado nos pesos médios
		    #N=bestN()
        
        #compute ranking criteria
        multiclassWeights = multiclassWeights * multiclassWeights
        rankingCriteria = 0
        
        
        for(i in 1:ncol(multiclassWeights))rankingCriteria[i] = mean(multiclassWeights[,i])


#         if (rankedFeatureIndex == n){
#           #write.matrix(multiclassWeights, file = "multiclassWeights.csv", sep = ",") 
#           auxm<-multiclassWeights
#           rankingCriteriaFull = rankingCriteria
#         }
        
        #rank the features
        (ranking = sort(rankingCriteria, index.return = TRUE)$ix)      

        #gene weight - which interaction
        featureRankedListWeights[[survivingFeaturesIndexes[ranking[1]]]] = multiclassWeights[,ranking[1]]
        rankingCriteriaFull[survivingFeaturesIndexes[ranking[1]]] = mean(multiclassWeights[,ranking[1]])
        
		#ranking[i] é o index do index do atributo. em surviving pode ter restado o atributo de index 200 na posição 1. ou seja, surv[ranking[1]] = 200 -> atributo de index 200 da tabela original.
		
        #update feature ranked list
        (featureRankedList[rankedFeatureIndex] = survivingFeaturesIndexes[ranking[1]])
        rankedFeatureIndex = rankedFeatureIndex - 1        
        
        #eliminate the feature with smallest ranking criterion
        (survivingFeaturesIndexes = survivingFeaturesIndexes[-ranking[1]])
        #cat(length(survivingFeaturesIndexes),"\n")
    }
    
    result = list(featureRankedList=featureRankedList, mean_weights = rankingCriteriaFull, weights=featureRankedListWeights)
	return (result)
}

################################################
# This function gives the weights of the hiperplane
################################################
svm.weights<-function(model){
w=0
  if(model$nclasses==2){
       w=t(model$coefs)%*%model$SV
  }else{    #when we deal with OVO svm classification
      ## compute start-index
      start <- c(1, cumsum(model$nSV)+1)
      start <- start[-length(start)]

      calcw <- function (i,j) {
        ## ranges for class i and j:
        ri <- start[i] : (start[i] + model$nSV[i] - 1)
        rj <- start[j] : (start[j] + model$nSV[j] - 1)

      ## coefs for (i,j):
        coef1 <- model$coefs[ri, j-1]
        coef2 <- model$coefs[rj, i]
        ## return w values:
        w=t(coef1)%*%model$SV[ri,]+t(coef2)%*%model$SV[rj,]
        return(w)
      }

      W=NULL
      for (i in 1 : (model$nclasses - 1)){
        for (j in (i + 1) : model$nclasses){
          wi=calcw(i,j)
          W=rbind(W,wi)
        }
      }
      w=W
  }
  return(w)
}

balanced.folds <- function (y, nfolds = min(min(table(y)), 10)) 
{
  totals <- table(y)
  fmax <- max(totals)
  nfolds <- min(nfolds, fmax)
  nfolds = max(nfolds, 2)
  folds <- as.list(seq(nfolds))
  yids <- split(seq(y), y)
  bigmat <- matrix(NA, ceiling(fmax/nfolds) * nfolds, length(totals))
  for (i in seq(totals)) {
    cat(i)
    if (length(yids[[i]]) > 1) {
      bigmat[seq(totals[i]), i] <- sample(yids[[i]])
    }
    if (length(yids[[i]]) == 1) {
      bigmat[seq(totals[i]), i] <- yids[[i]]
    }
  }
  smallmat <- matrix(bigmat, nrow = nfolds)
  smallmat <- permute.rows(t(smallmat))
  res <- vector("list", nfolds)
  for (j in 1:nfolds) {
    jj <- !is.na(smallmat[, j])
    res[[j]] <- smallmat[jj, j]
  }
  return(res)
}

permute.rows <- function (x)
{
  dd <- dim(x)
  n <- dd[1]
  p <- dd[2]
  mm <- runif(length(x)) + rep(seq(n) * 10, rep(p, n))
  matrix(t(x)[order(mm)], n, p, byrow = TRUE)
}

getLowerNLowerError <- function(errors, Ns){
  min_error = min(errors)
  index_min_error = which(errors == min_error, arr.ind=TRUE)
  min_N = min(Ns[index_min_error])
  result = list(error = min_error, N = min_N)
  return (result)
}






#=========================== SETUP YOUR INPUT AND PARAMETERS ==============================

# PARAMETERS:    file.txt
# SEE THE train.txt AND FOLLOW THE PATTERN
db <- read.table("./dataset/current/train.txt", header=TRUE,sep="\t")


global_tune = FALSE        # DECIDE IF THE SCRIPT SHOULD TUNE THE MODEL PARAMETERS IN EVERY LOOP OF EACH REPETITION
                           # OR USE DEFAULT COST AND GAMMA VALUES: 10 AND 10^-4
                           # NOTE THAT WITH SMALL NUMBER OF SAMPLES THE TUNE FUNCTION CAN'T RUN 
                           # THE K-FOLD CROSS VALIDATION IF THE FOLDS CONTAIN ONLY ONE OR 2 SAMPLES, OR WORST: ZERO SAMPLES.                           
                           

repetition <- 2           # NUMBER OF DOUBLE CROSS VALIDATION THAT WILL BE EXECUTED
                           # default: 100
                                                      

#nCluster = detectCores()                           
nCluster = 2               # NUMBER OF CORES/NUCLEUS FOR PARALLEL PROCESSING. 
                           # FOR INSTANCE, Intel i7 processors might have 8 virtual nucleus.
                           # IF YOU SET nCluster = 8, THE SCRIPT WILL RUN
                           # 8 DOUBLE CROSS VALIDATION AT THE SAME TIME, ONE PER CORE.
                           # I PREFERE TO SET nCluster = <NUMBER OF REAL CORES>
                           # In general, 8 virtual nucleus are actually 4 real cores.
                           # You can use the command "nCluster = detectCores()" to
                           # detect the max. number of cores (including virtual ones).
                           # In 8 virtual cores processor, the function detectCores() will
                           # return: 8.






#diminuir a dimensionalidade apenas para testes rápidos do código
#db2 <- t(db[1:500,])
db2 <-t(db)

secretome.x <- as.matrix(db2[3:nrow(db2),2:(ncol(db2))])
class(secretome.x) <- "numeric"
#colnames(x)<-1:ncol(x)
colnames(secretome.x)<-db2[2,2:ncol(db2)]
#yvector<-as.vector(unlist(db2[3:ncol(db2),1]))
secretome.y <- as.factor(db2[3:nrow(db2),1])


# z-score? Models like SVM have padronization by z-score by default.
# Do not use this z-score code, unless you really want and know what
# you are doing.
# part of the code padronizes by lines.
# the other part padronizes by columns.
#x = log(x)
#scale samples with mean zero and standard deviation one
#for(i in 1:nrow(x))x[i,] = (x[i,]-mean(x[i,]))/sd(x[i,])
#scale features with mean zero and standard deviation one
#for(i in 1:ncol(x))x[,i] = (x[,i]-mean(x[,i]))/sd(x[,i])
#x = 2*atan(x/2)

#featureRankedList = svmrfeFeatureRankingForMulticlass(secretome.x, secretome.y)


selected_outer = list()
final_outer_error = vector()
final_outer_error_min = vector()
final_outer_error_max = vector()
final_N = vector()
logdouble = list()

cl = makeCluster(nCluster,type="SOCK")
registerDoSNOW(cl) 
getDoParWorkers() 
allrep = list()

stime <- system.time({
allrep <- foreach (rep = 1:repetition, .combine="cbind", .packages = c("e1071","dismo","class","rpart","MASS")) %dopar%{  
#for (rep in 1:repetition){
  
  folds <- balanced.folds(y=secretome.y)
  nfold <- length(folds)

  selected_inner = list()
  outer_error = vector()
  mean_error_by_N =  list()
  Nf = vector()
  logouter = list()
  outer_N = vector()
  
  for (i in 1:nfold){
    #definindo conjunto de treinamento i
    secretome2.x <- secretome.x[-folds[[i]], ,drop=FALSE]
    secretome2.y <- secretome.y[-folds[[i]]]
    secretome2.testX <- secretome.x[folds[[i]],,drop=FALSE]
    secretome2.testY <- secretome.y[folds[[i]]]
    
    #inner loop - k-fold-cross-validation
    folds_inner <- balanced.folds(y=secretome2.y)
    nfold_inner <- length(folds_inner)
     
    inner_rank = list()
    inners_errors = list()
    inner_rankCriteria = list()
    for (ii in 1:nfold_inner){
      inner_error = vector() 
      
      #definindo conjunto de treinamento ii      
      secretome2_inner.x <- secretome2.x[-folds_inner[[ii]], ,drop=FALSE]
      secretome2_inner.y <- secretome2.y[-folds_inner[[ii]]]
      secretome2_inner.testX <- secretome2.x[folds_inner[[ii]],,drop=FALSE]
      secretome2_inner.testY <- secretome2.y[folds_inner[[ii]]]
      
      
      #treino inner
      tuned_inner = NULL
      if (global_tune == FALSE){        
        tuned_inner$best.parameters[1] = 10^(-4)
        tuned_inner$best.parameters[2] = 10        
      }else{      
        cross <- length(balanced.folds(y=secretome2_inner.y))
        tuned_inner <- tune.svm(x = secretome2_inner.x, y=secretome2_inner.y, gamma = 10^(-6:-1), cost = 10^(1:2), tunecontrol = tune.control(cross=cross))
      }
      
      aux = svmrfeFeatureRankingForMulticlass(secretome2_inner.x, secretome2_inner.y, tuned_inner)
      inner_rank = aux$featureRankedList
      inner_rankCriteria[[ii]] = aux$mean_weights
      
  
      #test inner
      iierror = rep(1,length(inner_rank))
      for(pow in 1:length(inner_rank)){
        nfeatures = pow                
        svmModel = svm(secretome2_inner.x[, inner_rank[1:nfeatures]], secretome2_inner.y, cost = tuned_inner$best.parameters[2], gamma = tuned_inner$best.parameters[1], scale = T, type="C-classification", kernel="linear")
        pred <- predict(svmModel, secretome2_inner.testX[, inner_rank[1:nfeatures]])
        inner_error[pow] = length(which(!as.logical(pred == secretome2_inner.testY)))/length(secretome2_inner.testY)
      }           
           
      #armazena os erros das voltas dos inners
      inners_errors[[ii]] <- inner_error      
    }#fim inner
    
    siz=length(inners_errors[[1]])
    
    rankingCriteriaInners = rep(1,siz)
    for (fo in 1:nfold_inner){
      rankingCriteriaInners = rankingCriteriaInners * inner_rankCriteria[[fo]]      
    }    
    
    #rank the features
    (ranking = rev(sort(rankingCriteriaInners, index.return = TRUE)$ix))
    
    mean_error_by_N[[i]] = vector()
    aux_matrix = list()
    for(N in 1:siz){
      aux_matrix[[N]] = vector()
    }      
    for(ii in 1:nfold_inner){
      for(N in 1:siz){
        aux_matrix[[N]][ii] = inners_errors[[ii]][N]
      }      
    }
    for(N in 1:siz){
      mean_error_by_N[[i]][N] = mean(aux_matrix[[N]])
    }  
    
    
    current_mean_error = mean_error_by_N[[i]]
    minerror = min(current_mean_error)
    indexN = 1
    for (N in 1:length(current_mean_error)){
      if(current_mean_error[N] == minerror){
        indexN = N
        break
      }
    }
    Nf[i] = indexN
    nfeatures = indexN
    
    #treino outer
    #aux = svmrfeFeatureRankingForMulticlass(secretome2.x, secretome2.y)
    #outer_rank = aux$featureRankedList
    outer_rank = ranking
    
    tuned_outer = NULL
    svmModel = NULL
    if (global_tune==TRUE){
    #seleciona apenas as N (escolhido no inner)
    cross_outer <- length(balanced.folds( y=secretome2.y))
    tuned_outer <- tune.svm(x = secretome2.x[, outer_rank[1:nfeatures]], y=secretome2.y, gamma = 10^(-6:-1), cost = 10^(1:2),tunecontrol = tune.control(cross=cross_outer))    
    }else{
      tuned_outer$best.parameters[1] = 10^(-4)
      tuned_outer$best.parameters[2] = 10
    }
    svmModel = svm(x = secretome2.x[, outer_rank[1:nfeatures]], y = secretome2.y, cost = tuned_outer$best.parameters[2], gamma = tuned_outer$best.parameters[1], scale = T, type="C-classification", kernel="linear")   
    
    #testa      
    pred <- predict(svmModel, secretome2.testX[, outer_rank[1:nfeatures]])
    outer_error[i] = length(which(!as.logical(pred == secretome2.testY)))/length(secretome2.testY)
    outer_N[i] = nfeatures
    logouter[[i]] = list(pred=pred, ref=secretome2.testY, outer_error = outer_error[i], fold = folds[[i]])
    
  }#fim outer
  
  logdouble[[rep]] = logouter
  
  siz = length(mean_error_by_N[[1]])
  mean_mean_inner_error = vector()
  aux_matrix = list()
  for(N in 1:siz){
    aux_matrix[[N]] = vector()
  } 
  for(fo in 1:nfold_inner){
    for(N in 1:siz){
      aux_matrix[[N]][fo] = mean_error_by_N[[fo]][N]
    }      
  }
  for(N in 1:siz){
    mean_mean_inner_error[[N]] = mean(aux_matrix[[N]])
  }  
  
  minerror = min(mean_mean_inner_error)
  for (N in 1:siz){
    if(minerror == mean_mean_inner_error[N]){
      final_outer_N = N
      break
    }
  }
  
  selected_outer[[rep]] = getLowerNLowerError(outer_error,outer_N)
  finalN_testado = selected_outer[[rep]]$N
  finalError_testado =selected_outer[[rep]]$error
  
  result = list(finalN = final_outer_N, finalN_testado = finalN_testado, finalError_testado = finalError_testado, final_outer_error = mean(outer_error), final_outer_error_min = min(outer_error), final_outer_error_max = max(outer_error), logdouble=logdouble)
  result
  
#   final_N[rep] = final_outer_N
#   final_outer_error[rep] = mean(outer_error)
#   final_outer_error_min[rep] = min(outer_error)
#   final_outer_error_max[rep] = max(outer_error)
#   
  #selected_outer[[rep]] = getLowerNLowerError(outer_error,outer_N)

  #sink()  
}#fim repeti??es
})[3]
stopCluster(cl)
stime
  
double_error = vector()
all_N = vector()
for (rep in 1:repetition){  
  double_error[rep] = allrep[[4,rep]]
  all_N[rep] = allrep[1,rep]  
}

double_error = as.vector(unlist(double_error))
all_N = as.vector(unlist(all_N))
mean_double_error = mean(double_error)
min_double_error = min(double_error)
max_double_error = max(double_error)
dif_double_error = abs(double_error - mean_double_error)
min_dif = min(dif_double_error)
index_error = which(dif_double_error == min_dif)
min_N = min(all_N[index_error])
index_N = which(all_N == min_N)
index = intersect(index_N,index_error)
final_double_error = mean(double_error[index])

if(length(index) != 1){
  cat ("Mais de um modelo pode ser utilizado para calcular as métricas!")
  index = index[1]
}

index2 = (7*(index))
size = length(allrep[[index2]][[index]])
pred = list()
ref = list()
for (i in 1:size){
  obj = allrep[[index2]][[index]][[i]]
  pred[[i]] = obj$pred
  ref[[i]] = obj$ref
}
pred = as.vector(unlist(pred))
ref = as.vector(unlist(ref))


sink("./results/svm-rfe/caret_svm-rfe.txt")
cat("Average error of  double cross-validation repetitions: ")
cat(mean_double_error)
cat("\n\n")
cat("Maximum error between repetitions: ")
cat(max_double_error)
cat("\n\n")
cat("Minimum error between repetitions: ")
cat(min_double_error)
cat("\n\n")
cat("Best model indexes: ")
cat(index_error)
cat("\n\n")
cat("All possible N according to selected models:")
cat(all_N[index_error])
cat("\n\n")
cat("Minimum N (selected final N): ")
cat(min_N)
cat("\n\n")
cat("Double-Cross-Validation Error: ")    
cat(final_double_error)
cat("\n\n")

pred=factor(pred, levels=unique(secretome.y))
ref=factor(ref, levels=unique(secretome.y))
     
confusionMatrix(pred, ref, positive=NULL)
sink()

m = cbind(all_N, double_error)
write.csv(m,"./results/smv-rfe/all_n_and_double_error_repetition_svm-rfe.csv")

tuned = NULL
#Tuna para encontrar parâmetros com o conjunto completo
if (global_tune==TRUE){
  tuned <- tune.svm(x = secretome.x, y=secretome.y, gamma = 10^(-6:-1), cost = 10^(1:2))
}else{
  tuned$best.parameters[1] = 10^(-4)
  tuned$best.parameters[2] = 10
}

#Calcula SVM-RFE
aux = svmrfeFeatureRankingForMulticlass(secretome.x, secretome.y, tuned)
featureRankedList = aux$featureRankedList
mean_weights = aux$mean_weights[featureRankedList]
weights = aux$weights
weights = do.call(rbind, weights)
weights = weights[featureRankedList,]



# #6-fold
# nf <- {}
# Accuracy <- {}
# for(pow in 1550:1600){
#   nfeatures = pow
#   truePredictions = 0
#   svmModel = svm(x[, featureRankedList[1:nfeatures]], y, cost = 10, gamma=0.0001, cachesize=500,  scale=T, type="C-classification", kernel="linear",cross=5) 
#   nf<-rbind(nf,nfeatures)
#   Accuracy<-rbind(Accuracy,mean(svmModel$accuracies))  
# }
# plot(nf,Accuracy)

x = secretome.x
names = colnames(subset(x,select=featureRankedList)) # ou... colnames(x[,featureRankedList])
x_ordered = cbind(names,t(x[,featureRankedList]))
merged = cbind(names,featureRankedList,1:length(names),weights) #merge(names,featureRankedList, by = "row.names", all = TRUE)
#merged = as.matrix(merged[-1])

#colnames para 2 classes (dataset4):
#colnames(merged) <- c("gene name","index","rank","v1")

#colnames para 2 classes (secretoma):
#colnames(merged) <- c("gene name","index","rank","v1","v2")

#colnames para 3 classes (secretoma):
#colnames(merged) <- c("gene name","index","rank","v1","v2","v3")

#colnames para 4 classes:
colnames(merged) <- c("gene name","index","rank","v1","v2","v3","v4","v5","v6")

#colnames para 5 classes:
#colnames(merged) <- c("gene name","index","rank","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")

#colnames para 7 classes:
#colnames(merged) <- c("gene name","index","rank","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10","v11","v12","v13","v14","v15","v16","v17","v18","v19","v20","v21")

merged = cbind(merged,x_ordered)

#write.matrix(merged, file = "big_matrix_svm-rfe.csv", sep = ",")
write.csv(merged, file = "./results/smv-rfe/big_matrix_svm-rfe.csv")

#write.matrix(x_ordered, file = "matrix_by_rank_scale.csv", sep = ",")

#v = c(55,54,25,1578,108,24,1575,51,1528,26,23,1550,34,1521,23,27,13,18,1560,1554,1520,26,25,37,36,71,43,25,29,34,52,1547,29,42,1568,164)
#d = density(v)

#plot(d)
#axis(1,at=36,labels=c("36"))




