
#Installing required packages
list.of.packages <- c("MASS", "pamr", "parallel")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library("MASS")
library("pamr")
library("parallel")

# LOOK FOR THE SPECIFIC "TEXT" TO EDIT INPUT OR PARAMETERS:
# #=========================== SETUP YOUR INPUT AND PARAMETERS ==============================
# SERACH FOR THE ABOVE STRING



# FUNCTIONS, DO NOT EDIT, UNLESS YOU WANT TO CHANGE THE SCRIPT
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

leave_one_out <- function(y){
  folds=list(list())
  for (i in 1:length(y)){
    folds[[i]]=i;
  }  
  return (folds)
}

getMeanProbabilities <- function(probabilities){
  total = vector()
  for (i in 1:nrow(probabilities)){
    total[i] = max(probabilities[i,])    
  }
  return(mean(total))  
}

getHigherThresholdWithMaxProbabilitySum <- function(thresholds,probabilities){
  total_probability = vector()
  for (t in 1:(length(thresholds))){
    total_probability[t] = 0.0
    p <- probabilities[,,t,drop=FALSE]
    for (s in 1:(nrow(p))){  
      if (ncol(p)>1){
        total_probability[t] = sum(total_probability[t], max(p[,s,]))
      }else{
        total_probability[t] = sum(total_probability[t], p[s])
      }
    }
    total_probability[t] = total_probability[t]/(nrow(p))
  }
  
  max_probability = max(total_probability)
  higher_threshold = min(thresholds)
  for (t in 1:(length(thresholds))){
    #if (max_probability - total_probability[t] <= 0.00001 && thresholds[t] > higher_threshold){    
    if (max_probability == total_probability[t] && thresholds[t] > higher_threshold){ 
      higher_threshold = thresholds[t]
    }
  }
  
  result = list(higher_threshold = higher_threshold, max_probability = max_probability)
  return (result)
}



getHigherThresholdWithMaxProbability <- function(thresholds,total_probability){
  max_probability = max(total_probability)  
  higher_threshold = min(thresholds)
  index = match(c(higher_threshold),thresholds)
  for (t in 1:(length(thresholds))){
    if (max_probability == total_probability[t]  && thresholds[t] >= higher_threshold){
      higher_threshold = thresholds[t]
      index = t
    }
  }
  result = list(higher_threshold = higher_threshold, probability = total_probability[index])
  return(result)
}

getLowerErrorMaxThreshold <- function(errors,thresholds){
  min_error = min(errors)
  higher_threshold = min(thresholds)
  index = match(c(higher_threshold),thresholds)
  
  for (t in 1:(length(thresholds))){
    if (min_error == errors[t]  && thresholds[t] >= higher_threshold){
      higher_threshold = thresholds[t]
      index = t
    }
  }
  result = list(higher_threshold = higher_threshold, error = errors[index])
  return(result)
}


getLowerErrorMaxProbabilityMaxThreshold <- function(errors, thresholds, probabilities){
  min_error = min(errors)
  higher_threshold = min(thresholds)
  index = match(c(higher_threshold),thresholds)
  
  index_error = which(errors == min_error, arr.ind=TRUE)
  #entre os de erro m??nimo, qual a máxima probabilidade
  max_probability = max(probabilities[index_error])
  
  index_prob = which(probabilities == max_probability, arr.ind=TRUE)
  
  index_error_prob = intersect(index_prob,index_error)
  
  selected_thresholds = thresholds[index_error_prob]
  
  max_threshold = max(selected_thresholds)
  
  index_threshold = which(thresholds == max_threshold, arr.ind=TRUE)
  
  index = intersect(index_error_prob,index_threshold)[1]
  
  result = list(higher_threshold = max_threshold, error = min_error, probability = max_probability, index=index)
  return(result)
}


make_plot <- function(path,x,y,str){
  pdf(path)
  plot(x, y, type="n", main=str)
  lines(x, y, type="l", col="gray")
  #lines(x,y,type="b", lwd=1.5)
  points(x, y, pch=4, col="blue")
  dev.off()
}




#=========================== SETUP YOUR INPUT AND PARAMETERS ==============================

# PARAMETERS:    file.txt, number_of_columns
# SEE THE train.txt AND FOLLOW THE PATTERN

aux_db <- read.table("./dataset/current/train.txt", header=TRUE,sep="\t")
initial_number_of_columns = ncol(aux_db)
secretome.data <- pamr.from.excel("./dataset/current/train.txt", initial_number_of_columns, sample.labels=TRUE)

repetition = 100         # NUMBER OF DOUBLE CROSS VALIDATION THAT WILL BE EXECUTED
                         # default: 100
                         
thresholds_number = 90   # NUMBER OF CALCULATED THRESHOLDS  default: 90
                         # default: 90

# z-score? Models like NSC have padronization by z-score by default.
# Do not use this z-score code, unless you really want and know what
# you are doing.
# part of the code padronizes by lines.
# the other part padronizes by columns.
#x = secretome.data$x
#for(i in 1:ncol(x))x[,i] = (x[,i]-mean(x[,i]))
#for(i in 1:nrow(x))x[i,] = (x[i,]-mean(x[i,]))
#for(i in 1:ncol(x))x[,i] = (x[,i]-mean(x[,i]))/sd(x[,i])
#for(i in 1:nrow(x))x[i,] = (x[i,]-mean(x[i,]))/sd(x[i,])
#secretome.data$x = x


#bestThresholdOuterLoop = vector()
#errorFold = vector()
bestFinal = list()
probs=vector()
thres=vector()
bestFinalNinner=vector()
bestFinalN=vector()
double_cross_error = vector()
folds_list = list()
probs_final = vector()
inner_inner_rep = list()
outer_rep = list()
error_outers_rep  = list()
probabilities_outer_rep = list()
selected_inner_model_final = list()
selected_outer_model_final = list()
preds_outer_rep = list()
refs_outer_rep = list()
predictmany_rep = list()



secretome.train <- pamr.train(secretome.data)
secretome.scales <- pamr.adaptthresh(secretome.train)
secretome.train2 <- pamr.train(secretome.data, threshold.scale=secretome.scales,n.threshold=thresholds_number, scale.sd=TRUE)
fixed_thresholds = secretome.train2$threshold

stime <- system.time({
for (rep in 1:repetition){
  #sink("outputgarbage",append = FALSE)
  #outer_n = min(table(secretome.data$y))
  
  #backups
  #fold = list()
  #bestThreshold = vector()
  #outer-loop k-fold-cross-validation
  
  #true_prediction = list()
  #true_prediction_total = vector()
  #i = 1
  #ii = 1
  
  #outer loop
  #folds <- balanced.folds(y=secretome.data$y,nfolds=10)
  #folds <- leave_one_out(secretome.data$y)
  folds <- balanced.folds(y=secretome.data$y)
  
  folds_list[[rep]]<-folds
  
  nfold <- length(folds)  
  secretome2.data<- secretome.data  
  cat(paste("\nK for Outer Loop:",nfold,"\n",sep=" "))
  mean_outer_total_probability = vector()  
  selected_inner_values = list()
  selected_inner_model = list()
  inner_train =list()
  inner_train_data = list()
  bestN = vector()
  error_outers =vector()
  pred_outers = list()
  ref_outers = list()
  probabilities_outer = vector()  
  thresholds_outers = vector()
  inner_inner = list()
  outer_model = list()
  predictmany = list()
  ref=list()
  for (i in 1:nfold){
    #definindo conjunto de treinamento i
    secretome2.data$x <- secretome.data$x[, -folds[[i]],drop=FALSE]
    secretome2.data$y <- secretome.data$y[-folds[[i]]]
    secretome2.data$genenames <- secretome.data$genenames[-folds[[i]]]
    secretome2.data$geneid <- secretome.data$geneid[-folds[[i]]]
    secretome2.data$samplelabels <- secretome.data$samplelabels[-folds[[i]]]
    secretome2.testX <- secretome.data$x[,folds[[i]],drop=FALSE]
    secretome2.testY <- secretome.data$y[folds[[i]]]
      
        #inner loop - k-fold-cross-validation
        #folds_inner <- leave_one_out(secretome2.data$y)
        folds_inner <- balanced.folds(y=secretome2.data$y)
        #folds_inner <- balanced.folds(y=secretome2.data$y,nfolds=2)
        nfold_inner <- length(folds_inner)
        cat(paste("\nK for Inner Loop:",nfold_inner,"\n",sep=" "))
        secretome2_inner.data <- secretome2.data
        inner = list()
        max_total_probability = vector()
        higher_inner_threshold = vector()
        errors_inner = vector()
        inner_model = list()
        for (ii in 1:nfold_inner){
          #definindo conjunto de treinamento ii      
          secretome2_inner.data$x <- secretome2.data$x[, -folds_inner[[ii]],drop=FALSE]
          secretome2_inner.data$y <- secretome2.data$y[-folds_inner[[ii]]]
          secretome2_inner.data$genenames <- secretome2.data$genenames[-folds_inner[[ii]]]
          secretome2_inner.data$geneid <- secretome2.data$geneid[-folds_inner[[ii]]]
          secretome2_inner.data$samplelabels <- secretome2.data$samplelabels[-folds_inner[[ii]]]
          secretome2_inner.testX <- secretome2.data$x[,folds_inner[[ii]],drop=FALSE]
          secretome2_inner.testY <- secretome2.data$y[folds_inner[[ii]]]
          
          #treinamento
          secretome2_inner.train <- pamr.train(secretome2_inner.data)
          secretome2_inner.scales <- pamr.adaptthresh(secretome2_inner.train)
          secretome2_inner.train2 <- pamr.train(secretome2_inner.data, threshold.scale=secretome2_inner.scales,n.threshold=thresholds_number)
          
          #teste
          secretome2_inner.prob = pamr.predictmany(secretome2_inner.train2, secretome2_inner.testX)                
          
          preds_inner = secretome2_inner.prob$predclass
          thresholds = secretome2_inner.train2$threshold
          probabilities_inner = vector()     
          nex=0
          for (iie in 1:ncol(preds_inner)){
            errors_inner[iie] = length(which(!as.logical(preds_inner[,iie] == secretome2_inner.testY)))/length(secretome2_inner.testY)
            prob_matrix = secretome2_inner.prob$prob[,,iie,drop=FALSE]
            total_probability = 0
            if(ncol(prob_matrix)==1){       
                nex=1
                total_probability = total_probability + max(prob_matrix[,1,,drop=FALSE])              
            }else{
              nex = nrow(prob_matrix)
              for (iiie in 1:nrow(prob_matrix)){
                total_probability = total_probability + max(prob_matrix[iiie,,,drop=FALSE])
              }              
            }
            probabilities_inner[iie] = total_probability/nex
          }
                    
          inner[[ii]] <- getLowerErrorMaxProbabilityMaxThreshold(errors_inner,thresholds,probabilities_inner)
          inner_train[[ii]] <- secretome2_inner.train2
          inner_train_data[[ii]] <- secretome2_inner.data
          inner_model[[ii]] <- list(data=secretome2_inner.data,train=secretome2_inner.train2,testX=secretome2_inner.testX,testY=secretome2_inner.testY)
        }
    
    inner_inner[[i]] <- inner
    errors_inners = vector()
    thresholds_inners = vector()
    probabilities_inners = vector()
    for (ii in 1:nfold_inner){
      errors_inners[ii] <- inner[[ii]]$error
      thresholds_inners[ii] <- inner[[ii]]$higher_threshold
      probabilities_inners[ii] <- inner[[ii]]$probability
    }
    
    selected_inner_values[[i]] = getLowerErrorMaxProbabilityMaxThreshold(errors_inners, thresholds_inners, probabilities_inners)    
    current_threshold = selected_inner_values[[i]]$higher_threshold
    selected_inner_model[[i]] <- inner_model[[selected_inner_values[[i]]$index]]
        
#     #seleciona features do melhor modelo do inner e re-prepara o conjunto de treinamento e de teste do outer
#     index <- match(c(current_threshold),thresholds_inners)        
#     selected_features_matrix <- pamr.listgenes(inner_train[[index]], inner_train_data[[index]],  current_threshold, genenames=TRUE)      
#     #can't use just 1 variable in model.
#     while(nrow(selected_features_matrix)<2){
#         current_threshold= current_threshold-2
#         selected_features_matrix <- pamr.listgenes(inner_train[[index]], inner_train_data[[index]],  current_threshold, genenames=TRUE)
#     }    
#     selected_features <- selected_features_matrix[,2] 
#     bestN[i] = nrow(selected_features_matrix)
#     features_index <- which(secretome2.data$genenames %in% selected_features)
#     secretome2.data$x <- secretome2.data$x[features_index,,drop=FALSE]    
#     secretome2.data$genenames <- secretome2.data$genenames[features_index]
#     secretome2.data$geneid <- secretome2.data$geneid[features_index]    
#     secretome2.testX <- secretome2.testX[features_index,,drop=FALSE]
#      
    
    #treinamento outer loop
    secretome2.train <- pamr.train(secretome2.data)
    secretome2.scales <- pamr.adaptthresh(secretome2.train)
    secretome2.train2 <- pamr.train(secretome2.data, threshold.scale=secretome2.scales,n.threshold=thresholds_number)
    
    #teste
    secretome2.class_pred = pamr.predict(secretome2.train2, secretome2.testX, current_threshold, type=c("class"))  
    secretome2.class_prob = pamr.predict(secretome2.train2, secretome2.testX, current_threshold, type=c("posterior"))
    
    pred_outers[[i]] = secretome2.class_pred
    ref_outers[[i]] = secretome2.testY
    error_outers[i] = length(which(!as.logical(secretome2.class_pred == secretome2.testY)))/length(secretome2.testY)
                
    prob_matrix = secretome2.class_prob
    total_probability = 0          
    for (iiie in 1:nrow(prob_matrix)){
        total_probability = total_probability + max(prob_matrix[iiie,,drop=FALSE])
    }    
    probabilities_outer[i] = total_probability / nrow(prob_matrix)
    
    
    outer_model[[i]] <- list(data = secretome2.data, train = secretome2.train2, testX = secretome2.testX, testY = secretome2.testY)
    
#   predictmany[[i]] = pamr.predictmany(fit=secretome2.train2, newx=secretome2.testX, threshold=fixed_thresholds)
    ref[[i]] = secretome2.testY
  }
  
  
  
  preds_outer_rep[[rep]] = pred_outers
  refs_outer_rep[[rep]] = ref_outers
  
  for (tr in 1:nfold){
    thresholds_outers[tr] = selected_inner_values[[tr]]$higher_threshold
    #errors_outers[tr] = threshold_outer[[tr]]$error não usar esse "error" pois é do inner, não do outer    
  }
    
  inner_inner_rep[[rep]] = inner_inner
  error_outers_rep[[rep]]  = error_outers
  probabilities_outer_rep[[rep]] = probabilities_outer
  
  #treinamento final
  bestFinal[[rep]] = getLowerErrorMaxProbabilityMaxThreshold(error_outers,thresholds_outers,probabilities_outer)  
  secretome.train <- pamr.train(secretome.data)
  secretome.scales <- pamr.adaptthresh(secretome.train)
  
  selected_inner_model_final[[rep]] = selected_inner_model[[bestFinal[[rep]]$index]]
  selected_outer_model_final[[rep]] = outer_model[[bestFinal[[rep]]$index]]
  
  secretome.train2 <- pamr.train(secretome.data, threshold.scale=secretome.scales)  
  

  #round(bestFinal[[rep]]$higher_threshold, digits = 3)
  if (bestFinal[[rep]]$higher_threshold > secretome.train2$threshold[length(secretome.train2$threshold)]){    
    bestFinal[[rep]]$higher_threshold <- round(secretome.train2$threshold[length(secretome.train2$threshold)], digits=3)
  }
  
  selected_features_matrix_rep <- nrow(pamr.listgenes(secretome.train2, secretome.data,  bestFinal[[rep]]$higher_threshold))
  
  #secretome.class_pred = pamr.predict(secretome.train2, secretome.data$x, bestFinal[[rep]]$higher_threshold, type=c("class"))  
  double_cross_error[rep] = mean(error_outers)   
  bestFinalN[rep] = selected_features_matrix_rep 
  thres[rep] = bestFinal[[rep]]$higher_threshold
  probs_final[rep] = bestFinal[[rep]]$probability
  #sink()
}
})[3]
stime



#Entre todas as repetiCOes do Double-Cross-Validation, selecione a cujo erro mais se aproxima da mEdia
mean_double_cross_validation = mean(double_cross_error)
min_double_cross_validation = min(double_cross_error)
max_double_cross_validation = max(double_cross_error)
dif = abs(double_cross_error - mean_double_cross_validation)
min_dif = min(dif)
index_min_dif = which(dif == min_dif)
max_prob = max(probs_final[index_min_dif])
index_max_prob = which(probs_final == max_prob)
inter_index = intersect(index_min_dif, index_max_prob)
max_thre = max(thres[inter_index])
index_thres = which(thres == max_thre)
index = index_thres
if(length(index_thres) > 1){
  cat("H? mais de uma possibilidade de threshold final!")
  index = index_thres[1]
}

secretome.train <- pamr.train(secretome.data)
secretome.scales <- pamr.adaptthresh(secretome.train)
secretome.train2 <- pamr.train(secretome.data, threshold.scale=secretome.scales,n.threshold=thresholds_number, scale.sd=TRUE)  

list_genes_final=pamr.listgenes(secretome.train2, secretome.data, thres[index])
complete_genes_list = pamr.listgenes(secretome.train2, secretome.data, 0, genenames=TRUE)
write.csv(complete_genes_list,"./results/nsc/gene_ranking_nsc_complete.csv")
complete_genes_names_list = complete_genes_list[,2]

# identifica proteínas descartadas para por no fim do ranking
constant_proteins_for_nsc <- as.matrix(secretome.data$genenames[!(secretome.data$genenames %in% complete_genes_names_list )])
constant_proteins_for_nsc_index <- match(constant_proteins_for_nsc,secretome.data$genenames)

write.matrix(cbind(secretome.data$genenames[constant_proteins_for_nsc_index],secretome.data$x[constant_proteins_for_nsc_index,]), "./results/nsc/removed_proteins_by_NSC.csv",sep=",")

# esse ranking se refere ao conjunto de treinamento
ranking = as.matrix(match(complete_genes_names_list, secretome.data$genenames))
ranking <- rbind(ranking, as.matrix(constant_proteins_for_nsc_index))
write.matrix(ranking, file = "./results/nsc/rank_index_nsc.csv", sep = ",")




# final_results <- cbind(thres,probs_final,double_cross_error,bestFinalN)
finalN = nrow(list_genes_final)
write(finalN, file="./results/nsc/double_selected_N_nsc.txt")

final_inner_model_rep = selected_inner_model_final[[index]]
final_outer_model_rep = selected_outer_model_final[[index]]

final_pred_list = preds_outer_rep[[index]]
final_ref_list = refs_outer_rep[[index]]



#imprime todos os erros e valores de N
allpossibleN = numeric(0)
for( countn in 1:repetition){
  auxlist=pamr.listgenes(secretome.train2, secretome.data, thres[countn])
  allpossibleN <- append(allpossibleN, nrow(auxlist))
}
allNThresErrors <- cbind(thres, double_cross_error)
allNThresErrors <- cbind(allNThresErrors,allpossibleN)
colnames(allNThresErrors)<-c("thresholds","double_cross_error","N")
write.csv(allNThresErrors,"./results/nsc/errors_thresholds_and_Ns_nsc.csv")



library(caret)
pred  = as.vector(unlist( final_pred_list ))
ref = as.vector(unlist( final_ref_list ))

sink("./results/nsc/caret_nsc.txt")
cat("Average accuracy of double cross-validation repetitions: ")
cat(mean_double_cross_validation)
cat("\n\n")
cat("Maximum accuracy between repetitions: ")
cat(max_double_cross_validation)
cat("\n\n")
cat("Minimum accuracy between repetitions: ")
cat(min_double_cross_validation)
cat("\n\n")
cat("Best model index: ")
cat(index)
cat("\n\n")
cat("Selected final N: ")


cat(finalN)
cat("\n\n")
cat("Double-Cross-Validation Accuracy: ")    
cat(double_cross_error[index])
cat("\n\n")

pred=factor(pred, levels=unique(secretome.data$y))
ref=factor(ref, levels=unique(secretome.data$y))

confusionMatrix(pred, ref, positive=NULL)
pred
ref
sink()

pdf("./results/nsc/class_nsc_plot.pdf")
pamr.plotcen(secretome.train2, secretome.data, thres[index])
dev.off()

final_genes = pamr.listgenes(secretome.train2, secretome.data, thres[index], fitcv=NULL, genenames=TRUE)
write.csv(final_genes,"./results/nsc/gene_ranking_nsc_cutoff.csv")

png(filename="./results/nsc/gene_plot_nsc.png", width=4000, height=4000, pointsize=12)
#png(filename="gene_plot_NSC.png", width=1000, height=1000, pointsize=20)
pamr.geneplot(secretome.train2, secretome.data, 1)
dev.off()






# simple crossvalidation
# 
secretome.train <- pamr.train(secretome.data)
secretome.scales <- pamr.adaptthresh(secretome.train)
secretome.train2 <- pamr.train(secretome.data, threshold.scale=secretome.scales,n.threshold=thresholds_number, scale.sd=TRUE)  


cv_result = pamr.cv(secretome.train2, secretome.data)
pdf("./results/nsc/simple_cross-validation_nsc_plot.pdf")
pamr.plotcv(cv_result)
dev.off()


