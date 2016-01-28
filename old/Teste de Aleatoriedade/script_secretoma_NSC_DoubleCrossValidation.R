library("MASS")
library("pamr")

library(parallel)

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

secretome.data <- pamr.from.excel("./spectral_counts_no_zeros_input.txt", 20, sample.labels=TRUE)

original_complete_data <- secretome.data
# remove proteínas constantes = 0
#nonzeros <- rowSums(abs(secretome.data$x))>0
#secretome.data$x <- secretome.data$x[nonzeros,]
#secretome.data$genenames <- secretome.data$genenames[nonzeros]
#secretome.data$geneid <- secretome.data$geneid[nonzeros]

#scale rows with mean zero and standard deviation one
#for(i in 1:nrow(secretome.data$x))secretome.data$x[i,] = (secretome.data$x[i,]-mean(secretome.data$x[i,]))/sd(secretome.data$x[i,])

# remove o conjunto de teste após normalizar
test_index = c(3,11,13)
secretome.data$x <- secretome.data$x[,-test_index]
secretome.data$y <- secretome.data$y[-test_index]
secretome.data$samplelabels <- secretome.data$samplelabels[-test_index]

#nonzeros <- rowSums(abs(secretome.data$x))>0

#bestThresholdOuterLoop = vector()
#errorFold = vector()
bestFinal = list()
probs=vector()
thres=vector()
bestFinalNinner=vector()
bestFinalN=vector()
double_cross_error = vector()
repetition = 100
folds_list = list()
probs_final = vector()
thresholds_number = 90
inner_inner_rep = list()
outer_rep = list()
error_outers_rep  = list()
probabilities_outer_rep = list()
selected_inner_model_final = list()
selected_outer_model_final = list()
preds_outer_rep = list()
refs_outer_rep = list()
predictmany_rep = list()



secretome.train <- pamr.train(secretome.data,scale.sd=TRUE)
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
          secretome2_inner.train <- pamr.train(secretome2_inner.data,scale.sd=TRUE)
          secretome2_inner.scales <- pamr.adaptthresh(secretome2_inner.train)
          secretome2_inner.train2 <- pamr.train(secretome2_inner.data, threshold.scale=secretome2_inner.scales,n.threshold=thresholds_number,scale.sd=TRUE)
          
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
    secretome2.train <- pamr.train(secretome2.data,scale.sd=TRUE)
    secretome2.scales <- pamr.adaptthresh(secretome2.train)
    secretome2.train2 <- pamr.train(secretome2.data, threshold.scale=secretome2.scales,n.threshold=thresholds_number,scale.sd=TRUE)
    
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
  secretome.train <- pamr.train(secretome.data,scale.sd=TRUE)
  secretome.scales <- pamr.adaptthresh(secretome.train)
  
  selected_inner_model_final[[rep]] = selected_inner_model[[bestFinal[[rep]]$index]]
  selected_outer_model_final[[rep]] = outer_model[[bestFinal[[rep]]$index]]
  
  secretome.train2 <- pamr.train(secretome.data, threshold.scale=secretome.scales,scale.sd=TRUE)  
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


#Entre todas as repeti??es do Double-Cross-Validation, selecione a cujo erro mais se aproxima da m?dia
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

secretome.train <- pamr.train(secretome.data,scale.sd=TRUE)
secretome.scales <- pamr.adaptthresh(secretome.train)
secretome.train2 <- pamr.train(secretome.data, threshold.scale=secretome.scales,n.threshold=thresholds_number, scale.sd=TRUE)  

list_genes_final=pamr.listgenes(secretome.train2, secretome.data, thres[index])
complete_genes_list = pamr.listgenes(secretome.train2, secretome.data, 0.000000000000000000000000000000000000000000000000000000000000, genenames=TRUE)
write.csv(complete_genes_list,"complete_genes_list.csv")
complete_genes_names_list = complete_genes_list[,2]

# identifica proteínas descartadas para por no fim do ranking
constant_proteins_for_nsc <- as.matrix(secretome.data$genenames[!(secretome.data$genenames %in% complete_genes_names_list )])
constant_proteins_for_nsc_index <- match(constant_proteins_for_nsc,secretome.data$genenames)

write.matrix(cbind(original_complete_data$genenames[constant_proteins_for_nsc_index],original_complete_data$x[constant_proteins_for_nsc_index,]), "removed_proteins_by_NSC.csv",sep=",")

# esse ranking se refere ao conjunto de treinamento
ranking = as.matrix(match(complete_genes_names_list, secretome.data$genenames))
ranking <- rbind(ranking, as.matrix(constant_proteins_for_nsc_index))
write.matrix(ranking, file = "rank_index_nsc.csv", sep = ",")

# final_results <- cbind(thres,probs_final,double_cross_error,bestFinalN)
finalN = nrow(list_genes_final)
write(finalN, file="double_selected_N_nsc.txt")

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
write.csv(allNThresErrors,"errors_thresholds_and_Ns_NSC.csv")


library(caret)
pred  = as.vector(unlist( final_pred_list ))
ref = as.vector(unlist( final_ref_list ))

sink("caret_NSC_result.txt")
cat("Media dos erros dos doubles das 100 repeticoes: ")
cat(mean_double_cross_validation)
cat("\n\n")
cat("Erro maximo entre repeticoes: ")
cat(max_double_cross_validation)
cat("\n\n")
cat("Erro maximo entre repeticoes: ")
cat(min_double_cross_validation)
cat("\n\n")
cat("Index de todos os modelos cuja diferenca de seu erro e a media minima ")
cat(index)
cat("\n\n")
cat("N final escolhido: ")


cat(finalN)
cat("\n\n")
cat("Double-Cross-Validation Error: ")    
cat(double_cross_error[index])
cat("\n\n")

confusionMatrix(pred, ref, positive=NULL)
sink()


pamr.plotcen(secretome.train2, secretome.data, thres[index])

final_genes = pamr.listgenes(secretome.train2, secretome.data, thres[index], fitcv=NULL, genenames=TRUE)
write.csv(final_genes,"final_gene_listNSC.csv")

png(filename="gene_plot_NSC.png", width=4000, height=4000, pointsize=1/500)
pamr.geneplot(secretome.train2, secretome.data, 1)
dev.off()



#final_values = getLowerErrorMaxProbabilityMaxThreshold(double_cross_error,thres,probs_final)

#final_inner_model_rep = selected_inner_model_final[[final_values$index]]
#final_outer_model_rep = selected_outer_model_final[[final_values$index]]

#final_pred_list = preds_outer_rep[[final_values$index]]
#final_ref_list = refs_outer_rep[[final_values$index]]
# library(caret)
# pred  = as.vector(unlist( final_pred_list ))
# ref = as.vector(unlist( final_ref_list ))
# 
# sink("outers_pred_caret.txt")
# confusionMatrix(pred, ref, positive=NULL)
# sink()





# 
# secretome.train <- pamr.train(secretome.data)
# secretome.scales <- pamr.adaptthresh(secretome.train)
# secretome.train2 <- pamr.train(secretome.data, threshold.scale=secretome.scales,n.threshold=thresholds_number, scale.sd=TRUE)  
# 
# list_genes_final=pamr.listgenes(secretome.train2, secretome.data, final_values$higher_threshold)
# # final_results <- cbind(thres,probs_final,double_cross_error,bestFinalN)
# finalN = nrow(list_genes_final)
# 
# 
# thresholds <- final_inner_model_rep$train$threshold
# 
# train_inner <- final_inner_model_rep$train
# data_inner <- final_inner_model_rep$data
# test_x_inner <- final_inner_model_rep$testX
# test_y_inner <- final_inner_model_rep$testY
# size = length(test_y_inner)
# predicted_inner = list()
# predicted_prob_inner = list()
# error_inner = vector()
# N_inner = vector()
# for (i in 1:length(thresholds)){
#   threshold = thresholds[i]
#   predicted_inner[[i]] <- pamr.predict(fit=train_inner, newx=test_x_inner, threshold=threshold)
#   predicted_prob_inner[[i]] <- pamr.predict(fit=train_inner, newx=test_x_inner, threshold=threshold, type="posterior")
#   error_inner[i] = length(which(!as.logical(predicted_inner[[i]] == test_y_inner)))/size
#   N_inner[i] = tryCatch(nrow(pamr.listgenes(fit=train_inner, data=data_inner, threshold=threshold)),error=function(e) 0)
# }
# 
# 
# 
# 
# 
# 
# train_outer <- final_outer_model_rep$train
# data_outer <- final_outer_model_rep$data
# test_x_outer <- final_outer_model_rep$testX
# test_y_outer <- final_outer_model_rep$testY
# size = length(test_y_outer)
# predicted_outer = list()
# predicted_prob_outer = list()
# error_outer = vector()
# N_outer = vector()
# for (i in 1:length(thresholds)){
#   threshold = thresholds[i]
#   predicted_outer[[i]] <- pamr.predict(fit=train_outer, newx=test_x_outer, threshold=threshold)
#   predicted_prob_outer[[i]] <- pamr.predict(fit=train_outer, newx=test_x_outer, threshold=threshold, type="posterior")
#   
#   error_outer[i] = length(which(!as.logical(predicted_outer[[i]] == test_y_outer)))/size
#   N_outer[i] = tryCatch(nrow(pamr.listgenes(fit=train_outer, data=data_outer, threshold=threshold)),error=function(e) 0)
# }
# 
# 
# 
# library(caret)
# pred  = pamr.predict(final_outer_model_rep$train,secretome.data$x,threshold=final_values$higher_threshold)
# ref = secretome.data$y
# 
# sink("out.txt")
# confusionMatrix(pred, ref, positive=NULL)
# sink()
# 
# 
# 
# 
# 
# 
# 
# train <- secretome.train2
# data <- secretome.data
# test_x <- secretome.data$x
# test_y <- secretome.data$y
# size = length(test_y)
# predicted = list()
# predicted_prob = list()
# error = vector()
# N = vector()
# for (i in 1:length(thresholds)){
#   threshold = thresholds[i]
#   predicted[[i]] <- pamr.predict(fit=train, newx=test_x, threshold=threshold)
#   predicted_prob[[i]] <- pamr.predict(fit=train, newx=test_x, threshold=threshold, type="posterior")
#   error[i] = length(which(!as.logical(predicted[[i]] == test_y)))/size
#   
#   #tryCatch(pJohnson(.18, parms), error=function(e) alternativeFunction())
#   N[i] = tryCatch(nrow(pamr.listgenes(fit=train, data=data, threshold=threshold)),error=function(e) 0)  
# }
# 
# inner_outer_table <- cbind(thresholds,error_inner,N_inner,error_outer,N_outer,error,N)
# write.csv(inner_outer_table,"inner_outer_thresholds_errors_Ns_90thresholds.csv")
# 
# 
# path = "inner_outer_thresholds_erros_Ns_90thresholds_plot.pdf"
# pdf(path)
# str="Models' Errors by Thresholds"
# x=thresholds
# y1=error_inner
# y2=error_outer
# y3=error
# plot(x, y1, type="n", main=str)
# lines(x, y1, type="l", col="blue")
# lines(x, y2, type="l", col="red")
# lines(x, y3, type="l", col="green")
# axis(1,at=final_values$higher_threshold,labels=c("t"))
# # Legend Example  
# legend("bottomright", inset=.005, title="Intersection",
#        c("inner","outer","full"), fill=c("blue","red","green"), horiz=TRUE,cex=0.6)
# dev.off()
# 
# 
# genes = pamr.listgenes(secretome.train2,secretome.data,threshold=final_values$higher_threshold,genenames=TRUE)
# original_genes = secretome.data$genenames 
# index_centroids = which(original_genes %in% genes)
# delta = secretome.train2$centroids - secretome.train2$centroid.overall
# delta = cbind(secretome.data$genenames,delta)
# centroids = delta[index_centroids,]
# write.csv(delta,"d_completo.csv")
# write.csv(centroids,"d_centroids.csv")
# 
# pamr.plotcen()

# pdf("./analysis/plot_6fold.pdf")
# cv1 <- pamr.cv(fit=secretome.train2,data=secretome.data,folds=folds_list[[1]])
# pamr.plotcv(cv1)
# dev.off()


#pamr.plotcvprob(secretome.train2,secretome.data,2.509043845)
#2.509043845 6-fold
#3.507413642 loo




  
# 

#   
#   
# n_from_inner = bestFinalNinner[match(maxmax$higher_threshold,thres)]
# 
# 
# #data[sort.list(data[,2]), ]
# final_results_sortbythres <- final_results[sort.list(final_results[,1]),]
# final_results_sortbyprobs <- final_results[sort.list(final_results[,2]),]
# final_results_sortbyerror <- final_results[sort.list(final_results[,3]),]
# final_results_sortbyn <- final_results[sort.list(final_results[,4]),]
# 
# # path = paste("./analysis/prob_vs_threshold_kfold_",repetition,".pdf",sep="")
# # pdf(path)
# # ## Plot the class centroids   format(x, digits=0
# # str=paste("Higher thre:",format(maxmax$higher_threshold,nsmall=3),". maxProb:", format(maxmax$probability, nsmall=3), "N=", nrow(list_genes_final),sep=" ")
# # plot(final_results_sortbythres[,1], final_results_sortbythres[,2], type="n", main=str)
# # lines(final_results_sortbythres[,1], final_results_sortbythres[,2], type="l")
# # dev.off()
# 
# 
# 
# 
# str=paste("Higher thre:",format(maxmax$higher_threshold,nsmall=3),". error:", format(double_cross_error[match(maxmax$higher_threshold,thres)], nsmall=3), "N=", nrow(list_genes_final), "innerN=",n_from_inner,sep=" ")
# 
# path = paste("./analysis/threshold_vs_error_kfold_",repetition,".pdf",sep="")
# make_plot(path,final_results_sortbythres[,1], final_results_sortbythres[,3],str)
# 
# path = paste("./analysis/probs_vs_error_kfold_",repetition,".pdf",sep="")
# make_plot(path,final_results_sortbyprobs[,2], final_results_sortbyprobs[,3],str)
# 
# path = paste("./analysis/prob_vs_threshold_kfold_",repetition,".pdf",sep="")
# make_plot(path,final_results_sortbyprobs[,2], final_results_sortbyprobs[,1],str)
# 
# path = paste("./analysis/error_vs_probs_kfold_",repetition,".pdf",sep="")
# make_plot(path,final_results_sortbyerror[,3],final_results_sortbyerror[,2],str)
# 
# path = paste("./analysis/error_vs_thres_kfold_",repetition,".pdf",sep="")
# make_plot(path,final_results_sortbyerror[,3],final_results_sortbyerror[,1],str)
# 
# path = paste("./analysis/n_vs_error_kfold_",repetition,".pdf",sep="")
# make_plot(path,final_results_sortbyn[,5],final_results_sortbyn[,3],str)
# 
# path = paste("./analysis/n_vs_probs_kfold_",repetition,".pdf",sep="")
# make_plot(path,final_results_sortbyn[,5],final_results_sortbyn[,2],str)
# 
# 
# path = paste("./analysis/n_vs_thres_kfold_",repetition,".pdf",sep="")
# make_plot(path,final_results_sortbyn[,5],final_results_sortbyn[,1],str)
# 
# 
# csv_matrix <- final_results_sortbyprobs
# colnames(csv_matrix) <- c("Threshold","MaxMeanProbability","MeanKFoldError","innerN","N")
# path = paste("./analysis/results_thres_probs_erros_ns_",repetition,".csv",sep="")
# write.csv(csv_matrix,path)
# 
# list_genes_threshold = 0 #use it to get all genes
# list_genes_threshold = maxmax$higher_threshold
# #print each gene list from each k-train-model for best threshold from inner-loop
# text = pamr.listgenes(secretome.train2, secretome.data,  list_genes_threshold, genenames=TRUE)
# path = paste("./analysis/final_gene_list_",repetition,".csv",sep="")
# write.csv(text,path)
# 
# 
# path = paste("./analysis/density_N_",repetition,".pdf",sep="")
# pdf(path)
# d<-density(bestFinalN)
# str= paste("Best N = ",nrow(list_genes_final),". Highest density(h) =",format(d$x[match(max(d$y),d$y)],nsmall=1,digits=2), sep="")
# plot(d, main=str)
# axis(1,at=nrow(list_genes_final),labels=c("N"))
# axis(1,at=d$x[match(max(d$y),d$y)],labels=c("h"))
# dev.off()














# path = paste("./analysis/gene_plot_",repetition,".pdf",sep="")
# pdf(path)
# #pamr.plotcen(secretome.train2, secretome.data, threshold=maxmax$higher_threshold)
# pamr.plotcen(secretome.train2, secretome.data, threshold=1.936)
# dev.off()




#print each gene list from each k-train-model for best threshold from inner-loop
#text = pamr.listgenes(secretome.train2, secretome.data,  3.489524, genenames=TRUE)
#path = "./final_list/final_gene_list.csv"
#write.csv(text,path)






# newthres=sort(thres)
# unique(thres)
# # Kernel Density Plot
# d <- density(thres) # returns the density data
# plot(intervals=d,) # plots the results 



# plot (c(1968,2010),c(0,10),type="n", # sets the x and y axes scales
#       
#       xlab="Year",ylab="Expenditures/GDP (%)") # adds titles to the axes
# 
# lines(year,defense,col="red",lwd=2.5) # adds a line for defense expenditures
# 
# lines(year,health,col="blue",lwd=2.5) # adds a line for health expenditures
# 
# legend(2000,9.5, # places a legend at the appropriate place c("Health","Defense"), # puts text in the legend
#        
#        lty=c(1,1), # gives the legend appropriate symbols (lines)
#        
#        lwd=c(2.5,2.5),col=c("blue","red")) # gives the legend lines the correct color and width


# segments(x, y-sd,x, y+sd)
# epsilon = 0.02
# segments(x-epsilon,y-sd,x+epsilon,y-sd)
# segments(x-epsilon,y+sd,x+epsilon,y+sd)







#utilizado para fazer media de probs por threshold...
#   predictmany_rep[[rep]] = predictmany

#   mean_error_by_threhsold = rep(0, length(fixed_thresholds))
#   mean_prob_by_threshold = rep(0,length(fixed_thresholds))

#   for (o in 1:1){
#     predict = predictmany[[o]]
#     probs = predict$prob
#     pred = predict$predclass
#     
#     for (thre in 1:length(fixed_thresholds)){
#       current_thre = fixed_thresholds[thre]
#       current_pred = unlist(pred[,thre])
#       current_ref = unlist(ref[o])
#       m = vector()
#       current_prob = probs[,,thre]
#       for (l in 1:nrow(current_prob)){
#         m[l]= max(current_prob[,l])
#       }
#       mean_prob_by_threshold[thre] =  mean(m)
#       
#       mean_error_by_threhsold[thre] = length(which(!as.logical(current_pred == current_ref)))/length(current_ref)     
#     }
#   }
#   
#   for (o in 2:nfold){
#     predict = predictmany[[o]]
#     probs = predict$prob
#     pred = predict$predclass  
#     for (thre in 1:length(fixed_thresholds)){
#       current_thre = fixed_thresholds[thre]
#       current_pred = unlist(pred[,thre])
#       current_ref = unlist(ref[o])
#       m = vector()
#       current_prob = probs[,,thre]
#       for (l in 1:nrow(current_prob)){
#         m[l]= max(current_prob[l,])
#       }
#       mean_prob_by_threshold[thre] = mean(c(mean_prob_by_threshold[thre], mean(m)))
#       
#       mean_error_by_threhsold[thre] = mean(c(mean_error_by_threhsold[thre], (length(which(!as.logical(current_pred == current_ref)))/length(current_ref))))
#     }
#   }
#   

