





#     secretome2.train <- pamr.train(secretome2.data)
#     secretome2.scales <- pamr.adaptthresh(secretome2.train)
#     secretome2.train2 <- pamr.train(secretome2.data, threshold.scale=secretome2.scales)

#     secretome2.results_inner <- pamr.cv(fit = secretome2.train2, data = secretome2.data)      
#     
#     bestThreshold[i] = getBestThreshold(secretome2.results_inner)  
#     
#     
#     secretome2.results_outer = pamr.predict(secretome2.train2, secretome.data$x[,folds[[i]]], threshold=bestThreshold[i])
#   
#     true_prediction[[i]]=  secretome.data$y[folds[[i]]] == secretome2.results_outer
#     true_prediction_total[folds[[i]]] = secretome.data$y[folds[[i]]] == secretome2.results_outer
#     
#     secretome2 = list(data= secretome2.data,
#                       train=secretome2.train, 
#                       train2=secretome2.train2, 
#                       results_inner=secretome2.results_inner, 
#                       results_outer=secretome2.results_outer)  
#     
#     fold[[i]] <- secretome2  
#     
#     path = paste("error_inner_fold",i,"NSC.pdf",sep="_")
#     pdf(path)
#     pamr.plotcv(secretome2.results_inner)
#     dev.off()





bestFinalThreshold = finalThreshold(bestThresholdOuterLoop,errorFold)
meanThreshold = mean(as.numeric(bestThresholdOuterLoop))


#which threshold should I use?
#cThreshold = meanThreshold
cThreshold = bestFinalThreshold


secretome.train <- pamr.train(secretome.data)
secretome.scales <- pamr.adaptthresh(secretome.train)
secretome.train2 <- pamr.train(secretome.data, threshold.scale=secretome.scales)
secretome.results_inner <- pamr.cv(fit = secretome.train2, data = secretome.data)      

pamr.listgenes(secretome.train2, secretome.data,  2.966496, genenames=TRUE)

pdf("centroids_NSC_double_cross_validation.pdf")
## Plot the class centroids
pamr.plotcen(secretome.train2, secretome.data, threshold=cThreshold)
dev.off()


pdf("Generic_cross-validated_class_probabilities_by_class_fixed_threshold.pdf")
## Plot the cross-validated class probabilities by class
pamr.plotcvprob(secretome.results_inner, secretome.data, threshold=cThreshold)
dev.off()


pdf("the_most_significant_genes_NSC_double_cross_validation.pdf")
## Make a gene plot of the most significant genes
pamr.geneplot(secretome.train2, secretome.data, threshold=cThreshold) #número mínimo de genes
dev.off()


pdf("Generic_cross-validated_Missclassification_error_by_threshold.pdf")
pamr.plotcv(secretome.results_inner)
dev.off()

cat(paste("\nThreshold Médio:", meanThreshold, "\n",sep=" "))
cat(paste("\nThreshold Máximo:", bestFinalThreshold, "\n",sep=" "))
cat(paste("\nThreshold escolhido:", cThreshold, "\n", sep=" "))







bestThreshold

which(as.logical(true_prediction))

which(!as.logical(true_prediction))





t = pamr.train(data=secretome.data, sample.subset=1:12)

pamr.predict()


for (i in 1:30){
  cat(paste("max of threshold ",fold[[1]]$results$threshold[i],": ",sep=" "))
  cat(max(sum(fold[[1]]$results$prob[,,i][,1]),sum(fold[[1]]$results$prob[,,i][,2]),sum(fold[[1]]$results$prob[,,i][,3])))
  cat("\n")
}


k=3
for (i in 1:30){
  cat(paste("max of threshold ",fold[[k]]$results$threshold[i],": ",sep=" "))
  total = 0.0
  for (j in 1:15){
    total=sum(total,max(fold[[k]]$results$prob[,,i][j,1],fold[[k]]$results$prob[,,i][j,2],fold[[k]]$results$prob[,,i][j,3]))
  }
  cat(total)
  cat("\n")
}

# secretoma.train <- pamr.train(secretoma.data)
# secretoma.results <- pamr.cv(fit = secretoma.train, data = secretoma.data)


## Cross-validate the classifier
# secretoma.results <- pamr.cv(fit = secretoma.train, data = secretoma.data)

## Plot the cross-validated error curves
# pamr.plotcv(secretoma.results)




getBestThreshold <- function(results){
  min_error = 1
  threshold = results$threshold[(length(results$threshold))]
  for(ii in (length(results$threshold)):1){
    error = results$error[ii]
    if (error < min_error){
      min_error = error
      threshold = results$threshold[ii]
    }    
  }  
  return (threshold)
}

listBestThresholdsOuter <- function(threshold,trueClassified){  
  bestThreshold = max(threshold)
  error = list()
  thresholdList = list()
  for(j in 1:(length(trueClassified))){
    error[j] = length(which(!as.logical(trueClassified[[j]]))) / length(trueClassified[[j]])
    thresholdList[j] = threshold[j]
  }
  result = list(error = error,threshold = thresholdList)
  
  return (result)
}

getBestThresholdOuter <- function(threshold,trueClassified){  
  best = max(threshold)  
  min_error = 1
  for(j in 1:(length(trueClassified))){
    error = length(which(!as.logical(trueClassified[[j]]))) / length(trueClassified[[j]])
    if (error < min_error && threshold[j] > best){
      min_error = error
      best = threshold[j]
    }
  }
  return (best)
}

finalThreshold <- function(threshold,errorFold){
  best = max(threshold)  
  min_error = 1
  for(j in 1:(length(errorFold))){
    if (errorFold[j] < min_error && threshold[j] > best){
      min_error = errorFold[j]
      best = threshold[j]
    }
  }
  return (best)  
}





# argy <- y
# n.threshold <- length(threshold)
# yhat <- rep(list(y), n.threshold)
# names(yhat) <- paste(seq(n.threshold))
# yhat <- data.frame(yhat)
# n.class <- table(y)
# prob <- array(1, c(n, length(n.class), n.threshold))
# size <- double(n.threshold)
# hetero <- object$hetero
# cv.objects = vector("list", nfold)
# for (ii in 1:nfold) {
#   cat("Fold", ii, ":")
#   a <- nsc(x[, -folds[[ii]]], y = argy[-folds[[ii]]], 
#            x[, folds[[ii]], drop = FALSE], proby = proby[-folds[[ii]], 
#                                                          ], threshold = threshold, threshold.scale = threshold.scale, 
#            se.scale = se.scale, prior = prior, hetero = hetero, 
#            ..., remove.zeros = FALSE)
#   size <- size + a$nonzero
#   prob[folds[[ii]], , ] <- a$prob
#   yhat[folds[[ii]], ] <- a$yhat
#   cat("\n")
#   cv.objects[[ii]] = a
# }
# if (missing(object)) 
#   size <- round(size/nfold)
# else size <- object$nonzero
# error <- rep(NA, n.threshold)
# loglik <- error
# pvalue.survival <- error
# pvalue.survival.func <- function(group, survival.time, censoring.status, 
#                                  ngroup.survival) {
#   temp <- coxph(Surv(survival.time, censoring.status) ~ 
#                   as.factor(group))
#   loglik <- 2 * (temp$loglik[2] - temp$loglik[1])
#   return(1 - pchisq(loglik, ngroup.survival - 1))
# }
# if (!is.null(proby)) {
#   proby.temp <- proby
# }
# else if (!is.null(survival.time)) {
#   proby.temp <- pamr.surv.to.class2(survival.time, censoring.status, 
#                                     n.class = ngroup.survival)$prob
# }
# for (i in 1:n.threshold) {
#   if (is.null(survival.time) & is.null(proby)) {
#     error[i] <- sum(yhat[, i] != y)/n
#   }
#   if (!is.null(survival.time)) {
#     temp <- c(yhat[, i], names(table(y)))
#     Yhat <- model.matrix(~factor(temp) - 1, data = list(y = temp))
#     Yhat <- Yhat[1:length(yhat[[ii]]), ]
#     error[i] <- (length(yhat[, i]) - sum(Yhat * proby.temp))/n
#   }
#   if (is.null(survival.time)) {
#     loglik[i] <- sum(log(prob[, , i][cbind(seq(1, n), 
#                                            unclass(y))]))/n
#   }
#   if (!is.null(survival.time)) {
#     pvalue.survival[i] <- pvalue.survival.func(yhat[, 
#                                                     i], survival.time, censoring.status, ngroup.survival)
#   }
# }
# obj <- list(threshold = threshold, error = error, loglik = loglik, 
#             size = size, yhat = yhat, y = y, prob = prob, folds = folds, 
#             cv.objects = cv.objects, pvalue.survival = pvalue.survival, 
#             call = this.call)
# class(obj) <- "nsccv"
# obj






folds = 1:outer_n #implica treino de tamanho ~13a15 e teste de tamanho ~3a5

outer_k <- kfold(t(secretome.data$x),k=folds,secretome.data$y)

#outer loop
for (test_fold in folds){
  train_folds <- folds[-test_fold]
  outer_train.x <- x[which(outer_k %in% c(train_folds))]   # all(!= i)
  outer_train.y <- y[which(outer_k %in% c(train_folds))]
  outer_test.x <-  x[which(outer_k %in% c(test_fold))]   # all(== i)
  outer_test.y <-  y[which(outer_k %in% c(test_fold))]   # all(== i)
  
  inner_n = 6
  inner_folds = 1:inner_n #se outer =5, implica treino de tamanho ~11a12 e teste de tamanho ~2a3
  inner_k <- kfold(outer_train.x,k=inner_n,t(outer_train.y))
  
  #inner loop - 6-fold-cross-validation feito pela função pamr.cv()
  
  #calcula NSC
  #armazena accuracy
  #5 Ns
  
}
#5 Ns

}












secretoma.data <- pamr.from.excel("counts_antiga_sem_normalização_vertical_NSC.txt", 20, sample.labels=TRUE)

secretoma.train <- pamr.train(secretoma.data)


#ntries = 10, reduction.factor = 0.9
#new.scales <- pamr.adaptthresh(secretoma.train,ntries=10,reduction.factor = 0.9)
#secretoma.train2 <- pamr.train(secretoma.data, threshold.scale=new.scales)
#secretoma.results2 <- pamr.cv(secretoma.train2, secretoma.data)


## Cross-validate the classifier
secretoma.results <- pamr.cv(fit = secretoma.train, data = secretoma.data)

## Plot the cross-validated error curves
pamr.plotcv(secretoma.results)

## Compute the confusion matrix for a particular model (threshold=4.0) 

pamr.confusion(secretoma.results, threshold=2.000)

## Plot the cross-validated class probabilities by class
pamr.plotcvprob(secretoma.results, secretoma.data, threshold=2.00)

## Plot the class centroids
pamr.plotcen(secretoma.train, secretoma.data, threshold=2.00)

## Make a gene plot of the most significant genes
pamr.geneplot(secretoma.train, secretoma.data, threshold=2.00) #número mínimo de genes
pamr.geneplot(secretoma.train, secretoma.data, threshold=0.952) #número máximo de genes


# Estimate false discovery rates and plot them
fdr.obj<- pamr.fdr(secretoma.train, secretoma.data)

pamr.plotfdr(fdr.obj)

## List the significant genes
pamr.listgenes(secretoma.train, secretoma.data, threshold=2.00)

## Try heterogeneity analysis, with class "normal" taken to be the normal group
secretoma.train2 <- pamr.train(secretoma.data,hetero="normal")
secretoma.results2 <-  pamr.cv(secretoma.train2, secretoma.data)


## Look for better threshold scalings
secretoma.scales <- pamr.adaptthresh(secretoma.train)
secretoma.train3 <- pamr.train(secretoma.data, threshold.scale=secretoma.scales)
secretoma.results3 <-  pamr.cv(secretoma.train3, secretoma.data)

