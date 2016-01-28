library(caret)
library(pamr)

secretome.data <- pamr.from.excel("./dataset/spectral_counts_no_zeros_input.txt", 20, sample.labels=TRUE)
secretome.train <- pamr.train(secretome.data)

pred = pamr.predict(fit=secretome.train, newx=secretome.data$x, threshold=10)
ref = secretome.data$y

sink("out.txt")
confusionMatrix(pred, ref, positive=NULL)
sink()