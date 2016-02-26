
#To install ibb package, please follow the instructions: http://www.oncoproteomics.nl/software/BetaBinomial.html
library(ibb)

d <- read.delim("./dataset/current/train_beta-binomial.txt", header = TRUE)

# using all available CPU cores -> n.threads=0


#melanoma,  carcinoma, normal
#out <- bb.test(d[, 3:20], colSums(d[, 3:20]), c(rep("melanoma", 6), rep("carcinoma", 6), rep("normal", 6)), n.threads = 0)

# four groups
out <- bb.test(d[, 3:36], colSums(d[, 3:36]), c(rep("12.5fmol", 8), rep("25fmol", 9), rep("50fmol", 8), rep("100fmol",9)), n.threads = 0)


# write result to file
merged_full <- cbind(d, out$p.value)
merged <-cbind(d[,1:2], out$p.value)

merged_full <- merged_full[order(merged_full[,(ncol(merged_full))]),]
merged <- merged[order(merged[,3]),]

write.table(merged, file = "./results/beta-binomial/beta-binomial.sorted.short.csv" , sep = ",", row.names = FALSE)
write.table(merged_full, file = "./results/beta-binomial/beta-binomial.sorted.csv" , sep = ",", row.names = FALSE)

#write.table(cbind(d, out$p.value), file = "./results/beta-binomial/beta-binomial.csv" , sep = ",", row.names = FALSE)
#write.table(cbind(d[,1:2], out$p.value), file = "./results/beta-binomial/beta-binomial.short.csv" , sep = ",", row.names = FALSE)
