library(ibb)

d <- read.delim("spectral_counts_no_zeros_input_BetaBinomial.txt", header = TRUE)

# using all available CPU cores
#melanoma,  carcinoma, normal
out <- bb.test(d[, 3:20], colSums(d[, 3:20]), c(rep("melanoma", 6), rep("carcinoma", 6), rep("normal", 6)), n.threads = 0)


# write result to file


merged_full <- cbind(d, out$p.value)
merged <-cbind(d[,1:2], out$p.value)

merged_full <- merged_full[order(merged_full[,(ncol(merged_full))]),]
merged <- merged[order(merged[,3]),]

write.table(merged, file = "beta-binomial_Total-out_ordenado_resumido.csv" , sep = ",", row.names = FALSE)
write.table(merged_full, file = "beta-binomial_Total-out_ordenado.csv" , sep = ",", row.names = FALSE)

write.table(cbind(d, out$p.value), file = "beta-binomial_Total-out.csv" , sep = ",", row.names = FALSE)
write.table(cbind(d[,1:2], out$p.value), file = "beta-binomial_Total-out_resumido.csv" , sep = ",", row.names = FALSE)
