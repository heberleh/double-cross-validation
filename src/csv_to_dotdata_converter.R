

#input example
#lines are samples
#columns are proteins

#check the last column of the resulted .data file. It must be the classes of your samples.
#The "class" attribute must not be listed at the "attributes list" line of the .data file

# IPI00465248_ENO1  IPI00018274_EGFR	IPI00021439_ACTB	IPI00023673_LGALS3BP	IPI00018219_TGFBI
# A2058_1	87.63	0	125.83	0	0
# A2058_2	64.94	0	99.98	0	0
# A2058_3	83.51	0	153.7	0	0
# Skmel_1	85.78	0	53.24	0	0
# Skmel_2	81.18	0	47.55	0	0
# Skmel_3	87.31	0.25	65.8	0	0
# A431_1	47.63	1	32.5	0	1
# A431_2	44.94	2	33.5	0	1
# A431_3	43.51	2	45.6	0	1

inputname = "dataset/Teste/TESTE_csv_to_data.csv"
outputname = "dataset/Teste/teste_zscore.data"

csv = read.csv(file=inputname, header=TRUE)

csv[1:3,]
size = nrow(csv)
n_attributes = ncol(csv) -2
attributes_names = colnames(csv)



#normalização
endcol = ncol(csv)-1
begincol = 2
intercol =begincol:endcol

endrow = nrow(csv)
beginrow = 1
interrow = beginrow:endrow

#for(i in intercol) csv[,i] = (csv[,i] - mean(csv[,i]))
#for(i in interrow) csv[i,] = (csv[i,intercol] - colMeans(csv[i,intercol]))
#for(i in 1:ncol(x))x[,i] = (x[,i]-mean(x[,i]))/sd(x[,i])
#for(i in 1:nrow(x))x[i,] = (x[i,]-mean(x[i,]))/sd(x[i,])


#csv[,intercol] = csv[,intercol] - rowMeans(csv[,intercol])
#csv[interrow,intercol] = csv[interrow,intercol] - colMeans(csv[interrow,intercol])


#normaliza por linha
#library(matrixStats)
#sdrow = transform(csv, SD=rowSds(csv[interrow,intercol], na.rm=TRUE))$SD
#csv[interrow,intercol] = (csv[interrow,intercol] - rowMeans(csv[interrow,intercol])) / sdrow

#normaliza por coluna
for(i in intercol) csv[,i] = (csv[,i]-mean(csv[,i]))/sd(csv[,i])


cat("DY",file=outputname)
cat("\n",file=outputname,append=TRUE)

cat(size,file=outputname,append=TRUE)
cat("\n",file=outputname,append=TRUE)

cat(n_attributes,file=outputname,append=TRUE)
cat("\n",file=outputname,append=TRUE)

cat(attributes_names[2:(n_attributes+1)],file=outputname,append=TRUE,sep=";")
cat("\n",file=outputname,append=TRUE)


for(i in 1:nrow(csv)){
  cat(as.character(csv[i,1]),file=outputname,append=TRUE)  
  
  for(j in 2:ncol(csv)){
    cat(";",file=outputname,append=TRUE)
    cat(as.character(csv[i,j]), file=outputname, append=TRUE)    
  }
    
   cat("\n",file=outputname,append=TRUE)
}


