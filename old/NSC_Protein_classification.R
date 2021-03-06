

#csv = read.csv(file="NSC_gene_list.csv")
csv = read.csv(file="NSC_complete_genes_list.csv")

names = colnames(csv)

data = csv[,4:6]

mean_melanoma = data[,1]
mean_carcinoma = data[,2]
mean_normal = data[,3]

classes = vector()

get_class <- function(index){ 
  if (index == 1) return ("Melanoma");
  if (index == 2) return ("Carcinoma");
  if (index == 3) return ("Normal"); 
  if (index == 0) return ("????????")
}

get_border <- function(mean, border){
  dif1 = abs(border[1]-mean)
  dif2 = abs(border[2]-mean)
  
  if (dif1 < dif2){
    return (1)
  }else{
    if (dif2 < dif1){
      return (2)
    }else{
      return (0)
    }
  }  
}

each_border_l = vector()
for (i in 1:nrow(data)){
  
  v = c(mean_melanoma[i], mean_carcinoma[i], mean_normal[i])
  
  indexs = sort(v, index.return = TRUE)$ix
  
  border = c((v[indexs[1]] + v[indexs[2]])/2, (v[indexs[2]] + v[indexs[3]])/2)
  
  each_border = c(get_border(v[1], border), get_border(v[2], border), get_border(v[3], border))
  each_border_l[i] = paste(each_border,collapse=" ")
  
  counts1 = length(which(each_border == 1))
  counts2 = length(which(each_border == 2))
  
  if (counts1 < counts2){
    index = which(each_border == 1)
  }else{
    if (counts2 < counts1){
      index = which(each_border == 2)
    }else{
      index = 0
    }
  }
  
  classes[i] = get_class(index)
  
}



newcsv = cbind(csv,classes)
newcsv = cbind(newcsv,each_border_l)

colnames(newcsv) <- c(names, c("Class", "Borders"))
#write.csv(newcsv,"NSC_gene_list_CLASSIFIED.csv")
write.csv(newcsv,"NSC_complete_gene_list_CLASSIFIED.csv")



