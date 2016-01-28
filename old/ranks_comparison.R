

library(sets)
library(parallel)


db <- read.table("ranks_comparison.csv", header=TRUE, sep=",")


sizeABC = vector()
sizeAB = vector()
sizeAC = vector()
sizeBC = vector()

exA = vector()
exB = vector()
exC = vector()

max_size = nrow(db)
#max_size = 601

for(i in 1:max_size){
  N <- i
  
  NSC<-db[1:N,1] #A
  SVMRFE<-db[1:N,2] #B
  Beta<-db[1:N,3] #C  
  
  intersectionABC <- intersect(SVMRFE,Beta)
  intersectionABC <- intersect(intersectionABC,NSC)
  
  intersectionAB <- intersect(NSC,SVMRFE)
  intersectionAC <- intersect(NSC,Beta)
  intersectionBC <- intersect(SVMRFE,Beta) 
    
  unionABC <- union(NSC,SVMRFE)
  unionABC <- union(unionABC,Beta)
  
  unionAB <- union(NSC, SVMRFE)
  unionAC <- union(NSC,Beta)
  unionBC <- union(SVMRFE,Beta)  
  
  uBC <- as.set(unionBC)
  uAC <- as.set(unionAC)
  uAB <- as.set(unionAB)
  sNSC <- as.set(as.character(NSC))
  sSVMRFE <- as.set(as.character(SVMRFE))
  sBeta <- as.set(as.character(Beta))
  
  exA[N] = length(set_complement(uBC, sNSC))
  exB[N] = length(set_complement(uAC, sSVMRFE))
  exC[N] = length(set_complement(uAB, sBeta))
  
  
  sizeABC[N] = length(intersectionABC)/length(unionABC)
  sizeAB[N] = length(intersectionAB)/length(unionAB) 
  sizeAC[N] = length(intersectionAC)/length(unionAC)
  sizeBC[N] = length(intersectionBC)/length(unionBC)  
}



xlim = c(1,1600)
ylim = c(0,1)


x = 1:max_size
y = sizeABC
str = "Similarities between sets of size N"
pdf("Jaccard_similarity_coefficient__3juntos.pdf")
plot(x, y, xlim=xlim,ylim=ylim, xlab="N", ylab="Jaccard similarity coefficient", type="n", main=str) #log="x"
lines(x, y, type="l", col="black")    
legend("bottomright", inset=.005, title="Similarity between sets:",
       c("NSC, SVM-RFE, Beta-Binomial"), fill=c("black"), horiz=TRUE,cex=0.75)
dev.off()

x = 1:max_size
y1 = sizeAB
y2 = sizeAC
y3 = sizeBC
str = "Similarities between sets of size N"
pdf("Jaccard_similarity_coefficient__2a2.pdf")
plot(x, y1, xlim=xlim,ylim=ylim,  xlab="N", ylab="Jaccard similarity coefficient", type="n", main=str)
lines(x, y1, type="l", col="red")    
lines(x, y2, type="l", col="black")
lines(x, y3, type="l", col="blue")
legend("bottomright", inset=.005, title="Similarity between sets:",
       c("NSC, SVM-RFE", "NSC, Beta-Binomial", "SVM-RFE, Beta-Binomial"), fill=c("red","black","blue"), horiz=TRUE,cex=0.6)
dev.off()


x = 1:max_size
y1 = sizeAB
y2 = sizeAC
y3 = sizeBC
y4 = sizeABC
str = "Similarities between sets of size N"
pdf("Jaccard_similarity_coefficient__2a2e3.pdf")
plot(x, y1, xlim=xlim,ylim=ylim, xlab="N", ylab="Jaccard similarity coefficient", type="n", main=str)
lines(x, y1, type="l", col="red")    
lines(x, y2, type="l", col="black")
lines(x, y3, type="l", col="blue")
lines(x, y4, type="l", col="gray")
legend("bottomright", inset=.005, title="Similarity between sets:",
       c("NSC, SVM-RFE", "NSC, Beta-Binomial", "SVM-RFE, Beta-Binomial","All three sets"), fill=c("red","black","blue","gray"), horiz=TRUE,cex=0.5)
dev.off()


x = 1:max_size
y1 = exA
y2 = exB
y3 = exC
m1 = max(exA)
m2 = max(exB)
m3 = max(exC)
max = max(m1,m2,m3)
str = "Disjoint regions of sets of size N"
pdf("Disjoint_regions.pdf")
plot(x, y1, xlim=xlim, ylim=c(0,max)  xlab="N", ylab="Disjoint region size", type="n", main=str)
lines(x, y1, type="l", col="red")    
lines(x, y2, type="l", col="black")
lines(x, y3, type="l", col="blue")
legend("bottomright", inset=.005, title="Sets:",
       c("NSC", "SVM-RFE", "Beta-Binomial"), fill=c("red","black","blue"), horiz=TRUE,cex=0.6)
dev.off()


m <- cbind(sizeAB,sizeAC,sizeBC,sizeABC,exA,exB,exC)
colnames(m) <- c("NSC,SVMRFE","NSC,Beta","SVMRFE,Beta","NSC,SVMRFE,Beta","NSC exclusive", "SVM-RFE exclusive","Beta-Binomial exclusive")
write.csv(m,"./rank intersections analysis/ranks_intersections_list.csv")






