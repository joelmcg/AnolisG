#Load required packages

library(ape)
library(evolqg)
library(mgcv)


#Load G matrices

crisG<-as.matrix(read.table("data/cris.txt"))
everG<-as.matrix(read.table("data/ever.txt"))
grahG<-as.matrix(read.table("data/grah.txt"))
lineG<-as.matrix(read.table("data/line.txt"))
pulcG<-as.matrix(read.table("data/pulc.txt"))
sagrG<-as.matrix(read.table("data/sagr.txt"))
smarG<-as.matrix(read.table("data/smar.txt"))



#Random Skewers

Gs<-list(crisG,everG,grahG,lineG,pulcG,sagrG,smarG)
rsG<-RandomSkewers(Gs,num.vectors=10000)
rsG

#Correlation between genetic distance and random skewer correlations

geneticdist <- c(0.64127367, 0.987323604, 0.987323604, 0.541900856, 0.987323604, 1, 0.987323604, 0.987323604, 0.64127367, 
                 0.987323604, 1, 0.744421662, 0.987323604, 0.894946327, 1, 0.987323604, 0.894946327, 1, 0.987323604, 1, 1)
skews <- rsG$correlations[lower.tri(rsG$correlations, diag = FALSE)]
cor(geneticdist, skews)

#Partial correlation with ecomorph, controlling for genetic distance

geneticdist6 <- c(0.64127367, 0.987323604, 0.987323604, 0.987323604, 1, 0.987323604, 0.987323604, 0.987323604, 1, 0.744421662, 
                  0.89496327, 1, 0.987323604, 1, 1)

eco.mat <- c(0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0)

skews6 <- rsG$correlations[c(1:4, 6:7), c(1:4, 6:7)]
skews6v <- skews6[lower.tri(skews6, diag = FALSE)]
pcordat <- as.data.frame(cbind(eco.mat, geneticdist6, skews6v))
pcor(pcordat)$estimate[1, 3]


#Read trees and data for phylogenetic analyses. Trees as newick format, data as csv with one line per species, listing element and SE (from ASReml) separately

rge<-as.data.frame(read.csv("rge.csv"))
a6tree<-read.tree("pruned3.tre")
anoletree<-read.tree("pruned7.tre")
anole.data<-as.data.frame(read.table("data/anoldat.txt"))
eco <- as.integer(anole.data$ecomorph)
names(eco) <- rownames(anole.data)

#Construct phylogenetic covariance matrix
V<-corMatrix(Initialize(corBrownian(phy=a6tree),anole.data))


#Phylogenetic correlations

phylr<-c()
for(j in 1:28){
  rvec<-c(as.numeric(rge[j,2]),as.numeric(rge[j,3]),as.numeric(rge[j,4]),as.numeric(rge[j,5]),as.numeric(rge[j,7]),as.numeric(rge[j,8]))
  names(rvec)<-c("cristatellus","evermanni","grahami","lineatopus","sagrei","smaragdinus")
  a.Y <- matrix(1,1,6) %*% solve(V) %*% rvec/sum(solve(V));
  a.X <- matrix(1,1,6) %*% solve(V) %*% eco/sum(solve(V));
  rvec.gls<-(rvec-c(a.Y)) %*% solve(V) %*% (eco-c(a.X))/sqrt(((rvec-c(a.Y)) %*% solve(V) %*% (rvec-c(a.Y)))*((eco-c(a.X)) %*% solve(V) %*% (eco-c(a.X))));
  phylr[j]<-as.numeric(rvec.gls)
}

phylovar <- c()
for (j in 1:8) {
  vvec <- c(crisG[j,j], everG[j,j], grahG[j,j], lineG[j,j], 
            sagrG[j,j], smarG[j,j])
  names(vvec) <- c("cristatellus", "evermanni", "grahami", "lineatopus", "sagrei", "smaragdinus")
  a.Y <- matrix(1, 1, 6) %*% invV %*% vvec/sum(invV)
  a.X <- matrix(1, 1, 6) %*% invV %*% eco/sum(invV)
  vvec.gls <- (vvec - c(a.Y)) %*% invV %*% (eco - c(a.X))/sqrt(((vvec - c(a.Y)) %*% invV %*% (vvec - 
                                                                                                c(a.Y))) * ((eco - c(a.X)) %*% invV %*% (eco - c(a.X))))
  phylovar[j] <- as.numeric(vvec.gls)
}



#Phylogenetic signal

Kr <- c()  
for (j in 1:28) {
  rvec<-c(as.numeric(rge[j,2]),as.numeric(rge[j,3]),as.numeric(rge[j,4]),as.numeric(rge[j,5]),as.numeric(rge[j,6]),as.numeric(rge[j,7]),as.numeric(rge[j,8]))
  names(rvec) <- c("cristatellus", "evermanni", "grahami", "lineatopus", "pulchellus", "sagrei", "smaragdinus")
  Kr[j] <- phylosig(anoletree, rvec)
}



Kvar <- c()
for (j in 1:8) {
  vvec <- c(crisG[j,j], everG[j,j], grahG[j,j], lineG[j,j], pulcG[j,j],
            sagrG[j,j], smarG[j,j])
  names(vvec) <- c("cristatellus", "evermanni", "grahami", "lineatopus", "pulchellus", "sagrei", "smaragdinus")
  Kvar[j] <- phylosig(anoletree, vvec)
}

