# Load required packages

library(ape)
library(evolqg)
library(mgcv)
library(geomorph)
library(phytools)

# Load G matrices (from McGlothlin et al. 2018)

crisG<-as.matrix(read.table("data/cris.txt"))
everG<-as.matrix(read.table("data/ever.txt"))
grahG<-as.matrix(read.table("data/grah.txt"))
lineG<-as.matrix(read.table("data/line.txt"))
pulcG<-as.matrix(read.table("data/pulc.txt"))
sagrG<-as.matrix(read.table("data/sagr.txt"))
smarG<-as.matrix(read.table("data/smar.txt"))

# Calculation of random skewers correlations (rRS, figs. 2-3, table A2 cols 1-3)

Gs<-list(crisG, everG, grahG, lineG, pulcG, sagrG, smarG)
rsG<-RandomSkewers(Gs, num.vectors=10000)
rsG

# Read trees and data for phylogenetic analyses. Trees as newick format, 
# data as csv with one line per species, listing element and SE (from ASReml) 
# separately

rge<-as.data.frame(read.csv("data/rge.csv"))
a6tree<-read.tree("data/pruned3.tre")
anoletree<-read.tree("data/pruned7.tre")
anole.data<-as.data.frame(read.table("data/anoldat.txt"))
eco <- as.integer(anole.data$ecomorph)
names(eco) <- rownames(anole.data)

# Rough plot of genetic distance and rRS (fig. 3 plotted in GraphPad Prism)

geneticdist <- 41.454817 * c(0.64127367, 0.987323604, 0.987323604, 0.541900856, 
                             0.987323604, 1, 0.987323604, 0.987323604, 
                             0.64127367, 0.987323604, 1, 0.744421662, 
                             0.987323604, 0.894946327, 1, 0.987323604, 
                             0.894946327, 1, 0.987323604, 1, 1)
skews <- rsG$correlations[lower.tri(rsG$correlations, diag = FALSE)]
plot(geneticdist, skews)

# K_mult for rRS (text)

rsmatpcoa <- function(x, sp.names){
  skmat <-  matrix(rep(1,length(x)), nrow(x), ncol(x))
  skmat[lower.tri(skmat, diag = FALSE)] <- x[lower.tri(x, diag=FALSE)]
  skmat[upper.tri(skmat, diag = FALSE)] <- t(skmat)[upper.tri(skmat, 
                                                              diag = FALSE)]
  row.names(skmat) <- sp.names
  colnames(skmat) <- sp.names
  skdist <- sqrt(2*(1-skmat))
  pcoa(skdist)$vectors
}

species <- c("cristatellus", "evermanni", "grahami", "lineatopus", "pulchellus", 
             "sagrei", "smaragdinus")

physignal(rsmatpcoa(rsG$correlations, species), anoletree)

# Mantel with phylogenetically informed permuations for ecomorph effect (text)

species6 <- c("cristatellus", "evermanni", "grahami", "lineatopus", "sagrei", 
              "smaragdinus")
eco.mat <- c(0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0)
eco6mat <- matrix(rep(0,36), 6, 6)
eco6mat[lower.tri(eco6mat, diag = FALSE)] <- 1 - eco.mat
eco6mat[upper.tri(eco6mat, diag = FALSE)] <- t(eco6mat)[upper.tri(eco6mat, 
                                                                  diag = FALSE)]
colnames(eco6mat) <- species6
row.names(eco6mat) <- species6

PhyloMantel(a6tree, rs6dist, eco6mat)

# Construct phylogenetic covariance matrix for individual element analyses

V <- corMatrix(Initialize(corBrownian(phy=a6tree), anole.data))

# Phylogenetic signal (table A3)

Kr <- c()  
for (j in 1:28) {
  rvec<-c(as.numeric(rge[j,2]), as.numeric(rge[j,3]), as.numeric(rge[j,4]),
          as.numeric(rge[j,5]), as.numeric(rge[j,6]), as.numeric(rge[j,7]),
          as.numeric(rge[j,8]))
  names(rvec) <- species
  Kr[j] <- phylosig(anoletree, rvec)
}

Kvar <- c()
for (j in 1:8) {
  vvec <- c(crisG[j,j], everG[j,j], grahG[j,j], lineG[j,j], pulcG[j,j],
            sagrG[j,j], smarG[j,j])
  names(vvec) <- species
  Kvar[j] <- phylosig(anoletree, vvec)
}


# Evolutionary correlations (fig. 4, table A4)

phylovar <- c()
for (j in 1:8) {
  vvec <- c(crisG[j,j], everG[j,j], grahG[j,j], lineG[j,j], 
            sagrG[j,j], smarG[j,j])
  names(vvec) <- species6
  a.Y <- matrix(1, 1, 6) %*% invV %*% vvec/sum(invV)
  a.X <- matrix(1, 1, 6) %*% invV %*% eco/sum(invV)
  vvec.gls <- (vvec - c(a.Y)) %*% invV %*% 
    (eco - c(a.X))/sqrt(((vvec - c(a.Y)) %*% invV %*% (vvec - c(a.Y))) * 
                          ((eco - c(a.X)) %*% invV %*% (eco - c(a.X))))
  phylovar[j] <- as.numeric(vvec.gls)
}

phylr <- c()
for(j in 1:28){
  rvec <- c(as.numeric(rge[j,2]),as.numeric(rge[j,3]),as.numeric(rge[j,4]),
            as.numeric(rge[j,5]),as.numeric(rge[j,7]),as.numeric(rge[j,8]))
  names(rvec) <- species6
  a.Y <- matrix(1,1,6) %*% solve(V) %*% rvec/sum(solve(V));
  a.X <- matrix(1,1,6) %*% solve(V) %*% eco/sum(solve(V));
  rvec.gls<-(rvec-c(a.Y)) %*% solve(V) %*% (eco-c(a.X))/sqrt(((rvec-c(a.Y)) %*% 
      solve(V) %*% (rvec-c(a.Y)))*((eco-c(a.X)) %*% solve(V) %*% (eco-c(a.X))));
  phylr[j]<-as.numeric(rvec.gls)
}


