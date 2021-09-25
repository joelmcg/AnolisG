## Analyses using resampled G matrices

# Load required packages

library(evolqg)
library(mgcv)
library(geomorph)
library(ape)
library(phytools)
library(tidyverse)

# Species names
species <- c("cristatellus", "evermanni", "grahami", "lineatopus", "pulchellus", 
             "sagrei", "smaragdinus")
species6 <- c("cristatellus", "evermanni", "grahami", "lineatopus", "sagrei", 
              "smaragdinus")
spp <- c("cris", "ever", "grah", "line", "pulc", "sagr", "smar")

# Summary functions

summarize_cols <- function (x) {
  estimate <- apply(simplify2array(x), 2, median)
  ci.025 <- apply(simplify2array(x), 2, quantile, probs = 0.025)
  ci.975 <- apply(simplify2array(x), 2, quantile, probs = 0.975)
  data.frame(estimate, ci.025, ci.975)
}

summarize_reps <- function (x) {
  estimate <- apply(simplify2array(x), 1, median)
  ci.025 <- apply(simplify2array(x), 1, quantile, probs = 0.025)
  ci.975 <- apply(simplify2array(x), 1, quantile, probs = 0.975)
  data.frame(estimate, ci.025, ci.975)
}

summarize_vec <- function(x) {
  estimate <- median(x)
  ci.025 <- quantile(x, probs = 0.025)
  ci.975 <- quantile(x, probs = 0.975)
  df <- data.frame(estimate, ci.025, ci.975)
  row.names(df) <- NULL
  df
}


# Perform random skewers for all species on each set of replicates

rslist <- lapply(Gset, RandomSkewers, numvectors = 10000)
rslist <- transpose(rslist)$correlations


# Report mean and 95% confidence intervals across replicates for random skewers 
# correlations. Shows zeroes on diagonals and in upper triangle. (Table A2)

rslistmedian <- apply(simplify2array(rslist), 1:2, median)
rownames(rslistmedian) <- colnames(rslistmedian) <- spp

rslist025 <- apply(simplify2array(rslist), 1:2, quantile, probs = 0.025)
rownames(rslist025) <- colnames(rslist025) <- spp

rslist975 <- apply(simplify2array(rslist), 1:2, quantile, probs = 0.975)
rownames(rslist975) <- colnames(rslist975) <- spp


# Read trees and data for phylogenetic analyses. Trees as newick format, data as 
# csv with one line per species, listing element and SE (from ASReml) separately

rge <- as.data.frame(read.csv("data/rge.csv"))

a6tree <- read.tree("data/pruned3.tre")
anoletree <- read.tree("data/pruned7.tre")
anole.data <- as.data.frame(read.table("data/anoldat.txt"))
eco <- as.integer(anole.data$ecomorph)
names(eco) <- rownames(anole.data)

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

skewdist <- lapply(rslist, rsmatpcoa, species)

Kmult <- function(x, phy){
  physignal(x, phy, iter = 2)$phy.signal
}

Kmults <- sapply(skewdist, Kmult, anoletree)

summarize_vec(Kmults)

# Correlations between rRS and ecomorph (text)

eco.mat <- c(0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0)
eco6mat <- matrix(rep(0,36), 6, 6)
eco6mat[lower.tri(eco6mat, diag = FALSE)] <- 1 - eco.mat
eco6mat[upper.tri(eco6mat, diag = FALSE)] <- t(eco6mat)[upper.tri(eco6mat, 
                                                                  diag = FALSE)]
colnames(eco6mat) <- species6
row.names(eco6mat) <- species6

ecocors <- c()
for (i in 1:10000) {
  rs6 <- rslist[[i]][c(1:4, 6:7), c(1:4, 6:7)]
  skew1 <- sqrt(2*(1-as.dist(rs6)))
  ecocors[i] <- cor(skew1, as.dist(eco6mat))
}

summarize_vec(ecocors)


# Construct phylogenetic covariance matrix

V <- corMatrix(Initialize(corBrownian(phy = a6tree), anole.data))
invV <- solve(V)

# Resampled phylogenetic signal (table A3)

Kset <- list()
for (i in 1:10000) {
  Kset2 <- c()
  for (j in 1:28) {
    rvec <- c(rnorm(1, rge[j, 2], rge[j, 9]), rnorm(1, rge[j, 3], rge[j, 10]), 
              rnorm(1, rge[j, 4], rge[j, 11]), rnorm(1, rge[j, 5], rge[j, 12]), 
              rnorm(1, rge[j, 6], rge[j, 13]), rnorm(1, rge[j, 7], rge[j, 14]), 
              rnorm(1, rge[j, 8], rge[j, 15]))
    names(rvec) <- species
    Kset2[j] <- phylosig(anoletree, rvec)
  }
  Kset[[i]] <- Kset2
}
rm(rvec, Kset2)

Ksummary <- summarize_reps(Kset)

K.mat <- matrix(NA,8,8) 
K.mat[upper.tri(K.mat)] <- Ksummary[,1]
K.mat[lower.tri(K.mat)] <- t(K.mat)[lower.tri(t(K.mat))]


Kvarset <- list()
for (i in 1:10000) {
  Kset2 <- c ()
  for (j in 1:8) {
    vvec <- c(mycrisG[[i]][j,j], myeverG[[i]][j,j], mygrahG[[i]][j,j], 
              mylineG[[i]][j,j], mypulcG[[i]][j,j], mysagrG[[i]][j,j], 
              mysmarG[[i]][j,j])
    names(vvec) <- species
    Kset2[j] <- phylosig(anoletree, vvec)
  }
  Kvarset[[i]] <- Kset2
}
rm(vvec, Kset2)

Kvarsummary <- summarize_reps(Kvarset)



# Resample evolutionary correlation between ecomorph and genetic correlations
# (table A4)

rset <- c()
for (i in 1:10000) {
    rset2 <- c()
    for (j in 1:28) {
        rvec <- c(rnorm(1, rge[j, 2], rge[j, 9]), rnorm(1, rge[j, 3], 
                  rge[j, 10]), rnorm(1, rge[j, 4], rge[j, 11]), 
                  rnorm(1, rge[j, 5], rge[j, 12]), rnorm(1, rge[j, 7], 
                  rge[j, 14]), rnorm(1, rge[j, 8], rge[j, 15]))
        names(rvec) <- species6
        a.Y <- matrix(1, 1, 6) %*% invV %*% rvec/sum(invV)
        a.X <- matrix(1, 1, 6) %*% invV %*% eco/sum(invV)
        rvec.gls <- (rvec - c(a.Y)) %*% invV %*% 
          (eco - c(a.X))/sqrt(((rvec - c(a.Y)) %*% invV %*% (rvec - 
            c(a.Y))) * ((eco - c(a.X)) %*% invV %*% (eco - c(a.X))))
        rset2[j] <- as.numeric(rvec.gls)
    }
    rset[[i]] <- rset2
}

rm(a.Y, a.X, rvec, rvec.gls, rset2)


phylrsummary <- summarize_reps(rset)

phylr.mat <- matrix(NA,8,8) 
phylr.mat[upper.tri(phylr.mat)] <- phylrsummary[,1]
phylr.mat[lower.tri(phylr.mat)] <- t(phylr.mat)[lower.tri(t(phylr.mat))]

# Resample evolutionary correlation between ecomorph and genetic variances
# (table A4)

rvarset <- list()
for (i in 1:10000) {
  rset2 <- c()
  for (j in 1:8) {
    vvec <- c(mycrisG[[i]][j,j], myeverG[[i]][j,j], mygrahG[[i]][j,j], 
              mylineG[[i]][j,j], mysagrG[[i]][j,j], mysmarG[[i]][j,j])
    names(vvec) <- species6
    a.Y <- matrix(1, 1, 6) %*% invV %*% vvec/sum(invV)
    a.X <- matrix(1, 1, 6) %*% invV %*% eco/sum(invV)
    vvec.gls <- (vvec - c(a.Y)) %*% invV %*% (eco - 
                c(a.X))/sqrt(((vvec - c(a.Y)) %*% invV %*% (vvec - 
                c(a.Y))) * ((eco - c(a.X)) %*% invV %*% (eco - c(a.X))))
    rset2[j] <- as.numeric(vvec.gls)
  }
  rvarset[[i]] <- rset2
}
rm(vvec, vvec.gls, a.Y, a.X, rset2)

rvarsummary <- summarize_reps(rvarset)





