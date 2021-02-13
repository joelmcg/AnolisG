###########################################################################
#
# Resampling G matrices from McGlothlin et al. (2018) using the REML-MVN 
# method (Houle and Meyer 2015). REML-MVN is directly implemented in
# WOMBAT. This script uses output from ASReml to acheive the same results.
# G matrices are resampled 10000 times. Code can easily be modified to 
# change the number of replicates, the number of G matrices, etc.
#
###########################################################################


# Load packages -----------------------------------------------------------
library(mgcv)


# Generate pseudoreplicates of G matrices using REML-MVN ------------------


# Load G matrices, saved as a column from ASReml output in separate files

cris <- scan("data/crislist.txt")
ever <- scan("data/everlist.txt")
grah <- scan("data/grahlist.txt")
line <- scan("data/linelist.txt")
pulc <- scan("data/pulclist.txt")
sagr <- scan("data/sagrlist.txt")
smar <- scan("data/smarlist.txt")


# Import error structure from ASReml .vvp files 

# Each .vvp file has one triangle of the variance-covariance matrix of 
# parameters estimated by the model. Be sure to delete the first row in the
# raw file from ASReml-W before importing. Here, there are three submatrices 
# (E, G, and PE), along with covariances between submatrices. The last step
# extracts only error for G. Row/column numbers will vary based on size of
# analysis. Here there are 8 traits, so error for G starts on row 37. 


# Function for reading vvp files from ASReml

read_vvp <- function(file, traits, PE = TRUE) {
  
  # Use larger matrix if PE matrix is included 
  if (PE == FALSE) {
    size <- traits * (traits + 1)
  } else {
    size <- traits * (traits + 1) * 3 / 2
  }
  
  # Define matrix
  mat <- matrix(NA, size, size)
  
  # Read file and fit to matrix
  mat[upper.tri(mat, diag = TRUE)] <- scan(file)
  mat[lower.tri(mat, diag = FALSE)] <- t(mat)[lower.tri(t(mat), diag = FALSE)]
  
  # Return error matrix for G parameters only 
  x <- traits * (traits + 1) / 2 + 1
  y <- traits * (traits + 1)
  mat[x:y, x:y]
}


# Read .vvp files for each species

crisv <- read_vvp("data/cris.vvp", 8)
everv <- read_vvp("data/ever.vvp", 8)
grahv <- read_vvp("data/grah.vvp", 8)
linev <- read_vvp("data/line.vvp", 8)
pulcv <- read_vvp("data/pulc.vvp", 8)
sagrv <- read_vvp("data/sagr.vvp", 8)
smarv <- read_vvp("data/smar.vvp", 8)



# Function for REML-MVN sampling

remlMVN_G <- function (par, var, traits = 8, reps = 10000) {
  require(mgcv)
  
  # REML-MVN sampling of matrix
  G.samples <- rmvn(reps, par, var)
  
  # Reorganize samples into matrices
  G.samples.list <- list()

  for(i in 1:reps){
    G.i <- matrix(NA, traits, traits)
    G.i[upper.tri(G.i, diag = TRUE)] <- G.samples[i, ]
    G.i[lower.tri(G.i, diag = FALSE)] <- t(G.i)[lower.tri(t(G.i), diag = FALSE)]
    G.samples.list[[i]] <- G.i
  }  
  
  # Return samples as
  G.samples.list
} 


# Resample each G matrix

mycrisG <- remlMVN_G(cris, crisv)
myeverG <- remlMVN_G(ever, everv)
mygrahG <- remlMVN_G(grah, grahv)
mylineG <- remlMVN_G(line, linev)
mypulcG <- remlMVN_G(pulc, pulcv)
mysagrG <- remlMVN_G(sagr, sagrv)
mysmarG <- remlMVN_G(smar, smarv)


# List of resampled G matrices for use in analyses

Gset <- list()

for (i in 1:10000) {
  Gset[[i]] <- list(mycrisG[[i]], myeverG[[i]], mygrahG[[i]], mylineG[[i]], 
                    mypulcG[[i]], mysagrG[[i]], mysmarG[[i]])
}


# Save results for later

#saveRDS(mycrisG, file = "mycrisG.rds")
#saveRDS(myeverG, file = "myeverG.rds")
#saveRDS(mygrahG, file = "mygrahG.rds")
#saveRDS(mylineG, file = "mylineG.rds")
#saveRDS(mypulcG, file = "mypulcG.rds")
#saveRDS(mysagrG, file = "mysagrG.rds")
#saveRDS(mysmarG, file = "mysmarG.rds")
#saveRDS(Gset, file = "Gset.rds")

# Read in saved files (use these files to reproduce results in text)

#mycrisG <- readRDS("data/mycrisG.rds")
#myeverG <- readRDS("data/myeverG.rds")
#mygrahG <- readRDS("data/mygrahG.rds")
#mylineG <- readRDS("data/mylineG.rds")
#mypulcG <- readRDS("data/mypulcG.rds")
#mysagrG <- readRDS("data/mysagrG.rds")
#mysmarG <- readRDS("data/mysmarG.rds")
#Gset <- readRDS("data/Gset.rds")
