---
title: "Master Simulation Code"
output: html_document
date: "2025-05-07"
embed-resources: true
format: 
  html:
    code-fold: true
    toc: true
editor: source
editor_options: 
  chunk_output_type: console
execute: 
  echo: true
  warning: false
  message: false
  fig-height: 8
  fig-width: 8
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, eval=FALSE)
```

## **Simulation Code Chunks:**
#### Set up libraries needed throughout:
```{r}
#install.packages(c("MASS","mvSLOUCH","ape","matrixNormal","data.table","XML","methods","dplyr","effectsize","HelpersMG","tidyverse"))

library(MASS)
library(mvSLOUCH)
library(ape)
library(phytools)
library(matrixNormal)
library(data.table)
library(XML)
library(methods)
library(dplyr)
library(effectsize)
library(HelpersMG)
library(tidyverse)
```


### **"Sweet Spot" Scenario**
#### Set up the initial dataframe that will guide the simulations in Unity - use the same dataframe for all to keep things consistent
```{r}
# Step 1: Set parameters for the analysis
### Make a data matrix of all possible combinations of 1-50 taxa and 1-50 traits
### Each step of the analysis will call on a unique combination to fill in the data matrix with that combination of ntax and ntraits
mydat <- as.data.frame(
  expand.grid(
    ntraits = 2:50,
    ntax = 2:50)) %>%
  .[sample(nrow(.)),] %>%
  group_by(ntax) %>%
  mutate(run = 2:50) %>%
  ungroup()
### Add in columns for the number of traits we expect to be pulled out
mydat$iteration<-1:length(mydat$ntraits)
### Parameters for simulating the loadings matrix
pfac1 = 0.50 # percent of traits factor 1 influences
nfac1 = floor(pfac1*mydat$ntraits)
pfac2 = 0.50 # percent of trait factor 2 influences
nfac2 = floor(pfac2*mydat$ntraits)
### Add in more columns for the number of traits we expect to be pulled out
mydat$fac.expec<-mydat$ntraits*pfac1
mydat$iteration<-1:length(mydat$ntraits)
mydat$fac.expec<-mydat$ntraits*pfac1
mydat$fac.expec.up<-ceiling(mydat$fac.expec)
mydat$fac.expec.down<-floor(mydat$fac.expec)
# Save the dataframe for use in the simulations
write.csv(mydat, file="mydat.csv", row.names=F)
# Steps 2-7: this will happen in the following code chunks
### Step 2: simulate phylogeny
### Step 3: generate factor values at base of phylogeny
### Step 4: simulate factor values at tips under BM
### Step 5: generate loadings matrix (L = K x P = factors x traits)
### Step 6: simulate error and add into general model
### Step 7: write final trait data to file
```

#### Generalized code block for all simulations
```{r}
### Set up parameters for simulating the loadings matrix
ntax = 10 # the number of tips in the tree (btw 1-50)
ntraits = 10 # the number of traits (btw 1-50)
nfactor = 2 # the number of factors
### Parameters for simulating the loadings matrix
pfac1 = 0.50 # percent of traits factor 1 influences
nfac1 = floor(pfac1*ntraits)
pfac2 = 0.50 # percent of trait factor 2 influences
nfac2 = floor(pfac2*ntraits)

### Either load in phylogeny or simulate a new one
# Load in phylogeny
phyltree<-read.tree("path to directory/simtree.txt")
# Simulate new phylogeny 
phyltree<-ape::rtree(ntax)
phyltree$tip.label<-gsub("t","",as.character(phyltree$tip.label)) # drops the "t" from the tip names to allow for trimming
phyltree$tip.label<-paste0("taxon",as.character(phyltree$tip.label)) # add the taxon labels back in for easier ID
phyltree<-phyltree_paths(phyltree) # advisable for speed
write.tree(phyltree, file="simtree.txt")

### Simulate root node
kappa0 = 1.0
factor.varcov = diag(nfactor)  # this creates the identity matrix in R
factor.root = mvrnorm(mu=rep(0,nfactor),Sigma=(1/kappa0)*factor.varcov)

### Simulate factor values at tips
myfactors = simulBMProcPhylTree(phyltree,X0 = factor.root,Sigma=factor.varcov)

### Specify loadings for n factors
myloadings = matrix(c(rep(1.0,nfac1),rep(0,ntraits-nfac1),rep(0,ntraits-nfac2),rep(1.0,nfac2)),nrow=nfactor,byrow=T)

### Simulate error
mean.mat = matrix(0,ncol=ntraits,nrow=ntax) # mean for error matrix
trait.errors = rep(1.0,ntraits) # variance of each trait
### Draw error from matrix multivariate normal distributions
myU = diag(ntax)
myV = solve(diag(trait.errors))
myerror = rmatnorm(s=1,M=mean.mat,U=myU,V=myV)

### Write down final data
simdat0 = as.data.table(myfactors %*% myloadings + myerror, keep.rownames="taxon") #gives the taxon column the proper name so Julia can recognize it
### Standardize final data - subtract mean & divide by SD on trait-by-trait basis using standardize function from effectsize package
simdat = standardize(simdat0,robust=F) 
### Write the final continuous trait dataset to csv
write.csv(simdat,"simdat.csv", row.names=F) #this will delete the first column of numbers so we don't have to edit the csv before input into julia
```

#### This code uses the same simulated tree (with 50 tips) to test the effect of keeping the same phylogenetic structure - simtree50 & simtree502
```{r}
#Create a new 50-tip simulated tree and write it to file - this is the tree we will use for all these simulations
phyltree<-ape::rtree(50)
write.tree(phyltree,file="OGsimtree50.txt")

#Now do the simulations
#Load in mydat matrix to guide simulations
mydat<-read.csv("mydat.csv")
#Run through for-loop to generate trait datasets
for (i in 1:nrow(mydat)){
  ### Set up parameters for simulating the loadings matrix
  ntax = mydat$ntax[i]
  ntraits = mydat$ntraits[i]
  nfactor = 2
  ### Parameters for simulating the loadings matrix
  pfac1 = 0.50
  nfac1 = floor(pfac1*ntraits)
  pfac2 = 0.50
  nfac2 = floor(pfac2*ntraits)
  ### Load in phylogeny and trim to number of taxa included in the analysis
  phyltree<-read.tree("OGsimtree50.txt")
  phyltree$tip.label<-gsub("t","",as.character(phyltree$tip.label))
  phyltree<-drop.tip(phyltree, sample(phyltree$tip.label)[0:(50-ntax)])
  phyltree$tip.label<-paste0("taxon",as.character(phyltree$tip.label))
  write.tree(phyltree, file=paste("simtree",mydat$iteration[i],".txt", sep=""))
  phyltree<-phyltree_paths(phyltree)
  ### Simulate root node
  kappa0 = 1.0
  factor.varcov = diag(nfactor)
  factor.root = mvrnorm(mu=rep(0,nfactor),Sigma=(1/kappa0)*factor.varcov)
  ### Simulate factor values at tips
  myfactors = simulBMProcPhylTree(phyltree,X0 = factor.root,Sigma=factor.varcov)
  ### Simulate loadings for n factors
  myloadings = matrix(c(rep(1.0,nfac1),rep(0,ntraits-nfac1),rep(0,ntraits-nfac2),rep(1.0,nfac2)),nrow=nfactor,byrow=T)
  ### Simulate error
  mean.mat = matrix(0,ncol=ntraits,nrow=ntax)
  trait.errors = rep(1.0,ntraits)
  ### Draw error from matrix multivariate normal distribution
  myU = diag(ntax)
  myV = solve(diag(trait.errors))
  myerror = rmatnorm(s=1,M=mean.mat,U=myU,V=myV)
  ### Write down final data
  simdat0 = as.data.table(myfactors %*% myloadings + myerror, keep.rownames="taxon")
  ### Standardize the trait data - mean=0 and SD=1
  simdat = standardize(simdat0,robust=F) 
  ### Save final data
  write.csv(simdat,paste("simdat",mydat$iteration[i],".csv", sep=""), row.names=F)
}
```

#### This code simulates a new tree each time to test the effect of different phylogenetic structures on the analysis - newsimtree & newsimtree2
```{r}
#Load in mydat matrix to guide simulations
mydat<-read.csv("mydat.csv")
#Simulate trait datasets
for (i in 1:nrow(mydat)){
  ### Set up parameters for simulating the loadings matrix
  ntax = mydat$ntax[i]
  ntraits = mydat$ntraits[i]
  nfactor = 2
  ### Parameters for simulating the loadings matrix
  pfac1 = 0.50
  nfac1 = floor(pfac1*ntraits)
  pfac2 = 0.50
  nfac2 = floor(pfac2*ntraits)
  ### Generate a random phylogeny
  phyltree<-ape::rtree(ntax)
  ### Write Newick tree to file so we have record of which tree goes with which iteration
  phyltree$tip.label<-gsub("t","",as.character(phyltree$tip.label))
  phyltree$tip.label<-paste0("taxon",as.character(phyltree$tip.label))
  write.tree(phyltree, file=paste("simtree",mydat$iteration[i],".txt", sep=""))
  phyltree<-phyltree_paths(phyltree)
  ### Simulate root node
  kappa0 = 1.0
  factor.varcov = diag(nfactor)
  factor.root = mvrnorm(mu=rep(0,nfactor),Sigma=(1/kappa0)*factor.varcov)
  ### Simulate factor values at tips
  myfactors = simulBMProcPhylTree(phyltree,X0 = factor.root,Sigma=factor.varcov)
  ### Simulate loadings for n factors
  myloadings = matrix(c(rep(1.0,nfac1),rep(0,ntraits-nfac1),rep(0,ntraits-nfac2),rep(1.0,nfac2)),nrow=nfactor,byrow=T)
  ### Simulate error
  mean.mat = matrix(0,ncol=ntraits,nrow=ntax)
  trait.errors = rep(1.0,ntraits)
  ### Draw error from matrix multivariate normal distribution
  myU = diag(ntax)
  myV = solve(diag(trait.errors))
  myerror = rmatnorm(s=1,M=mean.mat,U=myU,V=myV)
  ### Write down final data
  simdat0 = as.data.table(myfactors %*% myloadings + myerror, keep.rownames="taxon")
  ### Standardize the trait data - mean=0 and SD=1
  simdat = standardize(simdat0,robust=F)
  ### Write final data to file
  write.csv(simdat,paste("simdat",mydat$iteration[i],".csv", sep=""), row.names=F)
}
```


### **Error Scenario**
#### Simulate a base tree to run different errors with
```{r}
### Set up parameters for simulating the loadings matrix
ntax = 20 # the number of tips in the tree (btw 1-50)
ntraits = 20 # the number of traits (btw 1-50)
nfactor = 2 # the number of factors (stick with 2 for now)
### Parameters for simulating the loadings matrix
pfac1 = 0.50 # percent of traits factor 1 influences
nfac1 = floor(pfac1*ntraits)
pfac2 = 0.50 # percent of trait factor 2 influences
nfac2 = floor(pfac2*ntraits)
### code for simulating and saving new trees each time - add this into final for loop for trying other types of simulation
# generate a random phylogeny
phyltree<-ape::rtree(ntax)
# write Newick tree to file so we have record of which tree goes with which iteration
phyltree$tip.label<-gsub("t","",as.character(phyltree$tip.label)) # drops the "t" from the tip names to allow for trimming
phyltree$tip.label<-paste0("taxon",as.character(phyltree$tip.label)) # add the taxon labels back in for easier ID 
write.tree(phyltree, file=paste("errorsimtree.txt", sep=""))
```

#### Give half the traits a bigger variance to test how it affects the factors recognizing them - errorsimdat1
```{r}
phyltree<-read.tree("errorsimtree.txt")
phyltree<-phyltree_paths(phyltree) # advisable for speed
### Simulate root node
kappa0 = 1.0
factor.varcov = diag(nfactor)  # this creates the identity matrix in R
factor.root = mvrnorm(mu=rep(0,nfactor),Sigma=(1/kappa0)*factor.varcov)
### Simulate factor values at tips
myfactors = simulBMProcPhylTree(phyltree,X0 = factor.root,Sigma=factor.varcov)
### Simulate loadings for n factors
myloadings = matrix(c(rep(1.0,nfac1),rep(0,ntraits-nfac1),rep(0,ntraits-nfac2),rep(1.0,nfac2)),nrow=nfactor,byrow=T)
### Simulate error
mean.mat = matrix(0,ncol=ntraits,nrow=ntax) # mean for error matrix
trait.errors = c(rep(1, ntraits/2),rep(5, ntraits/2)) # variance of each trait - change this to modify the error the same
### Draw error from matrix multivariate normal distribution
myU = diag(ntax)
myV = solve(diag(trait.errors))
myerror = rmatnorm(s=1,M=mean.mat,U=myU,V=myV)
### Write down final data
simdat0 = as.data.table(myfactors %*% myloadings + myerror, keep.rownames="taxon")
### Standardize the trait data - mean=0 and SD=1
simdat = standardize(simdat0,robust=F)

write.csv(simdat, paste("errorsimdat1.csv", sep=""), row.names=F) #this will delete the first column of numbers so we don't have to edit the csv
```

#### Give the last couple traits of each half a huge variance - errorsimdat2
```{r}
phyltree<-read.tree("errorsimtree.txt")
phyltree<-phyltree_paths(phyltree) # advisable for speed
### Simulate root node
kappa0 = 1.0
factor.varcov = diag(nfactor)  # this creates the identity matrix in R
factor.root = mvrnorm(mu=rep(0,nfactor),Sigma=(1/kappa0)*factor.varcov)
### Simulate factor values at tips
myfactors = simulBMProcPhylTree(phyltree,X0 = factor.root,Sigma=factor.varcov)
### Simulate loadings for n factors
myloadings = matrix(c(rep(1.0,nfac1),rep(0,ntraits-nfac1),rep(0,ntraits-nfac2),rep(1.0,nfac2)),nrow=nfactor,byrow=T)
### Simulate error
mean.mat = matrix(0,ncol=ntraits,nrow=ntax) # mean for error matrix
trait.errors = c(rep(1, (ntraits/2)-2),rep(50, 2),rep(1, (ntraits/2)-2),rep(50,2)) # variance of each trait - change this to modify the error the same
### Draw error from matrix multivariate normal distribution
myU = diag(ntax)
myV = solve(diag(trait.errors))
myerror = rmatnorm(s=1,M=mean.mat,U=myU,V=myV)
### Write down final data
simdat0 = as.data.table(myfactors %*% myloadings + myerror, keep.rownames="taxon")
### Standardize the trait data - mean=0 and SD=1
simdat = standardize(simdat0,robust=F)

write.csv(simdat, paste("errorsimdat2.csv", sep=""), row.names=F) #this will delete the first column of numbers so we don't have to edit the csv
```

#### Give the last couple traits of each half a bigger variance - errorsimdat3
```{r}
phyltree<-read.tree("errorsimtree.txt")
phyltree<-phyltree_paths(phyltree) # advisable for speed
### Simulate root node
kappa0 = 1.0
factor.varcov = diag(nfactor)  # this creates the identity matrix in R
factor.root = mvrnorm(mu=rep(0,nfactor),Sigma=(1/kappa0)*factor.varcov)
### Simulate factor values at tips
myfactors = simulBMProcPhylTree(phyltree,X0 = factor.root,Sigma=factor.varcov)
### Simulate loadings for n factors
myloadings = matrix(c(rep(1.0,nfac1),rep(0,ntraits-nfac1),rep(0,ntraits-nfac2),rep(1.0,nfac2)),nrow=nfactor,byrow=T)
### Simulate error
mean.mat = matrix(0,ncol=ntraits,nrow=ntax) # mean for error matrix
trait.errors = c(rep(1, (ntraits/2)-2),rep(100, 2),rep(1, (ntraits/2)-2),rep(100,2)) # variance of each trait - change this to modify the error the same
### Draw error from matrix multivariate normal distribution
myU = diag(ntax)
myV = solve(diag(trait.errors))
myerror = rmatnorm(s=1,M=mean.mat,U=myU,V=myV)
### Write down final data
simdat0 = as.data.table(myfactors %*% myloadings + myerror, keep.rownames="taxon")
### Standardize the trait data - mean=0 and SD=1
simdat = standardize(simdat0,robust=F)

write.csv(simdat, paste("errorsimdat3.csv", sep=""), row.names=F) #this will delete the first column of numbers so we don't have to edit the csv
```

#### Replace the last couple traits of each half with random numbers - errorsimdat4
```{r}
phyltree<-read.tree("errorsimtree.txt")
phyltree<-phyltree_paths(phyltree) # advisable for speed
### Simulate root node
kappa0 = 1.0
factor.varcov = diag(nfactor)  # this creates the identity matrix in R
factor.root = mvrnorm(mu=rep(0,nfactor),Sigma=(1/kappa0)*factor.varcov)
### Simulate factor values at tips
myfactors = simulBMProcPhylTree(phyltree,X0 = factor.root,Sigma=factor.varcov)
### Simulate loadings for n factors
myloadings = matrix(c(rep(1.0,nfac1),rep(0,ntraits-nfac1),rep(0,ntraits-nfac2),rep(1.0,nfac2)),nrow=nfactor,byrow=T)
### Simulate error
mean.mat = matrix(0,ncol=ntraits,nrow=ntax) # mean for error matrix
trait.errors = c(rep(1.0, ntraits)) # variance of each trait - change this to modify the error the same
### Draw error from matrix multivariate normal distribution
myU = diag(ntax)
myV = solve(diag(trait.errors))
myerror = rmatnorm(s=1,M=mean.mat,U=myU,V=myV)
### Write down final data
simdat0 = as.data.table(myfactors %*% myloadings + myerror, keep.rownames="taxon")
### Standardize the trait data - mean=0 and SD=1
simdat = standardize(simdat0,robust=F)
#replace the last two trait columns for each factor with 
simdat$V9 = sample(-100:100, ntraits, replace=T)
simdat$V10 = sample(-100:100, ntraits, replace=T)
simdat$V19 = sample(-100:100, ntraits, replace=T)
simdat$V20 = sample(-100:100, ntraits, replace=T)

write.csv(simdat, paste("errorsimdat4.csv", sep=""), row.names=F) #this will delete the first column of numbers so we don't have to edit the csv
```

#### Add more tips and traits and make several traits with huge variance - error simdat5
#### Simulate new phylogeny
```{r}
### Set up parameters for simulating the loadings matrix
ntax = 40 # the number of tips in the tree (btw 1-50)
ntraits = 40 # the number of traits (btw 1-50)
nfactor = 2 # the number of factors (stick with 2 for now)
### Parameters for simulating the loadings matrix
pfac1 = 0.50 # percent of traits factor 1 influences
nfac1 = floor(pfac1*ntraits)
pfac2 = 0.50 # percent of trait factor 2 influences
nfac2 = floor(pfac2*ntraits)
# generate a random phylogeny
phyltree<-ape::rtree(ntax)
# write Newick tree to file so we have record of which tree goes with which iteration
phyltree$tip.label<-gsub("t","",as.character(phyltree$tip.label)) # drops the "t" from the tip names to allow for trimming
phyltree$tip.label<-paste0("taxon",as.character(phyltree$tip.label)) # add the taxon labels back in for easier ID 
write.tree(phyltree, file=paste("errorsimtree40.txt", sep=""))
```

#### Simulate trait data
```{r}
#load in simulated tree from working dir
phyltree<-read.tree("errorsimtree40.txt")
phyltree<-phyltree_paths(phyltree) # advisable for speed
### Simulate root node
kappa0 = 1.0
factor.varcov = diag(nfactor)  # this creates the identity matrix in R
factor.root = mvrnorm(mu=rep(0,nfactor),Sigma=(1/kappa0)*factor.varcov)
### Simulate factor values at tips
myfactors = simulBMProcPhylTree(phyltree,X0 = factor.root,Sigma=factor.varcov)
### Simulate loadings for n factors
myloadings = matrix(c(rep(1.0,nfac1),rep(0,ntraits-nfac1),rep(0,ntraits-nfac2),rep(1.0,nfac2)),nrow=nfactor,byrow=T)
### Simulate error
mean.mat = matrix(0,ncol=ntraits,nrow=ntax) # mean for error matrix
trait.errors = c(rep(1.0, ntraits)) # variance of each trait - change this to modify the error the same
### Draw error from matrix multivariate normal distribution
myU = diag(ntax)
myV = solve(diag(trait.errors))
myerror = rmatnorm(s=1,M=mean.mat,U=myU,V=myV)
### Write down final data
simdat0 = as.data.table(myfactors %*% myloadings + myerror, keep.rownames="taxon")
### Standardize the trait data - mean=0 and SD=1
simdat = standardize(simdat0,robust=F)
### Replace a couple trait columns for each factor with random numbers to create a large variance
simdat$V8 = sample(-100:100, ntraits, replace=T)
simdat$V13 = sample(-100:100, ntraits, replace=T)
simdat$V18 = sample(-100:100, ntraits, replace=T)

simdat$V28 = sample(-100:100, ntraits, replace=T)
simdat$V33 = sample(-100:100, ntraits, replace=T)
simdat$V38 = sample(-100:100, ntraits, replace=T)

write.csv(simdat, paste("errorsimdat40.csv", sep=""), row.names=F) #this will delete the first column of numbers so we don't have to edit the csv
```


### **Discrete Trait Scenario**
#### Simulate a base tree - 20 tips
```{r}
### Set up parameters for simulating the loadings matrix
ntax = 20 # the number of tips in the tree (btw 1-50)
ntraits = 20 # the number of traits (btw 1-50)
nfactor = 2 # the number of factors (stick with 2 for now)
### Parameters for simulating the loadings matrix
pfac1 = 0.50 # percent of traits factor 1 influences
nfac1 = floor(pfac1*ntraits)
pfac2 = 0.50 # percent of trait factor 2 influences
nfac2 = floor(pfac2*ntraits)
### Generate a random 20 tip phylogeny
phyltree<-ape::rtree(ntax)
### Write Newick tree to file so we have record of which tree goes with which iteration
phyltree$tip.label<-gsub("t","",as.character(phyltree$tip.label)) # drops the "t" from the tip names to allow for trimming
phyltree$tip.label<-paste0("taxon",as.character(phyltree$tip.label)) # add the taxon labels back in for easier ID 
write.tree(phyltree, file=paste("simtree.txt", sep=""))
```

#### Simulate a base trait matrix - 20 traits
```{r}
phyltree<-read.tree("simtree.txt")
phyltree<-phyltree_paths(phyltree) # advisable for speed
### Simulate root node
kappa0 = 1.0
factor.varcov = diag(nfactor)  # this creates the identity matrix in R
factor.root = mvrnorm(mu=rep(0,nfactor),Sigma=(1/kappa0)*factor.varcov)
### Simulate factor values at tips
myfactors = simulBMProcPhylTree(phyltree,X0 = factor.root,Sigma=factor.varcov)
### Simulate loadings for n factors
myloadings = matrix(c(rep(1.0,nfac1),rep(0,ntraits-nfac1),rep(0,ntraits-nfac2),rep(1.0,nfac2)),nrow=nfactor,byrow=T)
### Simulate error
mean.mat = matrix(0,ncol=ntraits,nrow=ntax) # mean for error matrix
trait.errors = c(rep(1.0, ntraits)) # variance of each trait
### Draw error from matrix multivariate normal distribution
myU = diag(ntax)
myV = solve(diag(trait.errors))
myerror = rmatnorm(s=1,M=mean.mat,U=myU,V=myV)
### Write down final data
simdat0 = as.data.table(myfactors %*% myloadings + myerror, keep.rownames="taxon")
### Standardize the trait data - mean=0 and SD=1
simdat = standardize(simdat0,robust=F)
### Save trait data to file
write.csv(simdat, paste("simdat0.csv", sep=""), row.names=F)
```

#### Edit the base trait dataset to convert different parts to binary
```{r}
### Set up the function to apply to the trait datasets for converting columns to binary
func<-function(x){
  factor(ifelse(x>0,"1","0"))
}
### Read in base trait dataset
traits<-read.csv("simdat0.csv")
### Edit the trait datasets and save results to file
#### 2 traits per factor
traits1<-traits %>% mutate(across(c(V2,V9,V12,V19), func))
write.csv(traits1, paste("simdat1.csv", sep=""), row.names=F)
#### 5 traits per factor
traits2<-traits %>% mutate(across(c(V2,V3,V4,V8,V9,V12,V13,V14,V18,V19),func))
write.csv(traits2, paste("simdat2.csv", sep=""), row.names=F)
#### 8 traits per factor
traits3<-traits %>% mutate(across(c(V2,V3,V4,V5,V6,V7,V8,V9,V12,V13,V14,V15,V16,V17,V18,V19),func))
write.csv(traits3, paste("simdat3.csv", sep=""), row.names=F)
#### All traits in the dataset
traits4<-traits %>% mutate(across(c(V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12,V13,V14,V15,V16,V17,V18,V19,V20),func))
write.csv(traits4, paste("simdat.csv", sep=""), row.names=F)
#### 5 traits only on factor 1
traits5<-traits %>% mutate(across(c(V2,V3,V4,V8,V9),func))
write.csv(traits5, paste("simdat5.csv", sep=""), row.names=F)
#### 5 traits only on factor 2
traits6<-traits %>% mutate(across(c(V12,V13,V14,V18,V19),func))
write.csv(traits6, paste("simdat6.csv", sep=""), row.names=F)
#### All traits on factor 1
traits7<-traits %>% mutate(across(c(V1,V2,V3,V4,V5,V6,V7,V8,V9,V10),func))
write.csv(traits7, paste("simdat7.csv", sep=""), row.names=F)
#### All traits on factor 2
traits8<-traits %>% mutate(across(c(V11,V12,V13,V14,V15,V16,V17,V18,V19,V20),func))
write.csv(traits8, paste("simdat8.csv", sep=""), row.names=F)
```

#### Now rerun it all with 40 tips, 40 traits
#### Simulate a base tree - 40 tips
```{r}
### Set up parameters for simulating the loadings matrix
ntax = 40 # the number of tips in the tree (btw 1-50)
ntraits = 40 # the number of traits (btw 1-50)
nfactor = 2 # the number of factors (stick with 2 for now)
### Parameters for simulating the loadings matrix
pfac1 = 0.50 # percent of traits factor 1 influences
nfac1 = floor(pfac1*ntraits)
pfac2 = 0.50 # percent of trait factor 2 influences
nfac2 = floor(pfac2*ntraits)
### Generate random 40 tip phylogeny
phyltree<-ape::rtree(ntax)
### write Newick tree to file so we have record of which tree goes with which iteration
phyltree$tip.label<-gsub("t","",as.character(phyltree$tip.label)) # drops the "t" from the tip names to allow for trimming
phyltree$tip.label<-paste0("taxon",as.character(phyltree$tip.label)) # add the taxon labels back in for easier ID 
write.tree(phyltree, file=paste("simtree.txt", sep=""))
```

#### Simulate a base trait matrix - 40 traits
```{r}
phyltree<-read.tree("simtree.txt")
phyltree<-phyltree_paths(phyltree) # advisable for speed
### Simulate root node
kappa0 = 1.0
factor.varcov = diag(nfactor)  # this creates the identity matrix in R
factor.root = mvrnorm(mu=rep(0,nfactor),Sigma=(1/kappa0)*factor.varcov)
### Simulate factor values at tips
myfactors = simulBMProcPhylTree(phyltree,X0 = factor.root,Sigma=factor.varcov)
### Simulate loadings for n factors
myloadings = matrix(c(rep(1.0,nfac1),rep(0,ntraits-nfac1),rep(0,ntraits-nfac2),rep(1.0,nfac2)),nrow=nfactor,byrow=T)
### Simulate error
mean.mat = matrix(0,ncol=ntraits,nrow=ntax) # mean for error matrix
trait.errors = c(rep(1.0, ntraits)) # variance of each trait
### Draw error from matrix multivariate normal distribution
myU = diag(ntax)
myV = solve(diag(trait.errors))
myerror = rmatnorm(s=1,M=mean.mat,U=myU,V=myV)
### Write down final data
simdat0 = as.data.table(myfactors %*% myloadings + myerror, keep.rownames="taxon")
### Standardize the trait data - mean=0 and SD=1
simdat = standardize(simdat0,robust=F)
write.csv(simdat, paste("simdat0.csv", sep=""), row.names=F)
```

#### Edit the base trait dataset to convert different parts to binary
```{r}
### Set up the function to apply to the trait datasets for converting columns to binary
func<-function(x){
  factor(ifelse(x>0,"1","0"))
}
### Read in base trait dataset
traits<-read.csv("simdat0.csv")
### Edit the trait datasets and save results to file
#### 2 traits per factor
traits1<-traits %>% mutate(across(c(V5,V15,V25,V35), func))
write.csv(traits1, paste("simdat1.csv", sep=""), row.names=F)
#### 5 traits per factor
traits2<-traits %>% mutate(across(c(V5,V8,V10,V15,V18,V25,V28,V30,V35,V38),func))
write.csv(traits2, paste("simdat2.csv", sep=""), row.names=F)
#### 10 traits per factor
traits3<-traits %>% mutate(across(c(V5,V6,V7,V8,V10,V11,V15,V16,V17,V18,V25,V26,V27,V28,V30,V31,V35,V36,V37,V38),func))
write.csv(traits3, paste("simdat3.csv", sep=""), row.names=F)
#### All traits in the dataset
traits4<-traits %>% mutate(across(c(V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12,V13,V14,V15,V16,V17,V18,V19,V20,V21,V22,V23,V24,V25,V26,V27,V28,V29,V30,V31,V32,V33,V34,V35,V36,V37,V38,V39,V40),func))
write.csv(traits4, paste("simdat.csv", sep=""), row.names=F)
#### 10 traits only on factor 1
traits5<-traits %>% mutate(across(c(V5,V6,V7,V8,V10,V11,V15,V16,V17,V18),func))
write.csv(traits5, paste("simdat5.csv", sep=""), row.names=F)
#### 10 traits only on factor 2
traits6<-traits %>% mutate(across(c(V25,V26,V27,V28,V30,V31,V35,V36,V37,V38),func))
write.csv(traits6, paste("simdat6.csv", sep=""), row.names=F)
#### All traits on factor 1
traits7<-traits %>% mutate(across(c(V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12,V13,V14,V15,V16,V17,V18,V19,V20),func))
write.csv(traits7, paste("simdat7.csv", sep=""), row.names=F)
#### All traits on factor 2
traits8<-traits %>% mutate(across(c(V21,V22,V23,V24,V25,V26,V27,V28,V29,V30,V31,V32,V33,V34,V35,V36,V37,V38,V39,V40),func))
write.csv(traits8, paste("simdat8.csv", sep=""), row.names=F)
```


### **New "Build-a-Bat" Scenario**
##### This is a 3-factor case with 20 traits and 40 taxa

All Traits: Brain(g), BL(mm), BM(g), CBL, WL(Nm2), AR, WS(mm), FA(mm), FMature, MMature, Litter, Life, Forage, HB(km2), Lon(deg), Elev(km), Realm, Echo, Diet, Hib

Traits converted to binary: Litter, Forage, Realm, Echo, Diet, Hib

##### Traits 1-8: Physiological/Morphological
###### Trait 1: brain mass = Brain(g), quantitative
###### Trait 2: body length = BL(mm), quantitative
###### Trait 3: body mass = BM(g), quantitative
###### Trait 4: condylobasal length = (CBL), quantitative
###### Trait 5: wing loading = WL(Nm2), quantitative
###### Trait 6: wing aspect ratio = AR, quantitiative
###### Trait 7: wingspan = WS(mm), quantitative
###### Trait 8: forearm length = FA(mm), quantitative

##### Traits 9-12: Reproductive
###### Trait 9: female maturity age = (FMature), quantitative
###### Trait 10: male maturity age = (MMature), quantitative
###### Trait 11: litter size = (Litter), qualitative (1=more than 2, 0=less than 2)
###### Trait 12: max longevity = (Life), quantitative

##### Traits 13-18: Ecological
###### Trait 13: foraging type = (Forage), qualitative (1=scansorial, 0=arboreal)
###### Trait 14: habitat breadth = (HB(km2)), quantitative
###### Trait 15: longitudinal extent = (Lon(deg)), quantitative
###### Trait 16: elevation range = (Elev(km)), quantitative
###### Trait 17: biological realm = (Realm), qualitative (1=nearctic, 0=neotropical)
###### Trait 18: echolocation type = (Echo), qualitative (1=oral, 0=nasal)
###### Trait 19: diet type = (Diet), qualitative (1=insectivore, 0=frugivore)
###### Trait 20: hibernation = (Hib), qualitative (1=hibernates, 0=not hibernates)
```{r}
setwd("C:/Users/shelb/OneDrive/Desktop/Chapter 2 Stuff/Simulation/script files/buildabat") 
### Set up parameters for simulating the loadings matrix
ntax = 40 # the number of tips in the tree
ntraits = 20 # the number of traits
nfactor = 3 # the number of factors
### Parameters for simulating the loadings matrix
pfac1 = 0.40 # percent of traits factor 1 influences - should be traits 1-8
nfac1 = floor(pfac1*ntraits)
pfac2 = 0.20 # percent of trait factor 2 influences - should be traits 9-12
nfac2 = floor(pfac2*ntraits)
pfac3 = 0.40 # percent of trait factor 2 influences - should be traits 13-20
nfac3 = floor(pfac3*ntraits)
### Create phylogeny and save it so we can use it again
phyltree<-ape::rtree(ntax)
phyltree$tip.label<-gsub("t","",as.character(phyltree$tip.label)) # drops the "t" from the tip names to allow for trimming
phyltree$tip.label<-paste0("taxon",as.character(phyltree$tip.label)) # add the taxon labels back in for easier ID
phyltree<-phyltree_paths(phyltree) # advisable for speed
write.tree(phyltree, file="simtree.txt")
### Simulate root node
kappa0 = 1.0
factor.varcov = diag(nfactor)  # this creates the identity matrix in R
factor.root = mvrnorm(mu=rep(0,nfactor),Sigma=(1/kappa0)*factor.varcov)
### Simulate factor values at tips
myfactors = simulBMProcPhylTree(phyltree,X0 = factor.root,Sigma=factor.varcov)
### Specify loadings for n factors
myloadings = matrix(c(1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0, 
                      0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1),nrow=nfactor,byrow=T)
### Simulate error
mean.mat = matrix(0,ncol=ntraits,nrow=ntax) # mean for error matrix
trait.errors = rep(1.0,ntraits) # variance of each trait
### Set up variance-covariance matrix to include correlation for certain traits
myU = diag(ntax)
myV = solve(diag(trait.errors))
### Traits 1-8 are F1, traits 9-12 are F2, traits 13-20 are F3
### Set up correlation values to add into the V matrix
corr1=0.5
corr2=0.3
corr3=0.9
### Set up the new correlations matrix
### Make traits 1-8 correlated
myV[1:8,1:8]<-corr3
### Make traits 9-12 correlated
myV[9:12,9:12]<-corr2
### Make traits 13-20 correlated
myV[13:20,13:20]<-corr1

### Edit the V matrix to ensure it functions properly in the error simulation
diag(myV)<-1
all(eigen(myV)$values > 0) #this tests to ensure that the matrix works properly - if returns TRUE you're good to go
colnames(myV)<-c("Brain(g)", "BL(mm)", "BM(g)", "CBL", "WL(Nm2)", "AR", "WS(mm)", "FA(mm)", "FMature", "MMature", "Litter", "Life", "Forage", "HB(km2)", "Lon(deg)", "Elev(km)", "Realm", "Echo", "Diet", "Hib")
row.names(myV)<-c("Brain(g)", "BL(mm)", "BM(g)", "CBL", "WL(Nm2)", "AR", "WS(mm)", "FA(mm)", "FMature", "MMature", "Litter", "Life", "Forage", "HB(km2)", "Lon(deg)", "Elev(km)", "Realm", "Echo", "Diet", "Hib")

### Draw error from matrix multivariate normal distribution
myerror = rmatnorm(s=1,M=mean.mat,U=myU,V=myV)
### Write down final trait data
simdat0 = as.data.table(myfactors %*% myloadings + myerror, keep.rownames="taxon")
### Standardize the trait data - mean=0 and SD=1
simdat = standardize(simdat0,robust=F)
### Edit and save trait data to file
colnames(simdat)<-c("taxon","Brain(g)", "BL(mm)", "BM(g)", "CBL", "WL(Nm2)", "AR", "WS(mm)", "FA(mm)", "FMature", "MMature", "Litter", "Life", "Forage", "HB(km2)", "Lon(deg)", "Elev(km)", "Realm", "Echo", "Diet", "Hib")
write.csv(simdat,"simdat1.csv", row.names=F)

### Convert some of the variables to qualitative (in this case, binary)
func<-function(x){factor(ifelse(x>median(x),"1","0"))}
simdat2<-simdat %>% mutate(across(c(Litter, Forage, Realm, Echo, Diet, Hib), func))

simdat$Litter<-factor(ifelse(simdat$Litter>median(simdat$Litter),"1","0"))
simdat$Forage<-factor(ifelse(simdat$Forage>median(simdat$Forage),"1","0"))
simdat$Realm<-factor(ifelse(simdat$Realm>median(simdat$Realm),"1","0"))
simdat$Echo<-factor(ifelse(simdat$Echo>median(simdat$Echo),"1","0"))
simdat$Diet<-factor(ifelse(simdat$Diet>median(simdat$Diet),"1","0"))
simdat$Hib<-factor(ifelse(simdat$Hib>median(simdat$Hib),"1","0"))
### Write final binary trait data to file (run again in a new directory to create replicates if desired)
write.csv(simdat2,"simdat2.csv", row.names=F)
```


### **Organizing Results Files**
#### Compiling success measure information from pfa output files

#### Compiling information from simulation runs on Unity:
##### Set to home directory on Unity (cd "directory" for unix code)
### Navigating in R - shrinkage prior run
setwd("/home/moshier.12/Desktop/pfa/pfa_batchrun_ortho/sweet_spot/pfa_newsimtree/pfa_results") 
setwd("/home/moshier.12/Desktop/pfa/pfa_batchrun_ortho/sweet_spot/pfa_newsimtree2/pfa_results")
setwd("/home/moshier.12/Desktop/pfa/pfa_batchrun_ortho/sweet_spot/pfa_simtree50/pfa_results")
setwd("/home/moshier.12/Desktop/pfa/pfa_batchrun_ortho/sweet_spot/pfa_simtree502/pfa_results")

### Navigating in R - iid prior run
setwd("/home/moshier.12/Desktop/pfa/pfa_batchrun_new/newsimtree/pfa_results") 
setwd("/home/moshier.12/Desktop/pfa/pfa_batchrun_new/newsimtree2/pfa_results")
setwd("/home/moshier.12/Desktop/pfa/pfa_batchrun_new/simtree50/pfa_results")
setwd("/home/moshier.12/Desktop/pfa/pfa_batchrun_new/simtree502/pfa_results")
```{r}
#load in necessary libraries
library(tidyr)
library(dplyr)

#load in csv files
loadings<-list.files(path=".", pattern="loadingsStats")
means<-list.files(path=".", pattern="factorMeans")

myfunc<-function(x){
  a = as.data.frame(read.csv(x, header=T))
  iteration = gsub("loadingsStats|.csv","",x)
  
  fac1 = I(list(a$trait[a$factor==1 & a$L>=0.3]))
  fac2 = I(list(a$trait[a$factor==2 & a$L>=0.3]))
  
  fac1L30 = I(list(a$trait[a$factor==1 & a$L>=0.30]))
  fac1L40 = I(list(a$trait[a$factor==1 & a$L>=0.40]))
  fac1L60 = I(list(a$trait[a$factor==1 & a$L>=0.60]))
  
  fac1hpdl = I(list(a$trait[a$factor==1 & a$hpdl>=0]))
  fac1both30 = I(list(a$trait[a$factor==1 & a$L>=0.30 & a$hpdl>=0]))
  fac1both40 = I(list(a$trait[a$factor==1 & a$L>=0.40 & a$hpdl>=0]))
  fac1both60 = I(list(a$trait[a$factor==1 & a$L>=0.60 & a$hpdl>=0]))
  
  fac2L30 = I(list(a$trait[a$factor==2 & a$L>=0.30]))
  fac2L40 = I(list(a$trait[a$factor==2 & a$L>=0.40]))
  fac2L60 = I(list(a$trait[a$factor==2 & a$L>=0.60]))
  
  fac2hpdl = I(list(a$trait[a$factor==2 & a$hpdl>=0]))
  fac2both30 = I(list(a$trait[a$factor==2 & a$L>=0.30 & a$hpdl>=0]))
  fac2both40 = I(list(a$trait[a$factor==2 & a$L>=0.40 & a$hpdl>=0]))
  fac2both60 = I(list(a$trait[a$factor==2 & a$L>=0.60 & a$hpdl>=0]))
  
  mydat_final = data.frame(iteration, fac1, fac2,
                               fac1L30, fac1L40, fac1L60,
                               fac1hpdl, fac1both30, fac1both40, fac1both60, 
                               fac2L30, fac2L40, fac2L60,
                               fac2hpdl, fac2both30, fac2both40, fac2both60) 
  return(mydat_final)
}

mydat_final<-lapply(loadings, FUN = myfunc) %>% bind_rows(.id = NULL)

#write final data and output - these are for shrinkage prior results
saveRDS(mydat_final, file="pfa_newsimtree_results.rds")

saveRDS(mydat_final, file="pfa_newsimtree2_results.rds")

saveRDS(mydat_final, file="pfa_simtree50_results.rds")

saveRDS(mydat_final, file="pfa_simtree502_results.rds")

#write final data and output - these are for iid prior results
saveRDS(mydat_final, file="pfa_newsimtree_resultsiid.rds")

saveRDS(mydat_final, file="pfa_newsimtree2_resultsiid.rds")

saveRDS(mydat_final, file="pfa_simtree50_resultsiid.rds")

saveRDS(mydat_final, file="pfa_simtree502_resultsiid.rds")
```
#### Unix commands for copying the csv files to another directory on Unity for download
cp /home/moshier.12/Desktop/pfa/pfa_batchrun_ortho/sweet_spot/pfa_newsimtree/pfa_results/pfa_newsimtree_results.rds /home/moshier.12/Desktop/pfa/results_filesNov25
cp /home/moshier.12/Desktop/pfa/pfa_batchrun_ortho/sweet_spot/pfa_newsimtree2/pfa_results/pfa_newsimtree2_results.rds /home/moshier.12/Desktop/pfa/results_filesNov25
cp /home/moshier.12/Desktop/pfa/pfa_batchrun_ortho/sweet_spot/pfa_simtree50/pfa_results/pfa_simtree50_results.rds /home/moshier.12/Desktop/pfa/results_filesNov25
cp /home/moshier.12/Desktop/pfa/pfa_batchrun_ortho/sweet_spot/pfa_simtree502/pfa_results/pfa_simtree502_results.rds /home/moshier.12/Desktop/pfa/results_filesNov25

cp /home/moshier.12/Desktop/pfa/pfa_batchrun_new/newsimtree/pfa_results/pfa_newsimtree_resultsiid.rds /home/moshier.12/Desktop/pfa/results_filesNov25
cp /home/moshier.12/Desktop/pfa/pfa_batchrun_new/newsimtree2/pfa_results/pfa_newsimtree2_resultsiid.rds /home/moshier.12/Desktop/pfa/results_filesNov25
cp /home/moshier.12/Desktop/pfa/pfa_batchrun_new/simtree50/pfa_results/pfa_simtree50_resultsiid.rds /home/moshier.12/Desktop/pfa/results_filesNov25
cp /home/moshier.12/Desktop/pfa/pfa_batchrun_new/simtree502/pfa_results/pfa_simtree502_resultsiid.rds /home/moshier.12/Desktop/pfa/results_filesNov25


#### Code for compiling error_sim success measures from results files
```{r}
#load in necessary libraries
library(tidyr)
library(dplyr)

#load in csv files
loadings<-list.files(path=".", pattern="loadingsStats")
means<-list.files(path=".", pattern="factorMeans")

myfunc<-function(x){
  a = as.data.frame(read.csv(x, header=T))
  iteration = gsub("loadingsStats|.csv","",x)
  
  fac1 = I(list(a$trait[a$factor==1 & a$L>=0.3]))
  fac2 = I(list(a$trait[a$factor==2 & a$L>=0.3]))
  
  fac1L30 = I(list(a$trait[a$factor==1 & a$L>=0.30]))
  fac1L40 = I(list(a$trait[a$factor==1 & a$L>=0.40]))
  fac1L60 = I(list(a$trait[a$factor==1 & a$L>=0.60]))
  
  fac1hpdl = I(list(a$trait[a$factor==1 & a$hpdl>=0]))
  fac1both30 = I(list(a$trait[a$factor==1 & a$L>=0.30 & a$hpdl>=0]))
  fac1both40 = I(list(a$trait[a$factor==1 & a$L>=0.40 & a$hpdl>=0]))
  fac1both60 = I(list(a$trait[a$factor==1 & a$L>=0.60 & a$hpdl>=0]))
  
  fac2L30 = I(list(a$trait[a$factor==2 & a$L>=0.30]))
  fac2L40 = I(list(a$trait[a$factor==2 & a$L>=0.40]))
  fac2L60 = I(list(a$trait[a$factor==2 & a$L>=0.60]))
  
  fac2hpdl = I(list(a$trait[a$factor==2 & a$hpdl>=0]))
  fac2both30 = I(list(a$trait[a$factor==2 & a$L>=0.30 & a$hpdl>=0]))
  fac2both40 = I(list(a$trait[a$factor==2 & a$L>=0.40 & a$hpdl>=0]))
  fac2both60 = I(list(a$trait[a$factor==2 & a$L>=0.60 & a$hpdl>=0]))
  
  mydat_final = data.frame(iteration, fac1, fac2,
                               fac1L30, fac1L40, fac1L60,
                               fac1hpdl, fac1both30, fac1both40, fac1both60, 
                               fac2L30, fac2L40, fac2L60,
                               fac2hpdl, fac2both30, fac2both40, fac2both60) 
  return(mydat_final)
}

mydat_final<-lapply(loadings, FUN = myfunc) %>% bind_rows(.id = NULL)

#write final data and output - make sure to be in the proper working directory!
saveRDS(mydat_final, file="error_sim_success.rds")
```

#### Code for compiling discrete_sim success measures from results files
```{r}
#load in necessary libraries
library(tidyr)
library(dplyr)

#load in csv files
loadings<-list.files(path=".", pattern="loadingsStats")
means<-list.files(path=".", pattern="factorMeans")

myfunc<-function(x){
  a = as.data.frame(read.csv(x, header=T))
  iteration = gsub("loadingsStats|.csv","",x)
  
  fac1 = I(list(a$trait[a$factor==1 & a$L>=0.3]))
  fac2 = I(list(a$trait[a$factor==2 & a$L>=0.3]))
  
  fac1L30 = I(list(a$trait[a$factor==1 & a$L>=0.30]))
  fac1L40 = I(list(a$trait[a$factor==1 & a$L>=0.40]))
  fac1L60 = I(list(a$trait[a$factor==1 & a$L>=0.60]))
  
  fac1hpdl = I(list(a$trait[a$factor==1 & a$hpdl>=0]))
  fac1both30 = I(list(a$trait[a$factor==1 & a$L>=0.30 & a$hpdl>=0]))
  fac1both40 = I(list(a$trait[a$factor==1 & a$L>=0.40 & a$hpdl>=0]))
  fac1both60 = I(list(a$trait[a$factor==1 & a$L>=0.60 & a$hpdl>=0]))
  
  fac2L30 = I(list(a$trait[a$factor==2 & a$L>=0.30]))
  fac2L40 = I(list(a$trait[a$factor==2 & a$L>=0.40]))
  fac2L60 = I(list(a$trait[a$factor==2 & a$L>=0.60]))
  
  fac2hpdl = I(list(a$trait[a$factor==2 & a$hpdl>=0]))
  fac2both30 = I(list(a$trait[a$factor==2 & a$L>=0.30 & a$hpdl>=0]))
  fac2both40 = I(list(a$trait[a$factor==2 & a$L>=0.40 & a$hpdl>=0]))
  fac2both60 = I(list(a$trait[a$factor==2 & a$L>=0.60 & a$hpdl>=0]))
  
  hpd_width = I(list(a$hpdu + abs(a$hpdl)))
  
  mydat_final = data.frame(iteration, fac1, fac2,
                               fac1L30, fac1L40, fac1L60,
                               fac1hpdl, fac1both30, fac1both40, fac1both60, 
                               fac2L30, fac2L40, fac2L60,
                               fac2hpdl, fac2both30, fac2both40, fac2both60, hpd_width) 
  return(mydat_final)
}

mydat_final<-lapply(loadings, FUN = myfunc) %>% bind_rows(.id = NULL)

#write final data and output - 20tip runs
saveRDS(mydat_final, file="discrete20_results.rds")

saveRDS(mydat_final, file="discrete202_results.rds")

#write final data and output - 40tip runs
saveRDS(mydat_final, file="discrete40_results.rds")

saveRDS(mydat_final, file="discrete402_results.rds")
```

#### Code for newsimtree - done and fixed!
```{r}
library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)

setwd("C:/Users/shelb/OneDrive/Desktop/Chapter 2 Stuff/Simulation/script files/results_filesNov25") #with results files for both priors

#read in data files
mydat<-read.csv("mydat.csv")

mydat$fac.expec.up<-ceiling(mydat$fac.expec)
mydat$fac.expec.down<-floor(mydat$fac.expec)
#this works for making the list of expected stuff!!!
for (i in 1:length(mydat$ntraits)){
  #mydat$alltraits[[i]]<-paste0("V",1:mydat$ntraits[[i]],collapse=',')
  mydat$alltraits[[i]]<-1:mydat$ntraits[[i]]
  mydat$fac1traits[[i]]<-head(mydat$alltraits[[i]],n=(mydat$fac.expec.down[[i]]))
  mydat$fac2traits[[i]]<-tail(mydat$alltraits[[i]],n=(mydat$fac.expec.down[[i]]))
  mydat$fac1traits[[i]]<-paste0("V",mydat$fac1traits[[i]])
  mydat$fac2traits[[i]]<-paste0("V",mydat$fac2traits[[i]])
} 

#now read in the results files - this is for the shrinkage prior
#result<-readRDS("pfa_newsimtree_results.rds")
#dat1<-merge(mydat, result, by="iteration")

#now read in the results files - this is for the iid prior
result<-readRDS("pfa_newsimtree_resultsiid.rds")
dat1<-merge(mydat, result, by="iteration")

#success measure
for (i in 1:length(dat1$iteration)){
  #fac true <- the proportion of traits we expected in factor 1 vs what we actually got
  #fac label switch <- the traits just got switched between factors
  #sorting just by fac1 and fac2
  dat1$fac1true[[i]]<-dat1$fac1traits[[i]]%in%dat1$fac1[[i]]
  dat1$fac1lblsw[[i]]<-dat1$fac1traits[[i]]%in%dat1$fac2[[i]]
  dat1$fac2true[[i]]<-dat1$fac2traits[[i]]%in%dat1$fac2[[i]]
  dat1$fac2lblsw[[i]]<-dat1$fac2traits[[i]]%in%dat1$fac1[[i]]
  
  #account for traits loading onto multiple factors - fac1 vs fac2
  dat1$fac1overlap[[i]]<-`if`(any(dat1$fac1true[[i]]|dat1$fac1lblsw[[i]]==TRUE),(dat1$fac1true[[i]]&dat1$fac1lblsw[[i]]==TRUE),0)
  dat1$fac2overlap[[i]]<-`if`(any(dat1$fac2true[[i]]|dat1$fac2lblsw[[i]]==TRUE),(dat1$fac2true[[i]]&dat1$fac2lblsw[[i]]==TRUE),0)
  
  dat1$fac1success[[i]]<-mean(dat1$fac1true[[i]])
  dat1$fac1wrong[[i]]<-mean(dat1$fac1overlap[[i]])
  dat1$fac1lblswsuccess[[i]]<-mean(dat1$fac1lblsw[[i]])
  dat1$fac2success[[i]]<-mean(dat1$fac2true[[i]])
  dat1$fac2wrong[[i]]<-mean(dat1$fac2overlap[[i]])
  dat1$fac2lblswsuccess[[i]]<-mean(dat1$fac2lblsw[[i]])
  
  #when L=0.3 (our lowest cutoff)
  dat1$fac1L30true[[i]]<-dat1$fac1traits[[i]]%in%dat1$fac1L30[[i]]
  dat1$fac1L30lblsw[[i]]<-dat1$fac1traits[[i]]%in%dat1$fac2L30[[i]]
  dat1$fac2L30true[[i]]<-dat1$fac2traits[[i]]%in%dat1$fac2L30[[i]]
  dat1$fac2L30lblsw[[i]]<-dat1$fac2traits[[i]]%in%dat1$fac1L30[[i]]
  
  #account for traits loading onto multiple factors - L=0.3
  dat1$fac1L30overlap[[i]]<-`if`(any(dat1$fac1L30true[[i]]|dat1$fac1L30lblsw[[i]]==TRUE), (dat1$fac1L30true[[i]]&dat1$fac1L30lblsw[[i]]==TRUE), 0)
  dat1$fac2L30overlap[[i]]<-`if`(any(dat1$fac2L30true[[i]]|dat1$fac2L30lblsw[[i]]==TRUE), (dat1$fac2L30true[[i]]&dat1$fac2L30lblsw[[i]]==TRUE), 0)
  
  dat1$fac1L30success[[i]]<-mean(dat1$fac1L30true[[i]])
  dat1$fac1L30wrong[[i]]<-mean(dat1$fac1L30overlap[[i]])
  dat1$fac1L30lblswsuccess[[i]]<-mean(dat1$fac1L30lblsw[[i]])
  dat1$fac2L30success[[i]]<-mean(dat1$fac2L30true[[i]])
  dat1$fac2L30wrong[[i]]<-mean(dat1$fac2L30overlap[[i]])
  dat1$fac2L30lblswsuccess[[i]]<-mean(dat1$fac2L30lblsw[[i]])
  
  #when L=0.3 and hpdl>0 (our lowest cutoff and if the posterior dist is over 0)
  dat1$fac1both30true[[i]]<-dat1$fac1traits[[i]]%in%dat1$fac1both30[[i]]
  dat1$fac1both30lblsw[[i]]<-dat1$fac1traits[[i]]%in%dat1$fac2both30[[i]]
  dat1$fac2both30true[[i]]<-dat1$fac2traits[[i]]%in%dat1$fac2both30[[i]]
  dat1$fac2both30lblsw[[i]]<-dat1$fac2traits[[i]]%in%dat1$fac1both30[[i]]
  
  #account for traits loading onto multiple factors - L=0.3 & hpdl>0
  dat1$fac1both30overlap[[i]]<-`if`(any(dat1$fac1both30true[[i]]|dat1$fac1both30lblsw[[i]]==TRUE), (dat1$fac1both30true[[i]]&dat1$fac1both30lblsw[[i]]==TRUE), 0)
  dat1$fac2both30overlap[[i]]<-`if`(any(dat1$fac2both30true[[i]]|dat1$fac2both30lblsw[[i]]==TRUE), (dat1$fac2both30true[[i]]&dat1$fac2both30lblsw[[i]]==TRUE), 0)
  
  dat1$fac1both30success[[i]]<-mean(dat1$fac1both30true[[i]])
  dat1$fac1both30wrong[[i]]<-mean(dat1$fac1both30overlap[[i]])
  dat1$fac1both30lblswsuccess[[i]]<-mean(dat1$fac1both30lblsw[[i]])
  dat1$fac2both30success[[i]]<-mean(dat1$fac2both30true[[i]])
  dat1$fac2both30wrong[[i]]<-mean(dat1$fac2both30overlap[[i]])
  dat1$fac2both30lblswsuccess[[i]]<-mean(dat1$fac2both30lblsw[[i]])
  
  #when L=0.4 (our next lowest cutoff)
  dat1$fac1L40true[[i]]<-dat1$fac1traits[[i]]%in%dat1$fac1L40[[i]]
  dat1$fac1L40lblsw[[i]]<-dat1$fac1traits[[i]]%in%dat1$fac2L40[[i]]
  dat1$fac2L40true[[i]]<-dat1$fac2traits[[i]]%in%dat1$fac2L40[[i]]
  dat1$fac2L40lblsw[[i]]<-dat1$fac2traits[[i]]%in%dat1$fac1L40[[i]]
  
  #account for traits loading onto multiple factors - L=0.4
  dat1$fac1L40overlap[[i]]<-`if`(any(dat1$fac1L40true[[i]]|dat1$fac1L40lblsw[[i]]==TRUE), (dat1$fac1L40true[[i]]&dat1$fac1L40lblsw[[i]]==TRUE), 0)
  dat1$fac2L40overlap[[i]]<-`if`(any(dat1$fac2L40true[[i]]|dat1$fac2L40lblsw[[i]]==TRUE), (dat1$fac2L40true[[i]]&dat1$fac2L40lblsw[[i]]==TRUE), 0)
  
  dat1$fac1L40success[[i]]<-mean(dat1$fac1L40true[[i]])
  dat1$fac1L40wrong[[i]]<-mean(dat1$fac1L40overlap[[i]])
  dat1$fac1L40lblswsuccess[[i]]<-mean(dat1$fac1L40lblsw[[i]])
  dat1$fac2L40success[[i]]<-mean(dat1$fac2L40true[[i]])
  dat1$fac2L40wrong[[i]]<-mean(dat1$fac2L40overlap[[i]])
  dat1$fac2L40lblswsuccess[[i]]<-mean(dat1$fac2L40lblsw[[i]])
  
  #when L=0.4 and hpdl>0 (our next lowest cutoff and if the posterior dist is over 0)
  dat1$fac1both40true[[i]]<-dat1$fac1traits[[i]]%in%dat1$fac1both40[[i]]
  dat1$fac1both40lblsw[[i]]<-dat1$fac1traits[[i]]%in%dat1$fac2both40[[i]]
  dat1$fac2both40true[[i]]<-dat1$fac2traits[[i]]%in%dat1$fac2both40[[i]]
  dat1$fac2both40lblsw[[i]]<-dat1$fac2traits[[i]]%in%dat1$fac1both40[[i]]
  
  #account for traits loading onto multiple factors - L=0.4 & hpdl>0
  dat1$fac1both40overlap[[i]]<-`if`(any(dat1$fac1both40true[[i]]|dat1$fac1both40lblsw[[i]]==TRUE), (dat1$fac1both40true[[i]]&dat1$fac1both40lblsw[[i]]==TRUE), 0)
  dat1$fac2both40overlap[[i]]<-`if`(any(dat1$fac2both40true[[i]]|dat1$fac2both40lblsw[[i]]==TRUE), (dat1$fac2both40true[[i]]&dat1$fac2both40lblsw[[i]]==TRUE), 0)
  
  dat1$fac1both40success[[i]]<-mean(dat1$fac1both40true[[i]])
  dat1$fac1both40wrong[[i]]<-mean(dat1$fac1both40overlap[[i]])
  dat1$fac1both40lblswsuccess[[i]]<-mean(dat1$fac1both40lblsw[[i]])
  dat1$fac2both40success[[i]]<-mean(dat1$fac2both40true[[i]])
  dat1$fac2both40wrong[[i]]<-mean(dat1$fac2both40overlap[[i]])
  dat1$fac2both40lblswsuccess[[i]]<-mean(dat1$fac2both40lblsw[[i]])
  
  #when L=0.6 (our next lowest cutoff)
  dat1$fac1L60true[[i]]<-dat1$fac1traits[[i]]%in%dat1$fac1L60[[i]]
  dat1$fac1L60lblsw[[i]]<-dat1$fac1traits[[i]]%in%dat1$fac2L60[[i]]
  dat1$fac2L60true[[i]]<-dat1$fac2traits[[i]]%in%dat1$fac2L60[[i]]
  dat1$fac2L60lblsw[[i]]<-dat1$fac2traits[[i]]%in%dat1$fac1L60[[i]]
  
  #account for traits loading onto multiple factors - L=0.6
  dat1$fac1L60overlap[[i]]<-`if`(any(dat1$fac1L60true[[i]]|dat1$fac1L60lblsw[[i]]==TRUE), (dat1$fac1L60true[[i]]&dat1$fac1L60lblsw[[i]]==TRUE), 0)
  dat1$fac2L60overlap[[i]]<-`if`(any(dat1$fac2L60true[[i]]|dat1$fac2L60lblsw[[i]]==TRUE), (dat1$fac2L60true[[i]]&dat1$fac2L60lblsw[[i]]==TRUE), 0)
  
  dat1$fac1L60success[[i]]<-mean(dat1$fac1L60true[[i]])
  dat1$fac1L60wrong[[i]]<-mean(dat1$fac1L60overlap[[i]])
  dat1$fac1L60lblswsuccess[[i]]<-mean(dat1$fac1L60lblsw[[i]])
  dat1$fac2L60success[[i]]<-mean(dat1$fac2L60true[[i]])
  dat1$fac2L60wrong[[i]]<-mean(dat1$fac2L60overlap[[i]])
  dat1$fac2L60lblswsuccess[[i]]<-mean(dat1$fac2L60lblsw[[i]])
  
  #when L=0.6 and hpdl>0 (our next lowest cutoff and if the posterior dist is over 0)
  dat1$fac1both60true[[i]]<-dat1$fac1traits[[i]]%in%dat1$fac1both60[[i]]
  dat1$fac1both60lblsw[[i]]<-dat1$fac1traits[[i]]%in%dat1$fac2both60[[i]]
  dat1$fac2both60true[[i]]<-dat1$fac2traits[[i]]%in%dat1$fac2both60[[i]]
  dat1$fac2both60lblsw[[i]]<-dat1$fac2traits[[i]]%in%dat1$fac1both60[[i]]
  
  #account for traits loading onto multiple factors - L=0.4 & hpdl>0
  dat1$fac1both60overlap[[i]]<-`if`(any(dat1$fac1both60true[[i]]|dat1$fac1both60lblsw[[i]]==TRUE), (dat1$fac1both60true[[i]]&dat1$fac1both60lblsw[[i]]==TRUE), 0)
  dat1$fac2both60overlap[[i]]<-`if`(any(dat1$fac2both60true[[i]]|dat1$fac2both60lblsw[[i]]==TRUE), (dat1$fac2both60true[[i]]&dat1$fac2both60lblsw[[i]]==TRUE), 0)
  
  dat1$fac1both60success[[i]]<-mean(dat1$fac1both60true[[i]])
  dat1$fac1both60wrong[[i]]<-mean(dat1$fac1both60overlap[[i]])
  dat1$fac1both60lblswsuccess[[i]]<-mean(dat1$fac1both60lblsw[[i]])
  dat1$fac2both60success[[i]]<-mean(dat1$fac2both60true[[i]])
  dat1$fac2both60wrong[[i]]<-mean(dat1$fac2both60overlap[[i]])
  dat1$fac2both60lblswsuccess[[i]]<-mean(dat1$fac2both60lblsw[[i]])
}

#make sure the success columns are numeric & replace NaN with 0
x<-c(33:38, 45:50, 57:62, 69:74, 81:86, 93:98, 105:110)
dat1[,x]<-lapply(dat1[,x], as.numeric)
dat1[is.na(dat1)]<-0 #need to do is.na when it's a data frame

#the overall number of things labeled either fac1 or fac2
dat1$fac1overall<-dat1$fac1success + dat1$fac1lblswsuccess - 2*dat1$fac1wrong
dat1$fac2overall<-dat1$fac2success + dat1$fac2lblswsuccess - 2*dat1$fac2wrong

#when L=0.3 (our next lowest cutoff)
dat1$fac1L30overall<-dat1$fac1L30success + dat1$fac1L30lblswsuccess - 2*dat1$fac1L30wrong
dat1$fac2L30overall<-dat1$fac2L30success + dat1$fac2L30lblswsuccess - 2*dat1$fac2L30wrong

#when L=0.3 and hpdl>0 (our next lowest cutoff and if the posterior dist is over 0)
dat1$fac1both30overall<-dat1$fac1both30success + dat1$fac1both30lblswsuccess - 2*dat1$fac1both30wrong
dat1$fac2both30overall<-dat1$fac2both30success + dat1$fac2both30lblswsuccess - 2*dat1$fac2both30wrong

#when L=0.4 (our next lowest cutoff)
dat1$fac1L40overall<-dat1$fac1L40success + dat1$fac1L40lblswsuccess - 2*dat1$fac1L40wrong
dat1$fac2L40overall<-dat1$fac2L40success + dat1$fac2L40lblswsuccess - 2*dat1$fac2L40wrong

#when L=0.4 and hpdl>0 (our next lowest cutoff and if the posterior dist is over 0)
dat1$fac1both40overall<-dat1$fac1both40success + dat1$fac1both40lblswsuccess - 2*dat1$fac1both40wrong
dat1$fac2both40overall<-dat1$fac2both40success + dat1$fac2both40lblswsuccess - 2*dat1$fac2both40wrong

#when L=0.6 (our next lowest cutoff)
dat1$fac1L60overall<-dat1$fac1L60success + dat1$fac1L60lblswsuccess - 2*dat1$fac1L60wrong
dat1$fac2L60overall<-dat1$fac2L60success + dat1$fac2L60lblswsuccess - 2*dat1$fac2L60wrong

#when L=0.6 and hpdl>0 (our next lowest cutoff and if the posterior dist is over 0)
dat1$fac1both60overall<-dat1$fac1both60success + dat1$fac1both60lblswsuccess - 2*dat1$fac1both60wrong
dat1$fac2both60overall<-dat1$fac2both60success + dat1$fac2both60lblswsuccess - 2*dat1$fac2both60wrong

#add successes together for an overall measure between 0-2
dat1$facoverall<-dat1$fac1overall + dat1$fac2overall
dat1$fac30overall<-dat1$fac1L30overall + dat1$fac2L30overall
dat1$facboth30overall<-dat1$fac1both30overall + dat1$fac2both30overall
dat1$fac40overall<-dat1$fac1L40overall + dat1$fac2L40overall
dat1$facboth40overall<-dat1$fac1both40overall + dat1$fac2both40overall
dat1$fac60overall<-dat1$fac1L60overall + dat1$fac2L60overall
dat1$facboth60overall<-dat1$fac1both60overall + dat1$fac2both60overall

#use custom function to scale data columns to between 0-1
minmax <- function(x) {
    return((x- min(x)) /(max(x)-min(x)))}
dat1[,125:131] <- minmax(dat1[,125:131])

#write final data to file
#saveRDS(dat1, file="dat1_shrink.rds") #shrinkage prior
saveRDS(dat1, file="dat1_iid.rds") #iid prior
```

#### Code for newsimtree2 - done and fixed!
```{r}
library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)

setwd("C:/Users/shelb/OneDrive/Desktop/Chapter 2 Stuff/Simulation/script files/results_filesNov25") #with results files for both priors

#read in data files
mydat<-read.csv("mydat.csv")

mydat$fac.expec.up<-ceiling(mydat$fac.expec)
mydat$fac.expec.down<-floor(mydat$fac.expec)
#this works for making the list of expected stuff!!!
for (i in 1:length(mydat$ntraits)){
  #mydat$alltraits[[i]]<-paste0("V",1:mydat$ntraits[[i]],collapse=',')
  mydat$alltraits[[i]]<-1:mydat$ntraits[[i]]
  mydat$fac1traits[[i]]<-head(mydat$alltraits[[i]],n=(mydat$fac.expec.down[[i]]))
  mydat$fac2traits[[i]]<-tail(mydat$alltraits[[i]],n=(mydat$fac.expec.down[[i]]))
  mydat$fac1traits[[i]]<-paste0("V",mydat$fac1traits[[i]])
  mydat$fac2traits[[i]]<-paste0("V",mydat$fac2traits[[i]])
} 

#now read in the results files - this is for the shrinkage prior
result<-readRDS("pfa_newsimtree2_results.rds")
dat2<-merge(mydat, result, by="iteration")

#now read in the results files - this is for the iid prior
#result<-readRDS("pfa_newsimtree2_resultsiid.rds")
#dat2<-merge(mydat, result, by="iteration")

#success measure
for (i in 1:length(dat2$iteration)){
  #fac true <- the proportion of traits we expected in factor 1 vs what we actually got
  #fac label switch <- the traits just got switched between factors
  #sorting just by fac1 and fac2
  dat2$fac1true[[i]]<-dat2$fac1traits[[i]]%in%dat2$fac1[[i]]
  dat2$fac1lblsw[[i]]<-dat2$fac1traits[[i]]%in%dat2$fac2[[i]]
  dat2$fac2true[[i]]<-dat2$fac2traits[[i]]%in%dat2$fac2[[i]]
  dat2$fac2lblsw[[i]]<-dat2$fac2traits[[i]]%in%dat2$fac1[[i]]
  
  #account for traits loading onto multiple factors - fac1 vs fac2
  dat2$fac1overlap[[i]]<-`if`(any(dat2$fac1true[[i]]|dat2$fac1lblsw[[i]]==TRUE),(dat2$fac1true[[i]]&dat2$fac1lblsw[[i]]==TRUE),0)
  dat2$fac2overlap[[i]]<-`if`(any(dat2$fac2true[[i]]|dat2$fac2lblsw[[i]]==TRUE),(dat2$fac2true[[i]]&dat2$fac2lblsw[[i]]==TRUE),0)
  
  dat2$fac1success[[i]]<-mean(dat2$fac1true[[i]])
  dat2$fac1wrong[[i]]<-mean(dat2$fac1overlap[[i]])
  dat2$fac1lblswsuccess[[i]]<-mean(dat2$fac1lblsw[[i]])
  dat2$fac2success[[i]]<-mean(dat2$fac2true[[i]])
  dat2$fac2wrong[[i]]<-mean(dat2$fac2overlap[[i]])
  dat2$fac2lblswsuccess[[i]]<-mean(dat2$fac2lblsw[[i]])
  
  #when L=0.3 (our lowest cutoff)
  dat2$fac1L30true[[i]]<-dat2$fac1traits[[i]]%in%dat2$fac1L30[[i]]
  dat2$fac1L30lblsw[[i]]<-dat2$fac1traits[[i]]%in%dat2$fac2L30[[i]]
  dat2$fac2L30true[[i]]<-dat2$fac2traits[[i]]%in%dat2$fac2L30[[i]]
  dat2$fac2L30lblsw[[i]]<-dat2$fac2traits[[i]]%in%dat2$fac1L30[[i]]
  
  #account for traits loading onto multiple factors - L=0.3
  dat2$fac1L30overlap[[i]]<-`if`(any(dat2$fac1L30true[[i]]|dat2$fac1L30lblsw[[i]]==TRUE), (dat2$fac1L30true[[i]]&dat2$fac1L30lblsw[[i]]==TRUE), 0)
  dat2$fac2L30overlap[[i]]<-`if`(any(dat2$fac2L30true[[i]]|dat2$fac2L30lblsw[[i]]==TRUE), (dat2$fac2L30true[[i]]&dat2$fac2L30lblsw[[i]]==TRUE), 0)
  
  dat2$fac1L30success[[i]]<-mean(dat2$fac1L30true[[i]])
  dat2$fac1L30wrong[[i]]<-mean(dat2$fac1L30overlap[[i]])
  dat2$fac1L30lblswsuccess[[i]]<-mean(dat2$fac1L30lblsw[[i]])
  dat2$fac2L30success[[i]]<-mean(dat2$fac2L30true[[i]])
  dat2$fac2L30wrong[[i]]<-mean(dat2$fac2L30overlap[[i]])
  dat2$fac2L30lblswsuccess[[i]]<-mean(dat2$fac2L30lblsw[[i]])
  
  #when L=0.3 and hpdl>0 (our lowest cutoff and if the posterior dist is over 0)
  dat2$fac1both30true[[i]]<-dat2$fac1traits[[i]]%in%dat2$fac1both30[[i]]
  dat2$fac1both30lblsw[[i]]<-dat2$fac1traits[[i]]%in%dat2$fac2both30[[i]]
  dat2$fac2both30true[[i]]<-dat2$fac2traits[[i]]%in%dat2$fac2both30[[i]]
  dat2$fac2both30lblsw[[i]]<-dat2$fac2traits[[i]]%in%dat2$fac1both30[[i]]
  
  #account for traits loading onto multiple factors - L=0.3 & hpdl>0
  dat2$fac1both30overlap[[i]]<-`if`(any(dat2$fac1both30true[[i]]|dat2$fac1both30lblsw[[i]]==TRUE), (dat2$fac1both30true[[i]]&dat2$fac1both30lblsw[[i]]==TRUE), 0)
  dat2$fac2both30overlap[[i]]<-`if`(any(dat2$fac2both30true[[i]]|dat2$fac2both30lblsw[[i]]==TRUE), (dat2$fac2both30true[[i]]&dat2$fac2both30lblsw[[i]]==TRUE), 0)
  
  dat2$fac1both30success[[i]]<-mean(dat2$fac1both30true[[i]])
  dat2$fac1both30wrong[[i]]<-mean(dat2$fac1both30overlap[[i]])
  dat2$fac1both30lblswsuccess[[i]]<-mean(dat2$fac1both30lblsw[[i]])
  dat2$fac2both30success[[i]]<-mean(dat2$fac2both30true[[i]])
  dat2$fac2both30wrong[[i]]<-mean(dat2$fac2both30overlap[[i]])
  dat2$fac2both30lblswsuccess[[i]]<-mean(dat2$fac2both30lblsw[[i]])
  
  #when L=0.4 (our next lowest cutoff)
  dat2$fac1L40true[[i]]<-dat2$fac1traits[[i]]%in%dat2$fac1L40[[i]]
  dat2$fac1L40lblsw[[i]]<-dat2$fac1traits[[i]]%in%dat2$fac2L40[[i]]
  dat2$fac2L40true[[i]]<-dat2$fac2traits[[i]]%in%dat2$fac2L40[[i]]
  dat2$fac2L40lblsw[[i]]<-dat2$fac2traits[[i]]%in%dat2$fac1L40[[i]]
  
  #account for traits loading onto multiple factors - L=0.4
  dat2$fac1L40overlap[[i]]<-`if`(any(dat2$fac1L40true[[i]]|dat2$fac1L40lblsw[[i]]==TRUE), (dat2$fac1L40true[[i]]&dat2$fac1L40lblsw[[i]]==TRUE), 0)
  dat2$fac2L40overlap[[i]]<-`if`(any(dat2$fac2L40true[[i]]|dat2$fac2L40lblsw[[i]]==TRUE), (dat2$fac2L40true[[i]]&dat2$fac2L40lblsw[[i]]==TRUE), 0)
  
  dat2$fac1L40success[[i]]<-mean(dat2$fac1L40true[[i]])
  dat2$fac1L40wrong[[i]]<-mean(dat2$fac1L40overlap[[i]])
  dat2$fac1L40lblswsuccess[[i]]<-mean(dat2$fac1L40lblsw[[i]])
  dat2$fac2L40success[[i]]<-mean(dat2$fac2L40true[[i]])
  dat2$fac2L40wrong[[i]]<-mean(dat2$fac2L40overlap[[i]])
  dat2$fac2L40lblswsuccess[[i]]<-mean(dat2$fac2L40lblsw[[i]])
  
  #when L=0.4 and hpdl>0 (our next lowest cutoff and if the posterior dist is over 0)
  dat2$fac1both40true[[i]]<-dat2$fac1traits[[i]]%in%dat2$fac1both40[[i]]
  dat2$fac1both40lblsw[[i]]<-dat2$fac1traits[[i]]%in%dat2$fac2both40[[i]]
  dat2$fac2both40true[[i]]<-dat2$fac2traits[[i]]%in%dat2$fac2both40[[i]]
  dat2$fac2both40lblsw[[i]]<-dat2$fac2traits[[i]]%in%dat2$fac1both40[[i]]
  
  #account for traits loading onto multiple factors - L=0.4 & hpdl>0
  dat2$fac1both40overlap[[i]]<-`if`(any(dat2$fac1both40true[[i]]|dat2$fac1both40lblsw[[i]]==TRUE), (dat2$fac1both40true[[i]]&dat2$fac1both40lblsw[[i]]==TRUE), 0)
  dat2$fac2both40overlap[[i]]<-`if`(any(dat2$fac2both40true[[i]]|dat2$fac2both40lblsw[[i]]==TRUE), (dat2$fac2both40true[[i]]&dat2$fac2both40lblsw[[i]]==TRUE), 0)
  
  dat2$fac1both40success[[i]]<-mean(dat2$fac1both40true[[i]])
  dat2$fac1both40wrong[[i]]<-mean(dat2$fac1both40overlap[[i]])
  dat2$fac1both40lblswsuccess[[i]]<-mean(dat2$fac1both40lblsw[[i]])
  dat2$fac2both40success[[i]]<-mean(dat2$fac2both40true[[i]])
  dat2$fac2both40wrong[[i]]<-mean(dat2$fac2both40overlap[[i]])
  dat2$fac2both40lblswsuccess[[i]]<-mean(dat2$fac2both40lblsw[[i]])
  
  #when L=0.6 (our next lowest cutoff)
  dat2$fac1L60true[[i]]<-dat2$fac1traits[[i]]%in%dat2$fac1L60[[i]]
  dat2$fac1L60lblsw[[i]]<-dat2$fac1traits[[i]]%in%dat2$fac2L60[[i]]
  dat2$fac2L60true[[i]]<-dat2$fac2traits[[i]]%in%dat2$fac2L60[[i]]
  dat2$fac2L60lblsw[[i]]<-dat2$fac2traits[[i]]%in%dat2$fac1L60[[i]]
  
  #account for traits loading onto multiple factors - L=0.6
  dat2$fac1L60overlap[[i]]<-`if`(any(dat2$fac1L60true[[i]]|dat2$fac1L60lblsw[[i]]==TRUE), (dat2$fac1L60true[[i]]&dat2$fac1L60lblsw[[i]]==TRUE), 0)
  dat2$fac2L60overlap[[i]]<-`if`(any(dat2$fac2L60true[[i]]|dat2$fac2L60lblsw[[i]]==TRUE), (dat2$fac2L60true[[i]]&dat2$fac2L60lblsw[[i]]==TRUE), 0)
  
  dat2$fac1L60success[[i]]<-mean(dat2$fac1L60true[[i]])
  dat2$fac1L60wrong[[i]]<-mean(dat2$fac1L60overlap[[i]])
  dat2$fac1L60lblswsuccess[[i]]<-mean(dat2$fac1L60lblsw[[i]])
  dat2$fac2L60success[[i]]<-mean(dat2$fac2L60true[[i]])
  dat2$fac2L60wrong[[i]]<-mean(dat2$fac2L60overlap[[i]])
  dat2$fac2L60lblswsuccess[[i]]<-mean(dat2$fac2L60lblsw[[i]])
  
  #when L=0.6 and hpdl>0 (our next lowest cutoff and if the posterior dist is over 0)
  dat2$fac1both60true[[i]]<-dat2$fac1traits[[i]]%in%dat2$fac1both60[[i]]
  dat2$fac1both60lblsw[[i]]<-dat2$fac1traits[[i]]%in%dat2$fac2both60[[i]]
  dat2$fac2both60true[[i]]<-dat2$fac2traits[[i]]%in%dat2$fac2both60[[i]]
  dat2$fac2both60lblsw[[i]]<-dat2$fac2traits[[i]]%in%dat2$fac1both60[[i]]
  
  #account for traits loading onto multiple factors - L=0.4 & hpdl>0
  dat2$fac1both60overlap[[i]]<-`if`(any(dat2$fac1both60true[[i]]|dat2$fac1both60lblsw[[i]]==TRUE), (dat2$fac1both60true[[i]]&dat2$fac1both60lblsw[[i]]==TRUE), 0)
  dat2$fac2both60overlap[[i]]<-`if`(any(dat2$fac2both60true[[i]]|dat2$fac2both60lblsw[[i]]==TRUE), (dat2$fac2both60true[[i]]&dat2$fac2both60lblsw[[i]]==TRUE), 0)
  
  dat2$fac1both60success[[i]]<-mean(dat2$fac1both60true[[i]])
  dat2$fac1both60wrong[[i]]<-mean(dat2$fac1both60overlap[[i]])
  dat2$fac1both60lblswsuccess[[i]]<-mean(dat2$fac1both60lblsw[[i]])
  dat2$fac2both60success[[i]]<-mean(dat2$fac2both60true[[i]])
  dat2$fac2both60wrong[[i]]<-mean(dat2$fac2both60overlap[[i]])
  dat2$fac2both60lblswsuccess[[i]]<-mean(dat2$fac2both60lblsw[[i]])
}

#make sure the success columns are numeric & replace NaN with 0
x<-c(33:38, 45:50, 57:62, 69:74, 81:86, 93:98, 105:110)
dat2[,x]<-lapply(dat2[,x], as.numeric)
dat2[is.na(dat2)]<-0 #need to do is.na when it's a data frame

#the overall number of things labeled either fac1 or fac2
dat2$fac1overall<-dat2$fac1success + dat2$fac1lblswsuccess - 2*dat2$fac1wrong
dat2$fac2overall<-dat2$fac2success + dat2$fac2lblswsuccess - 2*dat2$fac2wrong

#when L=0.3 (our next lowest cutoff)
dat2$fac1L30overall<-dat2$fac1L30success + dat2$fac1L30lblswsuccess - 2*dat2$fac1L30wrong
dat2$fac2L30overall<-dat2$fac2L30success + dat2$fac2L30lblswsuccess - 2*dat2$fac2L30wrong

#when L=0.3 and hpdl>0 (our next lowest cutoff and if the posterior dist is over 0)
dat2$fac1both30overall<-dat2$fac1both30success + dat2$fac1both30lblswsuccess - 2*dat2$fac1both30wrong
dat2$fac2both30overall<-dat2$fac2both30success + dat2$fac2both30lblswsuccess - 2*dat2$fac2both30wrong

#when L=0.4 (our next lowest cutoff)
dat2$fac1L40overall<-dat2$fac1L40success + dat2$fac1L40lblswsuccess - 2*dat2$fac1L40wrong
dat2$fac2L40overall<-dat2$fac2L40success + dat2$fac2L40lblswsuccess - 2*dat2$fac2L40wrong

#when L=0.4 and hpdl>0 (our next lowest cutoff and if the posterior dist is over 0)
dat2$fac1both40overall<-dat2$fac1both40success + dat2$fac1both40lblswsuccess - 2*dat2$fac1both40wrong
dat2$fac2both40overall<-dat2$fac2both40success + dat2$fac2both40lblswsuccess - 2*dat2$fac2both40wrong

#when L=0.6 (our next lowest cutoff)
dat2$fac1L60overall<-dat2$fac1L60success + dat2$fac1L60lblswsuccess - 2*dat2$fac1L60wrong
dat2$fac2L60overall<-dat2$fac2L60success + dat2$fac2L60lblswsuccess - 2*dat2$fac2L60wrong

#when L=0.6 and hpdl>0 (our next lowest cutoff and if the posterior dist is over 0)
dat2$fac1both60overall<-dat2$fac1both60success + dat2$fac1both60lblswsuccess - 2*dat2$fac1both60wrong
dat2$fac2both60overall<-dat2$fac2both60success + dat2$fac2both60lblswsuccess - 2*dat2$fac2both60wrong

#add successes together for an overall measure between 0-2
dat2$facoverall<-dat2$fac1overall + dat2$fac2overall
dat2$fac30overall<-dat2$fac1L30overall + dat2$fac2L30overall
dat2$facboth30overall<-dat2$fac1both30overall + dat2$fac2both30overall
dat2$fac40overall<-dat2$fac1L40overall + dat2$fac2L40overall
dat2$facboth40overall<-dat2$fac1both40overall + dat2$fac2both40overall
dat2$fac60overall<-dat2$fac1L60overall + dat2$fac2L60overall
dat2$facboth60overall<-dat2$fac1both60overall + dat2$fac2both60overall

#use custom function to scale data columns to between 0-1
minmax <- function(x) {
    return((x- min(x)) /(max(x)-min(x)))}
dat2[,125:131] <- minmax(dat2[,125:131])

#write final data to file
saveRDS(dat2, file="dat2_shrink.rds") #shrinkage prior
#saveRDS(dat2, file="dat2_iid.rds") #iid prior
```

#### Code for simtree50 - done and fixed!
```{r}
library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)

setwd("C:/Users/shelb/OneDrive/Desktop/Chapter 2 Stuff/Simulation/script files/results_filesNov25") #with results files for both priors

#read in data files
mydat<-read.csv("mydat.csv")

mydat$fac.expec.up<-ceiling(mydat$fac.expec)
mydat$fac.expec.down<-floor(mydat$fac.expec)
#this works for making the list of expected stuff!!!
for (i in 1:length(mydat$ntraits)){
  #mydat$alltraits[[i]]<-paste0("V",1:mydat$ntraits[[i]],collapse=',')
  mydat$alltraits[[i]]<-1:mydat$ntraits[[i]]
  mydat$fac1traits[[i]]<-head(mydat$alltraits[[i]],n=(mydat$fac.expec.down[[i]]))
  mydat$fac2traits[[i]]<-tail(mydat$alltraits[[i]],n=(mydat$fac.expec.down[[i]]))
  mydat$fac1traits[[i]]<-paste0("V",mydat$fac1traits[[i]])
  mydat$fac2traits[[i]]<-paste0("V",mydat$fac2traits[[i]])
} 

#now read in the results files - this is for the shrinkage prior
#result<-readRDS("pfa_simtree50_results.rds")
#dat3<-merge(mydat, result, by="iteration")

#now read in the results files - this is for the iid prior
result<-readRDS("pfa_simtree50_resultsiid.rds")
dat3<-merge(mydat, result, by="iteration")

#success measure
for (i in 1:length(dat3$iteration)){
  #fac true <- the proportion of traits we expected in factor 1 vs what we actually got
  #fac label switch <- the traits just got switched between factors
  #sorting just by fac1 and fac2
  dat3$fac1true[[i]]<-dat3$fac1traits[[i]]%in%dat3$fac1[[i]]
  dat3$fac1lblsw[[i]]<-dat3$fac1traits[[i]]%in%dat3$fac2[[i]]
  dat3$fac2true[[i]]<-dat3$fac2traits[[i]]%in%dat3$fac2[[i]]
  dat3$fac2lblsw[[i]]<-dat3$fac2traits[[i]]%in%dat3$fac1[[i]]
  
  #account for traits loading onto multiple factors - fac1 vs fac2
  dat3$fac1overlap[[i]]<-`if`(any(dat3$fac1true[[i]]|dat3$fac1lblsw[[i]]==TRUE),(dat3$fac1true[[i]]&dat3$fac1lblsw[[i]]==TRUE),0)
  dat3$fac2overlap[[i]]<-`if`(any(dat3$fac2true[[i]]|dat3$fac2lblsw[[i]]==TRUE),(dat3$fac2true[[i]]&dat3$fac2lblsw[[i]]==TRUE),0)
  
  dat3$fac1success[[i]]<-mean(dat3$fac1true[[i]])
  dat3$fac1wrong[[i]]<-mean(dat3$fac1overlap[[i]])
  dat3$fac1lblswsuccess[[i]]<-mean(dat3$fac1lblsw[[i]])
  dat3$fac2success[[i]]<-mean(dat3$fac2true[[i]])
  dat3$fac2wrong[[i]]<-mean(dat3$fac2overlap[[i]])
  dat3$fac2lblswsuccess[[i]]<-mean(dat3$fac2lblsw[[i]])
  
  #when L=0.3 (our lowest cutoff)
  dat3$fac1L30true[[i]]<-dat3$fac1traits[[i]]%in%dat3$fac1L30[[i]]
  dat3$fac1L30lblsw[[i]]<-dat3$fac1traits[[i]]%in%dat3$fac2L30[[i]]
  dat3$fac2L30true[[i]]<-dat3$fac2traits[[i]]%in%dat3$fac2L30[[i]]
  dat3$fac2L30lblsw[[i]]<-dat3$fac2traits[[i]]%in%dat3$fac1L30[[i]]
  
  #account for traits loading onto multiple factors - L=0.3
  dat3$fac1L30overlap[[i]]<-`if`(any(dat3$fac1L30true[[i]]|dat3$fac1L30lblsw[[i]]==TRUE), (dat3$fac1L30true[[i]]&dat3$fac1L30lblsw[[i]]==TRUE), 0)
  dat3$fac2L30overlap[[i]]<-`if`(any(dat3$fac2L30true[[i]]|dat3$fac2L30lblsw[[i]]==TRUE), (dat3$fac2L30true[[i]]&dat3$fac2L30lblsw[[i]]==TRUE), 0)
  
  dat3$fac1L30success[[i]]<-mean(dat3$fac1L30true[[i]])
  dat3$fac1L30wrong[[i]]<-mean(dat3$fac1L30overlap[[i]])
  dat3$fac1L30lblswsuccess[[i]]<-mean(dat3$fac1L30lblsw[[i]])
  dat3$fac2L30success[[i]]<-mean(dat3$fac2L30true[[i]])
  dat3$fac2L30wrong[[i]]<-mean(dat3$fac2L30overlap[[i]])
  dat3$fac2L30lblswsuccess[[i]]<-mean(dat3$fac2L30lblsw[[i]])
  
  #when L=0.3 and hpdl>0 (our lowest cutoff and if the posterior dist is over 0)
  dat3$fac1both30true[[i]]<-dat3$fac1traits[[i]]%in%dat3$fac1both30[[i]]
  dat3$fac1both30lblsw[[i]]<-dat3$fac1traits[[i]]%in%dat3$fac2both30[[i]]
  dat3$fac2both30true[[i]]<-dat3$fac2traits[[i]]%in%dat3$fac2both30[[i]]
  dat3$fac2both30lblsw[[i]]<-dat3$fac2traits[[i]]%in%dat3$fac1both30[[i]]
  
  #account for traits loading onto multiple factors - L=0.3 & hpdl>0
  dat3$fac1both30overlap[[i]]<-`if`(any(dat3$fac1both30true[[i]]|dat3$fac1both30lblsw[[i]]==TRUE), (dat3$fac1both30true[[i]]&dat3$fac1both30lblsw[[i]]==TRUE), 0)
  dat3$fac2both30overlap[[i]]<-`if`(any(dat3$fac2both30true[[i]]|dat3$fac2both30lblsw[[i]]==TRUE), (dat3$fac2both30true[[i]]&dat3$fac2both30lblsw[[i]]==TRUE), 0)
  
  dat3$fac1both30success[[i]]<-mean(dat3$fac1both30true[[i]])
  dat3$fac1both30wrong[[i]]<-mean(dat3$fac1both30overlap[[i]])
  dat3$fac1both30lblswsuccess[[i]]<-mean(dat3$fac1both30lblsw[[i]])
  dat3$fac2both30success[[i]]<-mean(dat3$fac2both30true[[i]])
  dat3$fac2both30wrong[[i]]<-mean(dat3$fac2both30overlap[[i]])
  dat3$fac2both30lblswsuccess[[i]]<-mean(dat3$fac2both30lblsw[[i]])
  
  #when L=0.4 (our next lowest cutoff)
  dat3$fac1L40true[[i]]<-dat3$fac1traits[[i]]%in%dat3$fac1L40[[i]]
  dat3$fac1L40lblsw[[i]]<-dat3$fac1traits[[i]]%in%dat3$fac2L40[[i]]
  dat3$fac2L40true[[i]]<-dat3$fac2traits[[i]]%in%dat3$fac2L40[[i]]
  dat3$fac2L40lblsw[[i]]<-dat3$fac2traits[[i]]%in%dat3$fac1L40[[i]]
  
  #account for traits loading onto multiple factors - L=0.4
  dat3$fac1L40overlap[[i]]<-`if`(any(dat3$fac1L40true[[i]]|dat3$fac1L40lblsw[[i]]==TRUE), (dat3$fac1L40true[[i]]&dat3$fac1L40lblsw[[i]]==TRUE), 0)
  dat3$fac2L40overlap[[i]]<-`if`(any(dat3$fac2L40true[[i]]|dat3$fac2L40lblsw[[i]]==TRUE), (dat3$fac2L40true[[i]]&dat3$fac2L40lblsw[[i]]==TRUE), 0)
  
  dat3$fac1L40success[[i]]<-mean(dat3$fac1L40true[[i]])
  dat3$fac1L40wrong[[i]]<-mean(dat3$fac1L40overlap[[i]])
  dat3$fac1L40lblswsuccess[[i]]<-mean(dat3$fac1L40lblsw[[i]])
  dat3$fac2L40success[[i]]<-mean(dat3$fac2L40true[[i]])
  dat3$fac2L40wrong[[i]]<-mean(dat3$fac2L40overlap[[i]])
  dat3$fac2L40lblswsuccess[[i]]<-mean(dat3$fac2L40lblsw[[i]])
  
  #when L=0.4 and hpdl>0 (our next lowest cutoff and if the posterior dist is over 0)
  dat3$fac1both40true[[i]]<-dat3$fac1traits[[i]]%in%dat3$fac1both40[[i]]
  dat3$fac1both40lblsw[[i]]<-dat3$fac1traits[[i]]%in%dat3$fac2both40[[i]]
  dat3$fac2both40true[[i]]<-dat3$fac2traits[[i]]%in%dat3$fac2both40[[i]]
  dat3$fac2both40lblsw[[i]]<-dat3$fac2traits[[i]]%in%dat3$fac1both40[[i]]
  
  #account for traits loading onto multiple factors - L=0.4 & hpdl>0
  dat3$fac1both40overlap[[i]]<-`if`(any(dat3$fac1both40true[[i]]|dat3$fac1both40lblsw[[i]]==TRUE), (dat3$fac1both40true[[i]]&dat3$fac1both40lblsw[[i]]==TRUE), 0)
  dat3$fac2both40overlap[[i]]<-`if`(any(dat3$fac2both40true[[i]]|dat3$fac2both40lblsw[[i]]==TRUE), (dat3$fac2both40true[[i]]&dat3$fac2both40lblsw[[i]]==TRUE), 0)
  
  dat3$fac1both40success[[i]]<-mean(dat3$fac1both40true[[i]])
  dat3$fac1both40wrong[[i]]<-mean(dat3$fac1both40overlap[[i]])
  dat3$fac1both40lblswsuccess[[i]]<-mean(dat3$fac1both40lblsw[[i]])
  dat3$fac2both40success[[i]]<-mean(dat3$fac2both40true[[i]])
  dat3$fac2both40wrong[[i]]<-mean(dat3$fac2both40overlap[[i]])
  dat3$fac2both40lblswsuccess[[i]]<-mean(dat3$fac2both40lblsw[[i]])
  
  #when L=0.6 (our next lowest cutoff)
  dat3$fac1L60true[[i]]<-dat3$fac1traits[[i]]%in%dat3$fac1L60[[i]]
  dat3$fac1L60lblsw[[i]]<-dat3$fac1traits[[i]]%in%dat3$fac2L60[[i]]
  dat3$fac2L60true[[i]]<-dat3$fac2traits[[i]]%in%dat3$fac2L60[[i]]
  dat3$fac2L60lblsw[[i]]<-dat3$fac2traits[[i]]%in%dat3$fac1L60[[i]]
  
  #account for traits loading onto multiple factors - L=0.6
  dat3$fac1L60overlap[[i]]<-`if`(any(dat3$fac1L60true[[i]]|dat3$fac1L60lblsw[[i]]==TRUE), (dat3$fac1L60true[[i]]&dat3$fac1L60lblsw[[i]]==TRUE), 0)
  dat3$fac2L60overlap[[i]]<-`if`(any(dat3$fac2L60true[[i]]|dat3$fac2L60lblsw[[i]]==TRUE), (dat3$fac2L60true[[i]]&dat3$fac2L60lblsw[[i]]==TRUE), 0)
  
  dat3$fac1L60success[[i]]<-mean(dat3$fac1L60true[[i]])
  dat3$fac1L60wrong[[i]]<-mean(dat3$fac1L60overlap[[i]])
  dat3$fac1L60lblswsuccess[[i]]<-mean(dat3$fac1L60lblsw[[i]])
  dat3$fac2L60success[[i]]<-mean(dat3$fac2L60true[[i]])
  dat3$fac2L60wrong[[i]]<-mean(dat3$fac2L60overlap[[i]])
  dat3$fac2L60lblswsuccess[[i]]<-mean(dat3$fac2L60lblsw[[i]])
  
  #when L=0.6 and hpdl>0 (our next lowest cutoff and if the posterior dist is over 0)
  dat3$fac1both60true[[i]]<-dat3$fac1traits[[i]]%in%dat3$fac1both60[[i]]
  dat3$fac1both60lblsw[[i]]<-dat3$fac1traits[[i]]%in%dat3$fac2both60[[i]]
  dat3$fac2both60true[[i]]<-dat3$fac2traits[[i]]%in%dat3$fac2both60[[i]]
  dat3$fac2both60lblsw[[i]]<-dat3$fac2traits[[i]]%in%dat3$fac1both60[[i]]
  
  #account for traits loading onto multiple factors - L=0.4 & hpdl>0
  dat3$fac1both60overlap[[i]]<-`if`(any(dat3$fac1both60true[[i]]|dat3$fac1both60lblsw[[i]]==TRUE), (dat3$fac1both60true[[i]]&dat3$fac1both60lblsw[[i]]==TRUE), 0)
  dat3$fac2both60overlap[[i]]<-`if`(any(dat3$fac2both60true[[i]]|dat3$fac2both60lblsw[[i]]==TRUE), (dat3$fac2both60true[[i]]&dat3$fac2both60lblsw[[i]]==TRUE), 0)
  
  dat3$fac1both60success[[i]]<-mean(dat3$fac1both60true[[i]])
  dat3$fac1both60wrong[[i]]<-mean(dat3$fac1both60overlap[[i]])
  dat3$fac1both60lblswsuccess[[i]]<-mean(dat3$fac1both60lblsw[[i]])
  dat3$fac2both60success[[i]]<-mean(dat3$fac2both60true[[i]])
  dat3$fac2both60wrong[[i]]<-mean(dat3$fac2both60overlap[[i]])
  dat3$fac2both60lblswsuccess[[i]]<-mean(dat3$fac2both60lblsw[[i]])
}

#make sure the success columns are numeric & replace NaN with 0
x<-c(33:38, 45:50, 57:62, 69:74, 81:86, 93:98, 105:110)
dat3[,x]<-lapply(dat3[,x], as.numeric)
dat3[is.na(dat3)]<-0 #need to do is.na when it's a data frame

#the overall number of things labeled either fac1 or fac2
dat3$fac1overall<-dat3$fac1success + dat3$fac1lblswsuccess - 2*dat3$fac1wrong
dat3$fac2overall<-dat3$fac2success + dat3$fac2lblswsuccess - 2*dat3$fac2wrong

#when L=0.3 (our next lowest cutoff)
dat3$fac1L30overall<-dat3$fac1L30success + dat3$fac1L30lblswsuccess - 2*dat3$fac1L30wrong
dat3$fac2L30overall<-dat3$fac2L30success + dat3$fac2L30lblswsuccess - 2*dat3$fac2L30wrong

#when L=0.3 and hpdl>0 (our next lowest cutoff and if the posterior dist is over 0)
dat3$fac1both30overall<-dat3$fac1both30success + dat3$fac1both30lblswsuccess - 2*dat3$fac1both30wrong
dat3$fac2both30overall<-dat3$fac2both30success + dat3$fac2both30lblswsuccess - 2*dat3$fac2both30wrong

#when L=0.4 (our next lowest cutoff)
dat3$fac1L40overall<-dat3$fac1L40success + dat3$fac1L40lblswsuccess - 2*dat3$fac1L40wrong
dat3$fac2L40overall<-dat3$fac2L40success + dat3$fac2L40lblswsuccess - 2*dat3$fac2L40wrong

#when L=0.4 and hpdl>0 (our next lowest cutoff and if the posterior dist is over 0)
dat3$fac1both40overall<-dat3$fac1both40success + dat3$fac1both40lblswsuccess - 2*dat3$fac1both40wrong
dat3$fac2both40overall<-dat3$fac2both40success + dat3$fac2both40lblswsuccess - 2*dat3$fac2both40wrong

#when L=0.6 (our next lowest cutoff)
dat3$fac1L60overall<-dat3$fac1L60success + dat3$fac1L60lblswsuccess - 2*dat3$fac1L60wrong
dat3$fac2L60overall<-dat3$fac2L60success + dat3$fac2L60lblswsuccess - 2*dat3$fac2L60wrong

#when L=0.6 and hpdl>0 (our next lowest cutoff and if the posterior dist is over 0)
dat3$fac1both60overall<-dat3$fac1both60success + dat3$fac1both60lblswsuccess - 2*dat3$fac1both60wrong
dat3$fac2both60overall<-dat3$fac2both60success + dat3$fac2both60lblswsuccess - 2*dat3$fac2both60wrong

#add successes together for an overall measure between 0-2
dat3$facoverall<-dat3$fac1overall + dat3$fac2overall
dat3$fac30overall<-dat3$fac1L30overall + dat3$fac2L30overall
dat3$facboth30overall<-dat3$fac1both30overall + dat3$fac2both30overall
dat3$fac40overall<-dat3$fac1L40overall + dat3$fac2L40overall
dat3$facboth40overall<-dat3$fac1both40overall + dat3$fac2both40overall
dat3$fac60overall<-dat3$fac1L60overall + dat3$fac2L60overall
dat3$facboth60overall<-dat3$fac1both60overall + dat3$fac2both60overall

#use custom function to scale data columns to between 0-1
minmax <- function(x) {
    return((x- min(x)) /(max(x)-min(x)))}
dat3[,125:131] <- minmax(dat3[,125:131])

#write final data to file
#saveRDS(dat3, file="dat3_shrink.rds") #shrinkage prior
saveRDS(dat3, file="dat3_iid.rds") #iid prior
```

#### Code for simtree502 - done and fixed!
```{r}
library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)

setwd("C:/Users/shelb/OneDrive/Desktop/Chapter 2 Stuff/Simulation/script files/results_filesNov25") #with results files for both priors

#read in data files
mydat<-read.csv("mydat.csv")

mydat$fac.expec.up<-ceiling(mydat$fac.expec)
mydat$fac.expec.down<-floor(mydat$fac.expec)
#this works for making the list of expected stuff!!!
for (i in 1:length(mydat$ntraits)){
  #mydat$alltraits[[i]]<-paste0("V",1:mydat$ntraits[[i]],collapse=',')
  mydat$alltraits[[i]]<-1:mydat$ntraits[[i]]
  mydat$fac1traits[[i]]<-head(mydat$alltraits[[i]],n=(mydat$fac.expec.down[[i]]))
  mydat$fac2traits[[i]]<-tail(mydat$alltraits[[i]],n=(mydat$fac.expec.down[[i]]))
  mydat$fac1traits[[i]]<-paste0("V",mydat$fac1traits[[i]])
  mydat$fac2traits[[i]]<-paste0("V",mydat$fac2traits[[i]])
} 

#now read in the results files - this is for the shrinkage prior
#result<-readRDS("pfa_simtree502_results.rds")
#dat4<-merge(mydat, result, by="iteration")

#now read in the results files - this is for the iid prior
result<-readRDS("pfa_simtree502_resultsiid.rds")
dat4<-merge(mydat, result, by="iteration")

#success measure
for (i in 1:length(dat4$iteration)){
  #fac true <- the proportion of traits we expected in factor 1 vs what we actually got
  #fac label switch <- the traits just got switched between factors
  #sorting just by fac1 and fac2
  dat4$fac1true[[i]]<-dat4$fac1traits[[i]]%in%dat4$fac1[[i]]
  dat4$fac1lblsw[[i]]<-dat4$fac1traits[[i]]%in%dat4$fac2[[i]]
  dat4$fac2true[[i]]<-dat4$fac2traits[[i]]%in%dat4$fac2[[i]]
  dat4$fac2lblsw[[i]]<-dat4$fac2traits[[i]]%in%dat4$fac1[[i]]
  
  #account for traits loading onto multiple factors - fac1 vs fac2
  dat4$fac1overlap[[i]]<-`if`(any(dat4$fac1true[[i]]|dat4$fac1lblsw[[i]]==TRUE),(dat4$fac1true[[i]]&dat4$fac1lblsw[[i]]==TRUE),0)
  dat4$fac2overlap[[i]]<-`if`(any(dat4$fac2true[[i]]|dat4$fac2lblsw[[i]]==TRUE),(dat4$fac2true[[i]]&dat4$fac2lblsw[[i]]==TRUE),0)
  
  dat4$fac1success[[i]]<-mean(dat4$fac1true[[i]])
  dat4$fac1wrong[[i]]<-mean(dat4$fac1overlap[[i]])
  dat4$fac1lblswsuccess[[i]]<-mean(dat4$fac1lblsw[[i]])
  dat4$fac2success[[i]]<-mean(dat4$fac2true[[i]])
  dat4$fac2wrong[[i]]<-mean(dat4$fac2overlap[[i]])
  dat4$fac2lblswsuccess[[i]]<-mean(dat4$fac2lblsw[[i]])
  
  #when L=0.3 (our lowest cutoff)
  dat4$fac1L30true[[i]]<-dat4$fac1traits[[i]]%in%dat4$fac1L30[[i]]
  dat4$fac1L30lblsw[[i]]<-dat4$fac1traits[[i]]%in%dat4$fac2L30[[i]]
  dat4$fac2L30true[[i]]<-dat4$fac2traits[[i]]%in%dat4$fac2L30[[i]]
  dat4$fac2L30lblsw[[i]]<-dat4$fac2traits[[i]]%in%dat4$fac1L30[[i]]
  
  #account for traits loading onto multiple factors - L=0.3
  dat4$fac1L30overlap[[i]]<-`if`(any(dat4$fac1L30true[[i]]|dat4$fac1L30lblsw[[i]]==TRUE), (dat4$fac1L30true[[i]]&dat4$fac1L30lblsw[[i]]==TRUE), 0)
  dat4$fac2L30overlap[[i]]<-`if`(any(dat4$fac2L30true[[i]]|dat4$fac2L30lblsw[[i]]==TRUE), (dat4$fac2L30true[[i]]&dat4$fac2L30lblsw[[i]]==TRUE), 0)
  
  dat4$fac1L30success[[i]]<-mean(dat4$fac1L30true[[i]])
  dat4$fac1L30wrong[[i]]<-mean(dat4$fac1L30overlap[[i]])
  dat4$fac1L30lblswsuccess[[i]]<-mean(dat4$fac1L30lblsw[[i]])
  dat4$fac2L30success[[i]]<-mean(dat4$fac2L30true[[i]])
  dat4$fac2L30wrong[[i]]<-mean(dat4$fac2L30overlap[[i]])
  dat4$fac2L30lblswsuccess[[i]]<-mean(dat4$fac2L30lblsw[[i]])
  
  #when L=0.3 and hpdl>0 (our lowest cutoff and if the posterior dist is over 0)
  dat4$fac1both30true[[i]]<-dat4$fac1traits[[i]]%in%dat4$fac1both30[[i]]
  dat4$fac1both30lblsw[[i]]<-dat4$fac1traits[[i]]%in%dat4$fac2both30[[i]]
  dat4$fac2both30true[[i]]<-dat4$fac2traits[[i]]%in%dat4$fac2both30[[i]]
  dat4$fac2both30lblsw[[i]]<-dat4$fac2traits[[i]]%in%dat4$fac1both30[[i]]
  
  #account for traits loading onto multiple factors - L=0.3 & hpdl>0
  dat4$fac1both30overlap[[i]]<-`if`(any(dat4$fac1both30true[[i]]|dat4$fac1both30lblsw[[i]]==TRUE), (dat4$fac1both30true[[i]]&dat4$fac1both30lblsw[[i]]==TRUE), 0)
  dat4$fac2both30overlap[[i]]<-`if`(any(dat4$fac2both30true[[i]]|dat4$fac2both30lblsw[[i]]==TRUE), (dat4$fac2both30true[[i]]&dat4$fac2both30lblsw[[i]]==TRUE), 0)
  
  dat4$fac1both30success[[i]]<-mean(dat4$fac1both30true[[i]])
  dat4$fac1both30wrong[[i]]<-mean(dat4$fac1both30overlap[[i]])
  dat4$fac1both30lblswsuccess[[i]]<-mean(dat4$fac1both30lblsw[[i]])
  dat4$fac2both30success[[i]]<-mean(dat4$fac2both30true[[i]])
  dat4$fac2both30wrong[[i]]<-mean(dat4$fac2both30overlap[[i]])
  dat4$fac2both30lblswsuccess[[i]]<-mean(dat4$fac2both30lblsw[[i]])
  
  #when L=0.4 (our next lowest cutoff)
  dat4$fac1L40true[[i]]<-dat4$fac1traits[[i]]%in%dat4$fac1L40[[i]]
  dat4$fac1L40lblsw[[i]]<-dat4$fac1traits[[i]]%in%dat4$fac2L40[[i]]
  dat4$fac2L40true[[i]]<-dat4$fac2traits[[i]]%in%dat4$fac2L40[[i]]
  dat4$fac2L40lblsw[[i]]<-dat4$fac2traits[[i]]%in%dat4$fac1L40[[i]]
  
  #account for traits loading onto multiple factors - L=0.4
  dat4$fac1L40overlap[[i]]<-`if`(any(dat4$fac1L40true[[i]]|dat4$fac1L40lblsw[[i]]==TRUE), (dat4$fac1L40true[[i]]&dat4$fac1L40lblsw[[i]]==TRUE), 0)
  dat4$fac2L40overlap[[i]]<-`if`(any(dat4$fac2L40true[[i]]|dat4$fac2L40lblsw[[i]]==TRUE), (dat4$fac2L40true[[i]]&dat4$fac2L40lblsw[[i]]==TRUE), 0)
  
  dat4$fac1L40success[[i]]<-mean(dat4$fac1L40true[[i]])
  dat4$fac1L40wrong[[i]]<-mean(dat4$fac1L40overlap[[i]])
  dat4$fac1L40lblswsuccess[[i]]<-mean(dat4$fac1L40lblsw[[i]])
  dat4$fac2L40success[[i]]<-mean(dat4$fac2L40true[[i]])
  dat4$fac2L40wrong[[i]]<-mean(dat4$fac2L40overlap[[i]])
  dat4$fac2L40lblswsuccess[[i]]<-mean(dat4$fac2L40lblsw[[i]])
  
  #when L=0.4 and hpdl>0 (our next lowest cutoff and if the posterior dist is over 0)
  dat4$fac1both40true[[i]]<-dat4$fac1traits[[i]]%in%dat4$fac1both40[[i]]
  dat4$fac1both40lblsw[[i]]<-dat4$fac1traits[[i]]%in%dat4$fac2both40[[i]]
  dat4$fac2both40true[[i]]<-dat4$fac2traits[[i]]%in%dat4$fac2both40[[i]]
  dat4$fac2both40lblsw[[i]]<-dat4$fac2traits[[i]]%in%dat4$fac1both40[[i]]
  
  #account for traits loading onto multiple factors - L=0.4 & hpdl>0
  dat4$fac1both40overlap[[i]]<-`if`(any(dat4$fac1both40true[[i]]|dat4$fac1both40lblsw[[i]]==TRUE), (dat4$fac1both40true[[i]]&dat4$fac1both40lblsw[[i]]==TRUE), 0)
  dat4$fac2both40overlap[[i]]<-`if`(any(dat4$fac2both40true[[i]]|dat4$fac2both40lblsw[[i]]==TRUE), (dat4$fac2both40true[[i]]&dat4$fac2both40lblsw[[i]]==TRUE), 0)
  
  dat4$fac1both40success[[i]]<-mean(dat4$fac1both40true[[i]])
  dat4$fac1both40wrong[[i]]<-mean(dat4$fac1both40overlap[[i]])
  dat4$fac1both40lblswsuccess[[i]]<-mean(dat4$fac1both40lblsw[[i]])
  dat4$fac2both40success[[i]]<-mean(dat4$fac2both40true[[i]])
  dat4$fac2both40wrong[[i]]<-mean(dat4$fac2both40overlap[[i]])
  dat4$fac2both40lblswsuccess[[i]]<-mean(dat4$fac2both40lblsw[[i]])
  
  #when L=0.6 (our next lowest cutoff)
  dat4$fac1L60true[[i]]<-dat4$fac1traits[[i]]%in%dat4$fac1L60[[i]]
  dat4$fac1L60lblsw[[i]]<-dat4$fac1traits[[i]]%in%dat4$fac2L60[[i]]
  dat4$fac2L60true[[i]]<-dat4$fac2traits[[i]]%in%dat4$fac2L60[[i]]
  dat4$fac2L60lblsw[[i]]<-dat4$fac2traits[[i]]%in%dat4$fac1L60[[i]]
  
  #account for traits loading onto multiple factors - L=0.6
  dat4$fac1L60overlap[[i]]<-`if`(any(dat4$fac1L60true[[i]]|dat4$fac1L60lblsw[[i]]==TRUE), (dat4$fac1L60true[[i]]&dat4$fac1L60lblsw[[i]]==TRUE), 0)
  dat4$fac2L60overlap[[i]]<-`if`(any(dat4$fac2L60true[[i]]|dat4$fac2L60lblsw[[i]]==TRUE), (dat4$fac2L60true[[i]]&dat4$fac2L60lblsw[[i]]==TRUE), 0)
  
  dat4$fac1L60success[[i]]<-mean(dat4$fac1L60true[[i]])
  dat4$fac1L60wrong[[i]]<-mean(dat4$fac1L60overlap[[i]])
  dat4$fac1L60lblswsuccess[[i]]<-mean(dat4$fac1L60lblsw[[i]])
  dat4$fac2L60success[[i]]<-mean(dat4$fac2L60true[[i]])
  dat4$fac2L60wrong[[i]]<-mean(dat4$fac2L60overlap[[i]])
  dat4$fac2L60lblswsuccess[[i]]<-mean(dat4$fac2L60lblsw[[i]])
  
  #when L=0.6 and hpdl>0 (our next lowest cutoff and if the posterior dist is over 0)
  dat4$fac1both60true[[i]]<-dat4$fac1traits[[i]]%in%dat4$fac1both60[[i]]
  dat4$fac1both60lblsw[[i]]<-dat4$fac1traits[[i]]%in%dat4$fac2both60[[i]]
  dat4$fac2both60true[[i]]<-dat4$fac2traits[[i]]%in%dat4$fac2both60[[i]]
  dat4$fac2both60lblsw[[i]]<-dat4$fac2traits[[i]]%in%dat4$fac1both60[[i]]
  
  #account for traits loading onto multiple factors - L=0.4 & hpdl>0
  dat4$fac1both60overlap[[i]]<-`if`(any(dat4$fac1both60true[[i]]|dat4$fac1both60lblsw[[i]]==TRUE), (dat4$fac1both60true[[i]]&dat4$fac1both60lblsw[[i]]==TRUE), 0)
  dat4$fac2both60overlap[[i]]<-`if`(any(dat4$fac2both60true[[i]]|dat4$fac2both60lblsw[[i]]==TRUE), (dat4$fac2both60true[[i]]&dat4$fac2both60lblsw[[i]]==TRUE), 0)
  
  dat4$fac1both60success[[i]]<-mean(dat4$fac1both60true[[i]])
  dat4$fac1both60wrong[[i]]<-mean(dat4$fac1both60overlap[[i]])
  dat4$fac1both60lblswsuccess[[i]]<-mean(dat4$fac1both60lblsw[[i]])
  dat4$fac2both60success[[i]]<-mean(dat4$fac2both60true[[i]])
  dat4$fac2both60wrong[[i]]<-mean(dat4$fac2both60overlap[[i]])
  dat4$fac2both60lblswsuccess[[i]]<-mean(dat4$fac2both60lblsw[[i]])
}

#make sure the success columns are numeric & replace NaN with 0
x<-c(33:38, 45:50, 57:62, 69:74, 81:86, 93:98, 105:110)
dat4[,x]<-lapply(dat4[,x], as.numeric)
dat4[is.na(dat4)]<-0 #need to do is.na when it's a data frame

#the overall number of things labeled either fac1 or fac2
dat4$fac1overall<-dat4$fac1success + dat4$fac1lblswsuccess - 2*dat4$fac1wrong
dat4$fac2overall<-dat4$fac2success + dat4$fac2lblswsuccess - 2*dat4$fac2wrong

#when L=0.3 (our next lowest cutoff)
dat4$fac1L30overall<-dat4$fac1L30success + dat4$fac1L30lblswsuccess - 2*dat4$fac1L30wrong
dat4$fac2L30overall<-dat4$fac2L30success + dat4$fac2L30lblswsuccess - 2*dat4$fac2L30wrong

#when L=0.3 and hpdl>0 (our next lowest cutoff and if the posterior dist is over 0)
dat4$fac1both30overall<-dat4$fac1both30success + dat4$fac1both30lblswsuccess - 2*dat4$fac1both30wrong
dat4$fac2both30overall<-dat4$fac2both30success + dat4$fac2both30lblswsuccess - 2*dat4$fac2both30wrong

#when L=0.4 (our next lowest cutoff)
dat4$fac1L40overall<-dat4$fac1L40success + dat4$fac1L40lblswsuccess - 2*dat4$fac1L40wrong
dat4$fac2L40overall<-dat4$fac2L40success + dat4$fac2L40lblswsuccess - 2*dat4$fac2L40wrong

#when L=0.4 and hpdl>0 (our next lowest cutoff and if the posterior dist is over 0)
dat4$fac1both40overall<-dat4$fac1both40success + dat4$fac1both40lblswsuccess - 2*dat4$fac1both40wrong
dat4$fac2both40overall<-dat4$fac2both40success + dat4$fac2both40lblswsuccess - 2*dat4$fac2both40wrong

#when L=0.6 (our next lowest cutoff)
dat4$fac1L60overall<-dat4$fac1L60success + dat4$fac1L60lblswsuccess - 2*dat4$fac1L60wrong
dat4$fac2L60overall<-dat4$fac2L60success + dat4$fac2L60lblswsuccess - 2*dat4$fac2L60wrong

#when L=0.6 and hpdl>0 (our next lowest cutoff and if the posterior dist is over 0)
dat4$fac1both60overall<-dat4$fac1both60success + dat4$fac1both60lblswsuccess - 2*dat4$fac1both60wrong
dat4$fac2both60overall<-dat4$fac2both60success + dat4$fac2both60lblswsuccess - 2*dat4$fac2both60wrong

#add successes together for an overall measure between 0-2
dat4$facoverall<-dat4$fac1overall + dat4$fac2overall
dat4$fac30overall<-dat4$fac1L30overall + dat4$fac2L30overall
dat4$facboth30overall<-dat4$fac1both30overall + dat4$fac2both30overall
dat4$fac40overall<-dat4$fac1L40overall + dat4$fac2L40overall
dat4$facboth40overall<-dat4$fac1both40overall + dat4$fac2both40overall
dat4$fac60overall<-dat4$fac1L60overall + dat4$fac2L60overall
dat4$facboth60overall<-dat4$fac1both60overall + dat4$fac2both60overall

#use custom function to scale data columns to between 0-1
minmax <- function(x) {
    return((x- min(x)) /(max(x)-min(x)))}
dat4[,125:131] <- minmax(dat4[,125:131])

#write final data to file
#saveRDS(dat4, file="dat4_shrink.rds") #shrinkage prior
saveRDS(dat4, file="dat4_iid.rds") #iid prior
```

#### Code for error_sim - done and fixed!
```{r}
library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)

setwd("C:/Users/shelb/OneDrive/Desktop/Chapter 2 Stuff/Simulation/script files/results_filesNov25")

#read in data files
mydat<-read.csv("mydat_error.csv")
#this works for making the list of expected stuff!!!
for (i in 1:length(mydat$ntraits)){
  #mydat$alltraits[[i]]<-paste0("V",1:mydat$ntraits[[i]],collapse=',')
  mydat$alltraits[[i]]<-1:mydat$ntraits[[i]]
  mydat$fac1traits[[i]]<-head(mydat$alltraits[[i]],n=(mydat$fac.expec[[i]])) 
  mydat$fac2traits[[i]]<-tail(mydat$alltraits[[i]],n=(mydat$fac.expec[[i]])) 
  mydat$fac1traits[[i]]<-paste0("V",mydat$fac1traits[[i]])
  mydat$fac2traits[[i]]<-paste0("V",mydat$fac2traits[[i]])
} 

#outdated code for expected traits
{#expected traits for 20-tip analyses
#mydat$fac1traits[[1]]<-c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10") #scenario 1 = second half of traits have higher variance
#mydat$fac2traits[[1]]<-c("V11","V12","V13","V14","V15","V16","V17","V18","V19","V20")

#mydat$fac1traits[[2]]<-c("V1","V2","V3","V4","V5","V6","V7","V8") #scenario 2 = last 2 traits of each half have even higher variance
#mydat$fac2traits[[2]]<-c("V11","V12","V13","V14","V15","V16","V17","V18")

#mydat$fac1traits[[3]]<-c("V1","V2","V3","V4","V5","V6","V7","V8") #scenario 3 = last 2 traits of each half have even higher variance
#mydat$fac2traits[[3]]<-c("V11","V12","V13","V14","V15","V16","V17","V18")

#mydat$fac1traits[[4]]<-c("V1","V2","V3","V4","V5","V6","V7","V8") #scenario 4 = last 2 traits of each half have random numbers
#mydat$fac2traits[[4]]<-c("V11","V12","V13","V14","V15","V16","V17","V18")



#expected traits for 40-tip analysis
#mydat$fac1traits[[5]]<-c("V1","V2","V3","V4","V5","V6","V7","V9","V10","V11","V12","V14","V15","V16","V17","V19","V20") #except for 8,13,18
#mydat$fac2traits[[5]]<-c("V21","V22","V23","V24","V25","V26","V27","V29","V30","V31","V32","V34","V35","V36","V37","V39","V40") #except for 28,33,38
}

#updated code for expected traits
result<-readRDS("error_sim_success.rds")
dat<-merge(mydat, result, by="iteration")

#success measure
for (i in 1:length(dat$iteration)){
  #fac true <- the proportion of traits we expected in factor 1 vs what we actually got
  #fac label switch <- the traits just got switched between factors
  #sorting just by fac1 and fac2
  dat$fac1true[[i]]<-dat$fac1traits[[i]]%in%dat$fac1[[i]]
  dat$fac1lblsw[[i]]<-dat$fac1traits[[i]]%in%dat$fac2[[i]]
  dat$fac2true[[i]]<-dat$fac2traits[[i]]%in%dat$fac2[[i]]
  dat$fac2lblsw[[i]]<-dat$fac2traits[[i]]%in%dat$fac1[[i]]
  
  #account for traits loading onto multiple factors - fac1 vs fac2
  dat$fac1overlap[[i]]<-`if`(any(dat$fac1true[[i]]|dat$fac1lblsw[[i]]==TRUE),(dat$fac1true[[i]]&dat$fac1lblsw[[i]]==TRUE),0)
  dat$fac2overlap[[i]]<-`if`(any(dat$fac2true[[i]]|dat$fac2lblsw[[i]]==TRUE),(dat$fac2true[[i]]&dat$fac2lblsw[[i]]==TRUE),0)
  
  dat$fac1success[[i]]<-mean(dat$fac1true[[i]])
  dat$fac1wrong[[i]]<-mean(dat$fac1overlap[[i]])
  dat$fac1lblswsuccess[[i]]<-mean(dat$fac1lblsw[[i]])
  dat$fac2success[[i]]<-mean(dat$fac2true[[i]])
  dat$fac2wrong[[i]]<-mean(dat$fac2overlap[[i]])
  dat$fac2lblswsuccess[[i]]<-mean(dat$fac2lblsw[[i]])
  
  #when L=0.3 (our lowest cutoff)
  dat$fac1L30true[[i]]<-dat$fac1traits[[i]]%in%dat$fac1L30[[i]]
  dat$fac1L30lblsw[[i]]<-dat$fac1traits[[i]]%in%dat$fac2L30[[i]]
  dat$fac2L30true[[i]]<-dat$fac2traits[[i]]%in%dat$fac2L30[[i]]
  dat$fac2L30lblsw[[i]]<-dat$fac2traits[[i]]%in%dat$fac1L30[[i]]
  
  #account for traits loading onto multiple factors - L=0.3
  dat$fac1L30overlap[[i]]<-`if`(any(dat$fac1L30true[[i]]|dat$fac1L30lblsw[[i]]==TRUE), (dat$fac1L30true[[i]]&dat$fac1L30lblsw[[i]]==TRUE), 0)
  dat$fac2L30overlap[[i]]<-`if`(any(dat$fac2L30true[[i]]|dat$fac2L30lblsw[[i]]==TRUE), (dat$fac2L30true[[i]]&dat$fac2L30lblsw[[i]]==TRUE), 0)
  
  dat$fac1L30success[[i]]<-mean(dat$fac1L30true[[i]])
  dat$fac1L30wrong[[i]]<-mean(dat$fac1L30overlap[[i]])
  dat$fac1L30lblswsuccess[[i]]<-mean(dat$fac1L30lblsw[[i]])
  dat$fac2L30success[[i]]<-mean(dat$fac2L30true[[i]])
  dat$fac2L30wrong[[i]]<-mean(dat$fac2L30overlap[[i]])
  dat$fac2L30lblswsuccess[[i]]<-mean(dat$fac2L30lblsw[[i]])
  
  #when L=0.3 and hpdl>0 (our lowest cutoff and if the posterior dist is over 0)
  dat$fac1both30true[[i]]<-dat$fac1traits[[i]]%in%dat$fac1both30[[i]]
  dat$fac1both30lblsw[[i]]<-dat$fac1traits[[i]]%in%dat$fac2both30[[i]]
  dat$fac2both30true[[i]]<-dat$fac2traits[[i]]%in%dat$fac2both30[[i]]
  dat$fac2both30lblsw[[i]]<-dat$fac2traits[[i]]%in%dat$fac1both30[[i]]
  
  #account for traits loading onto multiple factors - L=0.3 & hpdl>0
  dat$fac1both30overlap[[i]]<-`if`(any(dat$fac1both30true[[i]]|dat$fac1both30lblsw[[i]]==TRUE), (dat$fac1both30true[[i]]&dat$fac1both30lblsw[[i]]==TRUE), 0)
  dat$fac2both30overlap[[i]]<-`if`(any(dat$fac2both30true[[i]]|dat$fac2both30lblsw[[i]]==TRUE), (dat$fac2both30true[[i]]&dat$fac2both30lblsw[[i]]==TRUE), 0)
  
  dat$fac1both30success[[i]]<-mean(dat$fac1both30true[[i]])
  dat$fac1both30wrong[[i]]<-mean(dat$fac1both30overlap[[i]])
  dat$fac1both30lblswsuccess[[i]]<-mean(dat$fac1both30lblsw[[i]])
  dat$fac2both30success[[i]]<-mean(dat$fac2both30true[[i]])
  dat$fac2both30wrong[[i]]<-mean(dat$fac2both30overlap[[i]])
  dat$fac2both30lblswsuccess[[i]]<-mean(dat$fac2both30lblsw[[i]])
  
  #when L=0.4 (our next lowest cutoff)
  dat$fac1L40true[[i]]<-dat$fac1traits[[i]]%in%dat$fac1L40[[i]]
  dat$fac1L40lblsw[[i]]<-dat$fac1traits[[i]]%in%dat$fac2L40[[i]]
  dat$fac2L40true[[i]]<-dat$fac2traits[[i]]%in%dat$fac2L40[[i]]
  dat$fac2L40lblsw[[i]]<-dat$fac2traits[[i]]%in%dat$fac1L40[[i]]
  
  #account for traits loading onto multiple factors - L=0.4
  dat$fac1L40overlap[[i]]<-`if`(any(dat$fac1L40true[[i]]|dat$fac1L40lblsw[[i]]==TRUE), (dat$fac1L40true[[i]]&dat$fac1L40lblsw[[i]]==TRUE), 0)
  dat$fac2L40overlap[[i]]<-`if`(any(dat$fac2L40true[[i]]|dat$fac2L40lblsw[[i]]==TRUE), (dat$fac2L40true[[i]]&dat$fac2L40lblsw[[i]]==TRUE), 0)
  
  dat$fac1L40success[[i]]<-mean(dat$fac1L40true[[i]])
  dat$fac1L40wrong[[i]]<-mean(dat$fac1L40overlap[[i]])
  dat$fac1L40lblswsuccess[[i]]<-mean(dat$fac1L40lblsw[[i]])
  dat$fac2L40success[[i]]<-mean(dat$fac2L40true[[i]])
  dat$fac2L40wrong[[i]]<-mean(dat$fac2L40overlap[[i]])
  dat$fac2L40lblswsuccess[[i]]<-mean(dat$fac2L40lblsw[[i]])
  
  #when L=0.4 and hpdl>0 (our next lowest cutoff and if the posterior dist is over 0)
  dat$fac1both40true[[i]]<-dat$fac1traits[[i]]%in%dat$fac1both40[[i]]
  dat$fac1both40lblsw[[i]]<-dat$fac1traits[[i]]%in%dat$fac2both40[[i]]
  dat$fac2both40true[[i]]<-dat$fac2traits[[i]]%in%dat$fac2both40[[i]]
  dat$fac2both40lblsw[[i]]<-dat$fac2traits[[i]]%in%dat$fac1both40[[i]]
  
  #account for traits loading onto multiple factors - L=0.4 & hpdl>0
  dat$fac1both40overlap[[i]]<-`if`(any(dat$fac1both40true[[i]]|dat$fac1both40lblsw[[i]]==TRUE), (dat$fac1both40true[[i]]&dat$fac1both40lblsw[[i]]==TRUE), 0)
  dat$fac2both40overlap[[i]]<-`if`(any(dat$fac2both40true[[i]]|dat$fac2both40lblsw[[i]]==TRUE), (dat$fac2both40true[[i]]&dat$fac2both40lblsw[[i]]==TRUE), 0)
  
  dat$fac1both40success[[i]]<-mean(dat$fac1both40true[[i]])
  dat$fac1both40wrong[[i]]<-mean(dat$fac1both40overlap[[i]])
  dat$fac1both40lblswsuccess[[i]]<-mean(dat$fac1both40lblsw[[i]])
  dat$fac2both40success[[i]]<-mean(dat$fac2both40true[[i]])
  dat$fac2both40wrong[[i]]<-mean(dat$fac2both40overlap[[i]])
  dat$fac2both40lblswsuccess[[i]]<-mean(dat$fac2both40lblsw[[i]])
  
  #when L=0.6 (our next lowest cutoff)
  dat$fac1L60true[[i]]<-dat$fac1traits[[i]]%in%dat$fac1L60[[i]]
  dat$fac1L60lblsw[[i]]<-dat$fac1traits[[i]]%in%dat$fac2L60[[i]]
  dat$fac2L60true[[i]]<-dat$fac2traits[[i]]%in%dat$fac2L60[[i]]
  dat$fac2L60lblsw[[i]]<-dat$fac2traits[[i]]%in%dat$fac1L60[[i]]
  
  #account for traits loading onto multiple factors - L=0.6
  dat$fac1L60overlap[[i]]<-`if`(any(dat$fac1L60true[[i]]|dat$fac1L60lblsw[[i]]==TRUE), (dat$fac1L60true[[i]]&dat$fac1L60lblsw[[i]]==TRUE), 0)
  dat$fac2L60overlap[[i]]<-`if`(any(dat$fac2L60true[[i]]|dat$fac2L60lblsw[[i]]==TRUE), (dat$fac2L60true[[i]]&dat$fac2L60lblsw[[i]]==TRUE), 0)
  
  dat$fac1L60success[[i]]<-mean(dat$fac1L60true[[i]])
  dat$fac1L60wrong[[i]]<-mean(dat$fac1L60overlap[[i]])
  dat$fac1L60lblswsuccess[[i]]<-mean(dat$fac1L60lblsw[[i]])
  dat$fac2L60success[[i]]<-mean(dat$fac2L60true[[i]])
  dat$fac2L60wrong[[i]]<-mean(dat$fac2L60overlap[[i]])
  dat$fac2L60lblswsuccess[[i]]<-mean(dat$fac2L60lblsw[[i]])
  
  #when L=0.6 and hpdl>0 (our next lowest cutoff and if the posterior dist is over 0)
  dat$fac1both60true[[i]]<-dat$fac1traits[[i]]%in%dat$fac1both60[[i]]
  dat$fac1both60lblsw[[i]]<-dat$fac1traits[[i]]%in%dat$fac2both60[[i]]
  dat$fac2both60true[[i]]<-dat$fac2traits[[i]]%in%dat$fac2both60[[i]]
  dat$fac2both60lblsw[[i]]<-dat$fac2traits[[i]]%in%dat$fac1both60[[i]]
  
  #account for traits loading onto multiple factors - L=0.4 & hpdl>0
  dat$fac1both60overlap[[i]]<-`if`(any(dat$fac1both60true[[i]]|dat$fac1both60lblsw[[i]]==TRUE), (dat$fac1both60true[[i]]&dat$fac1both60lblsw[[i]]==TRUE), 0)
  dat$fac2both60overlap[[i]]<-`if`(any(dat$fac2both60true[[i]]|dat$fac2both60lblsw[[i]]==TRUE), (dat$fac2both60true[[i]]&dat$fac2both60lblsw[[i]]==TRUE), 0)
  
  dat$fac1both60success[[i]]<-mean(dat$fac1both60true[[i]])
  dat$fac1both60wrong[[i]]<-mean(dat$fac1both60overlap[[i]])
  dat$fac1both60lblswsuccess[[i]]<-mean(dat$fac1both60lblsw[[i]])
  dat$fac2both60success[[i]]<-mean(dat$fac2both60true[[i]])
  dat$fac2both60wrong[[i]]<-mean(dat$fac2both60overlap[[i]])
  dat$fac2both60lblswsuccess[[i]]<-mean(dat$fac2both60lblsw[[i]])
}

#make sure the success columns are numeric & replace NaN with 0
x<-c(31:36, 43:48, 55:60, 67:72, 79:84, 91:96, 103:108)
dat[,x]<-lapply(dat[,x], as.numeric)
dat[is.na(dat)]<-0 #need to do is.na when it's a data frame

#the overall number of things labeled either fac1 or fac2
dat$fac1overall<-dat$fac1success + dat$fac1lblswsuccess - 2*dat$fac1wrong
dat$fac2overall<-dat$fac2success + dat$fac2lblswsuccess - 2*dat$fac2wrong

#when L=0.3 (our next lowest cutoff)
dat$fac1L30overall<-dat$fac1L30success + dat$fac1L30lblswsuccess - 2*dat$fac1L30wrong
dat$fac2L30overall<-dat$fac2L30success + dat$fac2L30lblswsuccess - 2*dat$fac2L30wrong

#when L=0.3 and hpdl>0 (our next lowest cutoff and if the posterior dist is over 0)
dat$fac1both30overall<-dat$fac1both30success + dat$fac1both30lblswsuccess - 2*dat$fac1both30wrong
dat$fac2both30overall<-dat$fac2both30success + dat$fac2both30lblswsuccess - 2*dat$fac2both30wrong

#when L=0.4 (our next lowest cutoff)
dat$fac1L40overall<-dat$fac1L40success + dat$fac1L40lblswsuccess - 2*dat$fac1L40wrong
dat$fac2L40overall<-dat$fac2L40success + dat$fac2L40lblswsuccess - 2*dat$fac2L40wrong

#when L=0.4 and hpdl>0 (our next lowest cutoff and if the posterior dist is over 0)
dat$fac1both40overall<-dat$fac1both40success + dat$fac1both40lblswsuccess - 2*dat$fac1both40wrong
dat$fac2both40overall<-dat$fac2both40success + dat$fac2both40lblswsuccess - 2*dat$fac2both40wrong

#when L=0.6 (our next lowest cutoff)
dat$fac1L60overall<-dat$fac1L60success + dat$fac1L60lblswsuccess - 2*dat$fac1L60wrong
dat$fac2L60overall<-dat$fac2L60success + dat$fac2L60lblswsuccess - 2*dat$fac2L60wrong

#when L=0.6 and hpdl>0 (our next lowest cutoff and if the posterior dist is over 0)
dat$fac1both60overall<-dat$fac1both60success + dat$fac1both60lblswsuccess - 2*dat$fac1both60wrong
dat$fac2both60overall<-dat$fac2both60success + dat$fac2both60lblswsuccess - 2*dat$fac2both60wrong

#add successes together for an overall measure between 0-2
dat$facoverall<-dat$fac1overall + dat$fac2overall
dat$fac30overall<-dat$fac1L30overall + dat$fac2L30overall
dat$facboth30overall<-dat$fac1both30overall + dat$fac2both30overall
dat$fac40overall<-dat$fac1L40overall + dat$fac2L40overall
dat$facboth40overall<-dat$fac1both40overall + dat$fac2both40overall
dat$fac60overall<-dat$fac1L60overall + dat$fac2L60overall
dat$facboth60overall<-dat$fac1both60overall + dat$fac2both60overall
dat$minscale<-as.numeric('0')
dat$maxscale<-as.numeric('2')

#use custom function to scale data columns to between 0-1
minmax <- function(x) {
    return((x- min(x)) /(max(x)-min(x)))}
dat[,123:131] <- minmax(dat[,123:131])

#write final data to file
saveRDS(dat, file="dat_error.rds")
```

#### Code for discrete_sim - done and fixed!
```{r}
library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)

setwd("C:/Users/shelb/OneDrive/Desktop/Chapter 2 Stuff/Simulation/script files/results_filesNov25")

#read in data files
mydat20<-read.csv("mydat_discrete20.csv")
mydat40<-read.csv("mydat_discrete40.csv")
#make list of expected traits
for (i in 1:length(mydat20$ntraits)){
  mydat20$alltraits[[i]]<-1:mydat20$ntraits[[i]]
  mydat20$fac1traits[[i]]<-head(mydat20$alltraits[[i]],n=(mydat20$fac.expec[[i]])) 
  mydat20$fac2traits[[i]]<-tail(mydat20$alltraits[[i]],n=(mydat20$fac.expec[[i]]))
  mydat20$fac1traits[[i]]<-paste0("V",mydat20$fac1traits[[i]])
  mydat20$fac2traits[[i]]<-paste0("V",mydat20$fac2traits[[i]])
} 

for (i in 1:length(mydat40$ntraits)){
  mydat40$alltraits[[i]]<-1:mydat40$ntraits[[i]]
  mydat40$fac1traits[[i]]<-head(mydat40$alltraits[[i]],n=(mydat40$fac.expec[[i]])) 
  mydat40$fac2traits[[i]]<-tail(mydat40$alltraits[[i]],n=(mydat40$fac.expec[[i]]))
  mydat40$fac1traits[[i]]<-paste0("V",mydat40$fac1traits[[i]])
  mydat40$fac2traits[[i]]<-paste0("V",mydat40$fac2traits[[i]])
}

#now read in the results files
#result<-readRDS("discrete20_results.rds")
result<-readRDS("discrete40_results.rds")

dat<-merge(mydat20, result, by="iteration")
#dat<-merge(mydat40, result, by="iteration")

#success measure
for (i in 1:length(dat$iteration)){
  #fac true <- the proportion of traits we expected in factor 1 vs what we actually got
  #fac label switch <- the traits just got switched between factors
  #sorting just by fac1 and fac2
  dat$fac1true[[i]]<-dat$fac1traits[[i]]%in%dat$fac1[[i]]
  dat$fac1lblsw[[i]]<-dat$fac1traits[[i]]%in%dat$fac2[[i]]
  dat$fac2true[[i]]<-dat$fac2traits[[i]]%in%dat$fac2[[i]]
  dat$fac2lblsw[[i]]<-dat$fac2traits[[i]]%in%dat$fac1[[i]]
  
  #account for traits loading onto multiple factors - fac1 vs fac2
  dat$fac1overlap[[i]]<-`if`(any(dat$fac1true[[i]]|dat$fac1lblsw[[i]]==TRUE),(dat$fac1true[[i]]&dat$fac1lblsw[[i]]==TRUE),0)
  dat$fac2overlap[[i]]<-`if`(any(dat$fac2true[[i]]|dat$fac2lblsw[[i]]==TRUE),(dat$fac2true[[i]]&dat$fac2lblsw[[i]]==TRUE),0)
  
  dat$fac1success[[i]]<-mean(dat$fac1true[[i]])
  dat$fac1wrong[[i]]<-mean(dat$fac1overlap[[i]])
  dat$fac1lblswsuccess[[i]]<-mean(dat$fac1lblsw[[i]])
  dat$fac2success[[i]]<-mean(dat$fac2true[[i]])
  dat$fac2wrong[[i]]<-mean(dat$fac2overlap[[i]])
  dat$fac2lblswsuccess[[i]]<-mean(dat$fac2lblsw[[i]])
  
  #when L=0.3 (our lowest cutoff)
  dat$fac1L30true[[i]]<-dat$fac1traits[[i]]%in%dat$fac1L30[[i]]
  dat$fac1L30lblsw[[i]]<-dat$fac1traits[[i]]%in%dat$fac2L30[[i]]
  dat$fac2L30true[[i]]<-dat$fac2traits[[i]]%in%dat$fac2L30[[i]]
  dat$fac2L30lblsw[[i]]<-dat$fac2traits[[i]]%in%dat$fac1L30[[i]]
  
  #account for traits loading onto multiple factors - L=0.3
  dat$fac1L30overlap[[i]]<-`if`(any(dat$fac1L30true[[i]]|dat$fac1L30lblsw[[i]]==TRUE), (dat$fac1L30true[[i]]&dat$fac1L30lblsw[[i]]==TRUE), 0)
  dat$fac2L30overlap[[i]]<-`if`(any(dat$fac2L30true[[i]]|dat$fac2L30lblsw[[i]]==TRUE), (dat$fac2L30true[[i]]&dat$fac2L30lblsw[[i]]==TRUE), 0)
  
  dat$fac1L30success[[i]]<-mean(dat$fac1L30true[[i]])
  dat$fac1L30wrong[[i]]<-mean(dat$fac1L30overlap[[i]])
  dat$fac1L30lblswsuccess[[i]]<-mean(dat$fac1L30lblsw[[i]])
  dat$fac2L30success[[i]]<-mean(dat$fac2L30true[[i]])
  dat$fac2L30wrong[[i]]<-mean(dat$fac2L30overlap[[i]])
  dat$fac2L30lblswsuccess[[i]]<-mean(dat$fac2L30lblsw[[i]])
  
  #when L=0.3 and hpdl>0 (our lowest cutoff and if the posterior dist is over 0)
  dat$fac1both30true[[i]]<-dat$fac1traits[[i]]%in%dat$fac1both30[[i]]
  dat$fac1both30lblsw[[i]]<-dat$fac1traits[[i]]%in%dat$fac2both30[[i]]
  dat$fac2both30true[[i]]<-dat$fac2traits[[i]]%in%dat$fac2both30[[i]]
  dat$fac2both30lblsw[[i]]<-dat$fac2traits[[i]]%in%dat$fac1both30[[i]]
  
  #account for traits loading onto multiple factors - L=0.3 & hpdl>0
  dat$fac1both30overlap[[i]]<-`if`(any(dat$fac1both30true[[i]]|dat$fac1both30lblsw[[i]]==TRUE), (dat$fac1both30true[[i]]&dat$fac1both30lblsw[[i]]==TRUE), 0)
  dat$fac2both30overlap[[i]]<-`if`(any(dat$fac2both30true[[i]]|dat$fac2both30lblsw[[i]]==TRUE), (dat$fac2both30true[[i]]&dat$fac2both30lblsw[[i]]==TRUE), 0)
  
  dat$fac1both30success[[i]]<-mean(dat$fac1both30true[[i]])
  dat$fac1both30wrong[[i]]<-mean(dat$fac1both30overlap[[i]])
  dat$fac1both30lblswsuccess[[i]]<-mean(dat$fac1both30lblsw[[i]])
  dat$fac2both30success[[i]]<-mean(dat$fac2both30true[[i]])
  dat$fac2both30wrong[[i]]<-mean(dat$fac2both30overlap[[i]])
  dat$fac2both30lblswsuccess[[i]]<-mean(dat$fac2both30lblsw[[i]])
  
  #when L=0.4 (our next lowest cutoff)
  dat$fac1L40true[[i]]<-dat$fac1traits[[i]]%in%dat$fac1L40[[i]]
  dat$fac1L40lblsw[[i]]<-dat$fac1traits[[i]]%in%dat$fac2L40[[i]]
  dat$fac2L40true[[i]]<-dat$fac2traits[[i]]%in%dat$fac2L40[[i]]
  dat$fac2L40lblsw[[i]]<-dat$fac2traits[[i]]%in%dat$fac1L40[[i]]
  
  #account for traits loading onto multiple factors - L=0.4
  dat$fac1L40overlap[[i]]<-`if`(any(dat$fac1L40true[[i]]|dat$fac1L40lblsw[[i]]==TRUE), (dat$fac1L40true[[i]]&dat$fac1L40lblsw[[i]]==TRUE), 0)
  dat$fac2L40overlap[[i]]<-`if`(any(dat$fac2L40true[[i]]|dat$fac2L40lblsw[[i]]==TRUE), (dat$fac2L40true[[i]]&dat$fac2L40lblsw[[i]]==TRUE), 0)
  
  dat$fac1L40success[[i]]<-mean(dat$fac1L40true[[i]])
  dat$fac1L40wrong[[i]]<-mean(dat$fac1L40overlap[[i]])
  dat$fac1L40lblswsuccess[[i]]<-mean(dat$fac1L40lblsw[[i]])
  dat$fac2L40success[[i]]<-mean(dat$fac2L40true[[i]])
  dat$fac2L40wrong[[i]]<-mean(dat$fac2L40overlap[[i]])
  dat$fac2L40lblswsuccess[[i]]<-mean(dat$fac2L40lblsw[[i]])
  
  #when L=0.4 and hpdl>0 (our next lowest cutoff and if the posterior dist is over 0)
  dat$fac1both40true[[i]]<-dat$fac1traits[[i]]%in%dat$fac1both40[[i]]
  dat$fac1both40lblsw[[i]]<-dat$fac1traits[[i]]%in%dat$fac2both40[[i]]
  dat$fac2both40true[[i]]<-dat$fac2traits[[i]]%in%dat$fac2both40[[i]]
  dat$fac2both40lblsw[[i]]<-dat$fac2traits[[i]]%in%dat$fac1both40[[i]]
  
  #account for traits loading onto multiple factors - L=0.4 & hpdl>0
  dat$fac1both40overlap[[i]]<-`if`(any(dat$fac1both40true[[i]]|dat$fac1both40lblsw[[i]]==TRUE), (dat$fac1both40true[[i]]&dat$fac1both40lblsw[[i]]==TRUE), 0)
  dat$fac2both40overlap[[i]]<-`if`(any(dat$fac2both40true[[i]]|dat$fac2both40lblsw[[i]]==TRUE), (dat$fac2both40true[[i]]&dat$fac2both40lblsw[[i]]==TRUE), 0)
  
  dat$fac1both40success[[i]]<-mean(dat$fac1both40true[[i]])
  dat$fac1both40wrong[[i]]<-mean(dat$fac1both40overlap[[i]])
  dat$fac1both40lblswsuccess[[i]]<-mean(dat$fac1both40lblsw[[i]])
  dat$fac2both40success[[i]]<-mean(dat$fac2both40true[[i]])
  dat$fac2both40wrong[[i]]<-mean(dat$fac2both40overlap[[i]])
  dat$fac2both40lblswsuccess[[i]]<-mean(dat$fac2both40lblsw[[i]])
  
  #when L=0.6 (our next lowest cutoff)
  dat$fac1L60true[[i]]<-dat$fac1traits[[i]]%in%dat$fac1L60[[i]]
  dat$fac1L60lblsw[[i]]<-dat$fac1traits[[i]]%in%dat$fac2L60[[i]]
  dat$fac2L60true[[i]]<-dat$fac2traits[[i]]%in%dat$fac2L60[[i]]
  dat$fac2L60lblsw[[i]]<-dat$fac2traits[[i]]%in%dat$fac1L60[[i]]
  
  #account for traits loading onto multiple factors - L=0.6
  dat$fac1L60overlap[[i]]<-`if`(any(dat$fac1L60true[[i]]|dat$fac1L60lblsw[[i]]==TRUE), (dat$fac1L60true[[i]]&dat$fac1L60lblsw[[i]]==TRUE), 0)
  dat$fac2L60overlap[[i]]<-`if`(any(dat$fac2L60true[[i]]|dat$fac2L60lblsw[[i]]==TRUE), (dat$fac2L60true[[i]]&dat$fac2L60lblsw[[i]]==TRUE), 0)
  
  dat$fac1L60success[[i]]<-mean(dat$fac1L60true[[i]])
  dat$fac1L60wrong[[i]]<-mean(dat$fac1L60overlap[[i]])
  dat$fac1L60lblswsuccess[[i]]<-mean(dat$fac1L60lblsw[[i]])
  dat$fac2L60success[[i]]<-mean(dat$fac2L60true[[i]])
  dat$fac2L60wrong[[i]]<-mean(dat$fac2L60overlap[[i]])
  dat$fac2L60lblswsuccess[[i]]<-mean(dat$fac2L60lblsw[[i]])
  
  #when L=0.6 and hpdl>0 (our next lowest cutoff and if the posterior dist is over 0)
  dat$fac1both60true[[i]]<-dat$fac1traits[[i]]%in%dat$fac1both60[[i]]
  dat$fac1both60lblsw[[i]]<-dat$fac1traits[[i]]%in%dat$fac2both60[[i]]
  dat$fac2both60true[[i]]<-dat$fac2traits[[i]]%in%dat$fac2both60[[i]]
  dat$fac2both60lblsw[[i]]<-dat$fac2traits[[i]]%in%dat$fac1both60[[i]]
  
  #account for traits loading onto multiple factors - L=0.4 & hpdl>0
  dat$fac1both60overlap[[i]]<-`if`(any(dat$fac1both60true[[i]]|dat$fac1both60lblsw[[i]]==TRUE), (dat$fac1both60true[[i]]&dat$fac1both60lblsw[[i]]==TRUE), 0)
  dat$fac2both60overlap[[i]]<-`if`(any(dat$fac2both60true[[i]]|dat$fac2both60lblsw[[i]]==TRUE), (dat$fac2both60true[[i]]&dat$fac2both60lblsw[[i]]==TRUE), 0)
  
  dat$fac1both60success[[i]]<-mean(dat$fac1both60true[[i]])
  dat$fac1both60wrong[[i]]<-mean(dat$fac1both60overlap[[i]])
  dat$fac1both60lblswsuccess[[i]]<-mean(dat$fac1both60lblsw[[i]])
  dat$fac2both60success[[i]]<-mean(dat$fac2both60true[[i]])
  dat$fac2both60wrong[[i]]<-mean(dat$fac2both60overlap[[i]])
  dat$fac2both60lblswsuccess[[i]]<-mean(dat$fac2both60lblsw[[i]])
}

#make sure the success columns are numeric & replace NaN with 0
x<-c(32:37, 44:49, 56:61, 68:73, 80:85, 92:97, 104:109)
dat[,x]<-lapply(dat[,x], as.numeric)
dat[is.na(dat)]<-0 #need to do is.na when it's a data frame

#the overall number of things labeled either fac1 or fac2
dat$fac1overall<-dat$fac1success + dat$fac1lblswsuccess - 2*dat$fac1wrong
dat$fac2overall<-dat$fac2success + dat$fac2lblswsuccess - 2*dat$fac2wrong

#when L=0.3 (our next lowest cutoff)
dat$fac1L30overall<-dat$fac1L30success + dat$fac1L30lblswsuccess - 2*dat$fac1L30wrong
dat$fac2L30overall<-dat$fac2L30success + dat$fac2L30lblswsuccess - 2*dat$fac2L30wrong

#when L=0.3 and hpdl>0 (our next lowest cutoff and if the posterior dist is over 0)
dat$fac1both30overall<-dat$fac1both30success + dat$fac1both30lblswsuccess - 2*dat$fac1both30wrong
dat$fac2both30overall<-dat$fac2both30success + dat$fac2both30lblswsuccess - 2*dat$fac2both30wrong

#when L=0.4 (our next lowest cutoff)
dat$fac1L40overall<-dat$fac1L40success + dat$fac1L40lblswsuccess - 2*dat$fac1L40wrong
dat$fac2L40overall<-dat$fac2L40success + dat$fac2L40lblswsuccess - 2*dat$fac2L40wrong

#when L=0.4 and hpdl>0 (our next lowest cutoff and if the posterior dist is over 0)
dat$fac1both40overall<-dat$fac1both40success + dat$fac1both40lblswsuccess - 2*dat$fac1both40wrong
dat$fac2both40overall<-dat$fac2both40success + dat$fac2both40lblswsuccess - 2*dat$fac2both40wrong

#when L=0.6 (our next lowest cutoff)
dat$fac1L60overall<-dat$fac1L60success + dat$fac1L60lblswsuccess - 2*dat$fac1L60wrong
dat$fac2L60overall<-dat$fac2L60success + dat$fac2L60lblswsuccess - 2*dat$fac2L60wrong

#when L=0.6 and hpdl>0 (our next lowest cutoff and if the posterior dist is over 0)
dat$fac1both60overall<-dat$fac1both60success + dat$fac1both60lblswsuccess - 2*dat$fac1both60wrong
dat$fac2both60overall<-dat$fac2both60success + dat$fac2both60lblswsuccess - 2*dat$fac2both60wrong

#add successes together for an overall measure between 0-2
dat$facoverall<-dat$fac1overall + dat$fac2overall
dat$fac30overall<-dat$fac1L30overall + dat$fac2L30overall
dat$facboth30overall<-dat$fac1both30overall + dat$fac2both30overall
dat$fac40overall<-dat$fac1L40overall + dat$fac2L40overall
dat$facboth40overall<-dat$fac1both40overall + dat$fac2both40overall
dat$fac60overall<-dat$fac1L60overall + dat$fac2L60overall
dat$facboth60overall<-dat$fac1both60overall + dat$fac2both60overall
dat$minscale<-as.numeric('0')
dat$maxscale<-as.numeric('2')

#use custom function to scale data columns to between 0-1
minmax <- function(x) {
    return((x- min(x)) /(max(x)-min(x)))}
dat[,124:132] <- minmax(dat[,124:132])

#write final data to file
#saveRDS(dat, file="dat_discrete20.rds")
saveRDS(dat, file="dat_discrete40.rds")
```


## **Figure Generation Code Chunks:**
#### Basic simulation for 10 taxa/10 trait example used in the general methods section of the paper:
```{r}
### Set up parameters for simulating the loadings matrix
ntax = 10 # the number of tips in the tree (btw 1-50)
ntraits = 10 # the number of traits (btw 1-50)
nfactor = 2 # the number of factors (stick with 2 for now)
### Parameters for simulating the loadings matrix
pfac1 = 0.50 # percent of traits factor 1 influences
nfac1 = floor(pfac1*ntraits)
pfac2 = 0.50 # percent of trait factor 2 influences
nfac2 = floor(pfac2*ntraits)
### Either load in phylogeny or simulate a new one
# Load in phylogeny
phyltree<-read.tree("path to directory/simtree.txt")
# Simulate new phylogeny 
#phyltree<-ape::rtree(ntax)
#phyltree$tip.label<-gsub("t","",as.character(phyltree$tip.label)) # drops the "t" from the tip names to allow for trimming
#phyltree$tip.label<-paste0("taxon",as.character(phyltree$tip.label)) # add the taxon labels back in for easier ID
#phyltree<-phyltree_paths(phyltree) # advisable for speed
#write.tree(phyltree, file="simtree.txt")
### Simulate root node
kappa0 = 1.0
factor.varcov = diag(nfactor)  # this creates the identity matrix in R
factor.root = mvrnorm(mu=rep(0,nfactor),Sigma=(1/kappa0)*factor.varcov)
### Simulate factor values at tips
myfactors = simulBMProcPhylTree(phyltree,X0 = factor.root,Sigma=factor.varcov)
colnames(myfactors)<-c("factor_1", "factor_2")
order<-c("taxon1","taxon2","taxon3","taxon4","taxon5","taxon6","taxon7","taxon8","taxon9","taxon10")
myfactors<-myfactors[(match(order, rownames(myfactors))), ]
### Specify loadings for n factors
myloadings = matrix(c(rep(1.0,nfac1),rep(0,ntraits-nfac1),rep(0,ntraits-nfac2),rep(1.0,nfac2)),nrow=nfactor,byrow=T) #change to 
colnames(myloadings)<-c("trait_1","trait_2","trait_3","trait_4","trait_5","trait_6","trait_7","trait_8","trait_9","trait_10")
row.names(myloadings)<-c("factor_1","factor_2")
### Simulate error
mean.mat = matrix(0,ncol=ntraits,nrow=ntax) # mean for error matrix
trait.errors = rep(1.0,ntraits) # variance of each trait
### Draw error from matrix multivariate normal distributions
myU = diag(ntax)
myV = solve(diag(trait.errors))
myerror = rmatnorm(s=1,M=mean.mat,U=myU,V=myV)
colnames(myerror)<-c("trait_1","trait_2","trait_3","trait_4","trait_5","trait_6","trait_7","trait_8","trait_9","trait_10")
row.names(myerror)<-c("taxon_1","taxon_2","taxon_3","taxon_4","taxon_5","taxon_6","taxon_7","taxon_8","taxon_9","taxon_10")
### Write down final data
simdat0 = as.data.table(myfactors %*% myloadings + myerror, keep.rownames="taxon") #gives the taxon column the proper name so Julia can recognize it
colnames(simdat0)<-c("taxon","trait_1","trait_2","trait_3","trait_4","trait_5","trait_6","trait_7","trait_8","trait_9","trait_10")
simdat0<-simdat0 %>% arrange(match(taxon, c("taxon1","taxon2","taxon3","taxon4","taxon5","taxon6","taxon7","taxon8","taxon9","taxon10")))
### Standardize final data - subtract mean & divide by SD on trait-by-trait basis using standardize function from effectsize package
simdat = standardize(simdat0,robust=F) 
### Write the final continuous trait dataset to csv
#write.csv(simdat,"simdat_newL7.csv", row.names=F) #this will delete the first column of numbers so we don't have to edit the csv before input into julia
### Also try out binning to make binary categorical variables
traits<-simdat
traits$trait_2<-factor(ifelse(traits$trait_2>0,"1","0"))
traits$trait_3<-factor(ifelse(traits$trait_3>0,"1","0"))
traits$trait_8<-factor(ifelse(traits$trait_8>0,"1","0"))
traits$trait_9<-factor(ifelse(traits$trait_9>0,"1","0"))
### Write the final trait dataset including binary variables to csv
#write.csv(simdat,"simdat_newL7.csv", row.names=F)
```

#### Plotting all the success measures side by side
```{r}
library(ggplot2)
library(ggpubr)
library(NatParksPalettes)

setwd("C:/Users/shelb/OneDrive/Desktop/Chapter 2 Stuff/Simulation/script files/results_filesNov25")
par(mfrow=c(2,2))
#par(mfrow = c(1, 1)) #reset to normal layout

#factor 1 overall
plot1<-ggplot(dat1, aes(ntraits, ntax, fill= L40overall)) + 
  geom_tile() +
  scale_fill_natparks_c("Denali", direction=1) +
  guides(fill=guide_colourbar(barwidth=0.5, barheight=8, title="L=0.4")) +
  coord_fixed() +
  ggtitle("newsimtree1")
plot1

plot2<-ggplot(dat2, aes(ntraits, ntax, fill= L40overall)) + 
  geom_tile() +
  scale_fill_natparks_c("Denali", direction=1) +
  guides(fill=guide_colourbar(barwidth=0.5, barheight=8, title="L=0.4")) +
  coord_fixed() +
  ggtitle("newsimtree2")
plot2

plot3<-ggplot(dat3, aes(ntraits, ntax, fill= L40overall)) + 
  geom_tile() +
  scale_fill_natparks_c("Denali", direction=1) +
  guides(fill=guide_colourbar(barwidth=0.5, barheight=8, title="L=0.4")) +
  coord_fixed() +
  ggtitle("simtree501")
plot3

plot4<-ggplot(dat4, aes(ntraits, ntax, fill= L40overall)) + 
  geom_tile() +
  scale_fill_natparks_c("Denali", direction=1) +
  guides(fill=guide_colourbar(barwidth=0.5, barheight=8, title="L=0.4")) +
  coord_fixed() +
  ggtitle("simtree502")
plot4

figure<-ggarrange(plot1, plot2, plot3, plot4,
                  ncol=2, nrow=2)
figure

```

#### Success heatmaps - quad heatmap plot with single legend
```{r}
library(ggplot2)
library(ggpubr)
library(NatParksPalettes)
library(gridExtra)
library(dplyr)

setwd("C:/Users/shelb/OneDrive/Desktop/Chapter 2 Stuff/Simulation/script files/results_filesNov25")

#load in result datasets
dat1<-readRDS("dat1_shrink.rds")
dat2<-readRDS("dat2_shrink.rds")
dat3<-readRDS("dat3_shrink.rds")
dat4<-readRDS("dat4_shrink.rds")

dat1iid<-readRDS("dat1_iid.rds")
dat2iid<-readRDS("dat2_iid.rds")
dat3iid<-readRDS("dat3_iid.rds")
dat4iid<-readRDS("dat4_iid.rds")

#set up new small datasets - shrinkage prior
dat12<-left_join(dat1, dat2, by=c("iteration","ntraits","ntax"))
dat34<-left_join(dat3, dat4, by=c("iteration","ntraits","ntax"))

#grab averages from the two datasets - shrinkage prior
dat12$avg30<-(dat12$fac30overall.x+dat12$fac30overall.y)/2
dat12$max30<-pmax(dat12$fac30overall.x,dat12$fac30overall.y)
dat12$avg40<-(dat12$fac40overall.x+dat12$fac40overall.y)/2
dat12$max40<-pmax(dat12$fac40overall.x,dat12$fac40overall.y)
dat12$avg60<-(dat12$fac60overall.x+dat12$fac60overall.y)/2
dat12$max60<-pmax(dat12$fac60overall.x,dat12$fac60overall.y)

dat12$avg30both<-(dat12$facboth30overall.x+dat12$facboth30overall.y)/2
dat12$max30both<-pmax(dat12$facboth30overall.x,dat12$facboth30overall.y)
dat12$avg40both<-(dat12$facboth40overall.x+dat12$facboth40overall.y)/2
dat12$max40both<-pmax(dat12$facboth40overall.x,dat12$facboth40overall.y)
dat12$avg60both<-(dat12$facboth60overall.x+dat12$facboth60overall.y)/2
dat12$max60both<-pmax(dat12$facboth60overall.x,dat12$facboth60overall.y)

dat34$avg30<-(dat34$fac30overall.x+dat34$fac30overall.y)/2
dat34$max30<-pmax(dat34$fac30overall.x,dat34$fac30overall.y)
dat34$avg40<-(dat34$fac40overall.x+dat34$fac40overall.y)/2
dat34$max40<-pmax(dat34$fac40overall.x,dat34$fac40overall.y)
dat34$avg60<-(dat34$fac60overall.x+dat34$fac60overall.y)/2
dat34$max60<-pmax(dat34$fac60overall.x,dat34$fac60overall.y)

dat34$avg30both<-(dat34$facboth30overall.x+dat34$facboth30overall.y)/2
dat34$max30both<-pmax(dat34$facboth30overall.x,dat34$facboth30overall.y)
dat34$avg40both<-(dat34$facboth40overall.x+dat34$facboth40overall.y)/2
dat34$max40both<-pmax(dat34$facboth40overall.x,dat34$facboth40overall.y)
dat34$avg60both<-(dat34$facboth60overall.x+dat34$facboth60overall.y)/2
dat34$max60both<-pmax(dat34$facboth60overall.x,dat34$facboth60overall.y)

#set up new datasets - iid prior
dat12iid<-left_join(dat1iid, dat2iid, by=c("ntraits","ntax"))
dat34iid<-left_join(dat3iid, dat4iid, by=c("ntraits","ntax"))

#grab averages from the two datasets - iid prior
dat12iid$avg30<-(dat12iid$fac30overall.x+dat12iid$fac30overall.y)/2
dat12iid$max30<-pmax(dat12iid$fac30overall.x,dat12iid$fac30overall.y)
dat12iid$avg40<-(dat12iid$fac40overall.x+dat12iid$fac40overall.y)/2
dat12iid$max40<-pmax(dat12iid$fac40overall.x,dat12iid$fac40overall.y)
dat12iid$avg60<-(dat12iid$fac60overall.x+dat12iid$fac60overall.y)/2
dat12iid$max60<-pmax(dat12iid$fac60overall.x,dat12iid$fac60overall.y)

dat12iid$avg30both<-(dat12iid$facboth30overall.x+dat12iid$facboth30overall.y)/2
dat12iid$max30both<-pmax(dat12iid$facboth30overall.x,dat12iid$facboth30overall.y)
dat12iid$avg40both<-(dat12iid$facboth40overall.x+dat12iid$facboth40overall.y)/2
dat12iid$max40both<-pmax(dat12iid$facboth40overall.x,dat12iid$facboth40overall.y)
dat12iid$avg60both<-(dat12iid$facboth60overall.x+dat12iid$facboth60overall.y)/2
dat12iid$max60both<-pmax(dat12iid$facboth60overall.x,dat12iid$facboth60overall.y)

dat34iid$avg30<-(dat34iid$fac30overall.x+dat34iid$fac30overall.y)/2
dat34iid$max30<-pmax(dat34iid$fac30overall.x,dat34iid$fac30overall.y)
dat34iid$avg40<-(dat34iid$fac40overall.x+dat34iid$fac40overall.y)/2
dat34iid$max40<-pmax(dat34iid$fac40overall.x,dat34iid$fac40overall.y)
dat34iid$avg60<-(dat34iid$fac60overall.x+dat34iid$fac60overall.y)/2
dat34iid$max60<-pmax(dat34iid$fac60overall.x,dat34iid$fac60overall.y)

dat34iid$avg30both<-(dat34iid$facboth30overall.x+dat34iid$facboth30overall.y)/2
dat34iid$max30both<-pmax(dat34iid$facboth30overall.x,dat34iid$facboth30overall.y)
dat34iid$avg40both<-(dat34iid$facboth40overall.x+dat34iid$facboth40overall.y)/2
dat34iid$max40both<-pmax(dat34iid$facboth40overall.x,dat34iid$facboth40overall.y)
dat34iid$avg60both<-(dat34iid$facboth60overall.x+dat34iid$facboth60overall.y)/2
dat34iid$max60both<-pmax(dat34iid$facboth60overall.x,dat34iid$facboth60overall.y)

#plot the 4 heatmaps side by side - cutoff is the L=0.4
par(mfrow=c(2,2))

plot1<-ggplot(dat12, aes(ntraits, ntax, fill= max30both)) + 
  geom_tile() +
  scale_fill_natparks_c("Denali", direction=1) +
  guides(fill=guide_colourbar(barwidth=0.5, barheight=8, title=expression(paste("Success (", italic(" S"), ")")))) +
  coord_fixed() +
  ggtitle("a) Shrinkage Prior: \nNew N-tip Phylogeny Simulated Each Time") + xlab("Number of Traits (P)") + ylab("Number of Taxa (N)") +
  theme(plot.title = element_text(size=13), legend.position = "none")
#plot1

plot2<-ggplot(dat34, aes(ntraits, ntax, fill= max30both)) + 
  geom_tile() +
  scale_fill_natparks_c("Denali", direction=1) +
  guides(fill=guide_colourbar(barwidth=0.5, barheight=8, title=expression(paste("Success (", italic(" S"), ")")))) +
  coord_fixed() +
  ggtitle("b) Shrinkage Prior: \nBase 50-Tip Phylogeny Trimmed Each Time") + xlab("Number of Traits (P)") + ylab("Number of Taxa (N)") +
  theme(plot.title = element_text(size=13), legend.position = "none")
#plot2

plot3<-ggplot(dat12iid, aes(ntraits, ntax, fill= max30both)) + 
  geom_tile() +
  scale_fill_natparks_c("Denali", direction=1) +
  guides(fill=guide_colourbar(barwidth=0.5, barheight=8, title=expression(paste("Success (", italic(" S"), ")")))) +
  coord_fixed() +
  ggtitle("c) IID Prior: \nNew N-tip Phylogeny Simulated Each Time") + xlab("Number of Traits (P)") + ylab("Number of Taxa (N)") +
  theme(plot.title = element_text(size=13), legend.position = "none")
#plot3

plot4<-ggplot(dat34iid, aes(ntraits, ntax, fill= max30both)) + 
  geom_tile() +
  scale_fill_natparks_c("Denali", direction=1) +
  guides(fill=guide_colourbar(barwidth=8, barheight=0.5, label.position="bottom", title.position="top", title.hjust=0.5, title=expression(paste("Success (", italic(" S"), ")")))) +
  coord_fixed() +
  ggtitle("d) IID Prior: \nBase 50-Tip Phylogeny Trimmed Each Time") + xlab("Number of Traits (P)") + ylab("Number of Taxa (N)") +
  theme(plot.title = element_text(size=13), legend.position="none")
  #theme(plot.title = element_text(size=13), legend.position="bottom", legend.direction="horizontal")
#plot4

legend_plotv<-ggplot(dat34iid, aes(ntraits, ntax, fill= avg40)) +
  geom_tile() +
  scale_fill_natparks_c("Denali", direction=1) +
  guides(fill=guide_colourbar(barwidth=8, barheight=0.5, label.position="bottom", title.position="top", title.hjust=0.5, title=expression(paste("Success (", italic(" S"), ")")))) +
  theme(plot.title = element_text(size=13), legend.position="bottom", legend.direction="horizontal")
#legend_plotv

legend_ploth<-ggplot(dat34iid, aes(ntraits, ntax, fill= avg40)) +
  geom_tile() +
  scale_fill_natparks_c("Denali", direction=1) +
  guides(fill=guide_colourbar(barwidth=0.5, barheight=8, title.position="top", title.hjust=0.5, title=expression(paste("Success (", italic(" S"), ")")))) +
  theme(plot.title = element_text(size=13), legend.direction="vertical")
#legend_ploth

#extract the legend and add it back to all the plots
get_only_legend <- function(plot) {
# get tabular interpretation of plot
plot_table <- ggplot_gtable(ggplot_build(plot)) 
#  Mark only legend in plot
legend_plot <- which(sapply(plot_table$grobs, function(x) x$name) == "guide-box") 
# extract legend
legend <- plot_table$grobs[[legend_plot]]
# return legend
return(legend) 
}
legendv<-get_only_legend(legend_plotv)
legendh<-get_only_legend(legend_ploth)

combined_plot<-ggarrange(plot1, plot2, plot3, plot4, ncol=2, nrow=2)
figure2v<-grid.arrange(combined_plot, legendv, nrow=2, heights=c(10,1))
figure2h<-grid.arrange(combined_plot, legendh, ncol=2, widths=c(10,1.5))
figure2v
figure2h
```

#### Success measure for error simulation
```{r}
library(ggplot2)
library(ggpubr)
library(NatParksPalettes)
library(extrafont)
library(showtext)

#set up dataset for plotting
setwd("C:/Users/shelb/OneDrive/Desktop/Chapter 2 Stuff/Simulation/script files/results_filesNov25")
dat<-readRDS("dat_error.rds")
dat$Scenario<-c("A","B","C","D","E")

#plot barchart for different sub-scenarios
plot<-ggplot(dat, aes(x=Scenario, y=fac40overall)) + 
  geom_bar(stat = "identity", fill="#A4548B", width=0.55, alpha=0.75) + #"#977294" #F58551 #7F3C6A
  ylim(0, 1) +
  ggtitle("") +
  ylab("Success") +
  theme(text = element_text(family="Arial", face="bold", size=15, colour="black")) +
  theme(axis.title.x = element_text(vjust=-1)) 
  
plot #550x300
```

#### Success measure plot for discrete simulation
```{r}
library(ggplot2)
library(ggpubr)
library(NatParksPalettes)
library(extrafont)
library(showtext)

#set up datasets for plotting
setwd("C:/Users/shelb/OneDrive/Desktop/Chapter 2 Stuff/Simulation/script files/results_filesNov25")
dat20<-readRDS("dat_discrete20.rds")
dat20$Scenario<-c("A","B","C","D","E","F","G","H","I")
dat20$Dataset<-"Both N&P=20"
dat40<-readRDS("dat_discrete40.rds")
dat40$Scenario<-c("A","B","C","D","E","F","G","H","I")
dat40$Dataset<-"Both N&P=40"
datcombo<-rbind(dat20, dat40)

#plot side-by-side barcharts for the two sub-scenarios
plotcombo<-ggplot(datcombo, aes(x=Scenario, y=fac40overall, fill=Dataset)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.75, alpha=0.75) +
  scale_fill_manual(values=c("#34909A","#A4548B")) + #"#977294","#F58551" #"#2D6C73","#7F3C6A"
  ylim(0,1) +
  ggtitle("") +
  labs(fill="") +
  ylab("Success") +
  theme(legend.position="top", legend.key.spacing.x = unit(1, "cm")) +
  theme(text = element_text(family="Arial", face="bold", size=15, colour="black"))
 
plotcombo #550x350
```

#### Loadings plots (code adapted from Hassler et al 2022) - example loadings plots
```{r}
#### Loadings plots (code adapted from Hassler et al pfa paper)
library(ggplot2)
library(wesanderson)
library(colorspace)
library(tidyr)
library(extrafont)
library(showtext)

#load in fonts for use
#font_import()
#loadfonts(device="win")
#fonts()

font_add(family = "Bahnschrift", regular = "Bahnschrift.ttf")
showtext_auto()

#set up functions for plotting
f1 <- function(x) {
  return(x^2.5)
}

f2 <- function(x) {
  return(x^0.25)
}

custom_color_scale <- function(low, mid, high, f, n = 10) {
  lows = rev(half_scale(mid, low, f, n))
  highs = half_scale(mid, high, f, n)
  return(c(lows, hex(mid), highs))
}

half_scale <- function(c1, c2, f, n) {
  colors <- character(n)
  for (i in 1:n) {
    colors[i] <- hex(mixcolor(f(i / n), c1, c2))
  }
  return(colors)
}

##for blue & orange
#x <- 255
#c11 <- 100 #133
#c12 <- 190 #212
#c13 <- 200 #227
#c21 <- 250 #244
#c22 <- 181
#c23 <- 100 #189

##for green & purple
#x <- 255
#c11 <- 41 #133
#c12 <- 207 #212
#c13 <- 94 #227
#c21 <- 154 #244
#c22 <- 55
#c23 <- 184 #189

#for dark blue and light teal/green
x <- 255
c11 <- 0 #12 rgb(0% 21.9% 22.3%) new dark blue: 0, 52, 82
c12 <- 52 #23
c13 <- 82 #94
c21 <- 124 #28 rgb(0% 86.1% 71.1%) new light green: 131, 234, 167
c22 <- 226 #153
c23 <- 161 #101

load_colors <- custom_color_scale(sRGB(c11/x, c12/x, c13/x), sRGB(0.75, 0.75, 0.75), sRGB(c21/x, c22/x, c23/x), f1)
fac_colors <- custom_color_scale(sRGB(1, 0, 0), sRGB(1, 1, 1), sRGB(0, 0, 1), f2)


# colors <- c("#0000FF", "#1515EA", "#2B2BD5", "#4040C0", "#5555AA", "#6B6B95", "#808080", "#956B6B", )

## Loadings plot

plot_loadings <- function(csv_path, plot_name, labels_path = NA, factors = NA,
                          height_scale=1.0, width_scale=1.0,
                          verbose = FALSE,
                          lims = NA){
  
  if (verbose) {print("Starting loadings plot")}
  if (verbose) {print("Reading data")}
  
  df  <- read.csv(csv_path, header=TRUE, encoding="UTF-8")
  if (!all(is.na(factors))) {
    if (verbose) {print("Subsetting factors")}
    df <- df[df$factor %in% factors,]
  }
  
  trait_levs <- unique(df$trait)
  trait_labels <- c("Trait1","Trait2","Trait3","Trait4","Trait5","Trait6","Trait7","Trait8","Trait9","Trait10")
  
  if (verbose) {print("X")}
  
  if (!is.na(labels_path)) {
    labels_df <- read.csv(labels_path, header=TRUE, fileEncoding="UTF-8-BOM")
    for (i in 1:length(trait_levs)[[1]]) {
      trait_labels[labels_df$trait[i]] = labels_df$pretty[i]
    }
    
    #trait_levs <- labels_df$pretty
    cat_levs <- unique(labels_df$cat)
    
    cat_dict <- labels_df$cat
    names(cat_dict) <- labels_df$trait
    df$cat <- cat_dict[df$trait]
    
    # names(trait_levs) <- labels_df$trait
    
    # df$trait <- trait_levs[df$trait]
  } else {
    for (i in 1:length(trait_levs)[[1]]) {
      trait_labels[trait_levs[i]] = trait_levs[i]
    }
    df$cat <- rep("NA", nrow(df))
    cat_levs <- c("NA")
  }
  
  
  wes_pal <- wes_palette("GrandBudapest2")
  pal <- c(wes_pal[1], "grey", wes_pal[4])
  
  df$trait <- factor(df$trait, levels=trait_levs)
  df$cat <- factor(df$cat, levels=cat_levs)
  df$L <- sapply(df$L, as.numeric)
  df$sign_perc <- df$perc
  for (i in 1:length(df$perc)) {
    if (df$L[i] < 0.0) {
      df$sign_perc[i] <- 1.0 - df$sign_perc[i]
    }
    if (df$hpdu[i] == 0 && df$hpdl[i] == 0) {
      df$sign_perc[i] = NA
    }
  }
  
  df$sign <- sapply(df$L, sign)
  n <- dim(df)[[1]]
  for (i in 1:n) {
    if (df$hpdl[i] <= 0 && df$hpdu[i] >= 0) {
      df$sign[i] = 0
    }
  }
  df$sign <- factor(df$sign, levels=c(1,0, -1))
  ymin = min(df$hpdl)
  ymax = max(df$hpdu)
  if (is.na(lims)) {
    absmax = max(abs(ymin), abs(ymax))
    lims <- c(-absmax, absmax)
  }
  
  k <- max(df$factor)
  facet_labels <- character(k)
  facet_names <- integer(k)
  for (i in 1:k) {
    facet_labels[i] = paste("Factor", i)
    facet_names[i] = i
  }
  names(facet_labels) <- facet_names
  
  n_groups <- length(cat_levs)
  group_labels <- character(n_groups)
  if (length(n_groups) == 1) {
    names(group_labels) <- cat_levs
  }
  
  # ps <- c()
  # for (i in 1:max(df$row)) {
  #   ps <- c(ps, plot_single_row(df, 1, ymin, ymax))
  # }
  
  p <- ggplot(df) +
    geom_vline(xintercept=0, linetype="dashed") +
    geom_vline(xintercept=0.4, linetype="dashed", color='maroon') + #this is the cutoff for loadings that we want to show
    #geom_vline(xintercept=-0.4, linetype="dashed", color='maroon') +
    geom_point(aes(y=trait, x=L, color=sign_perc), size=2.5) +
    geom_errorbarh(aes(y=trait, xmin=hpdl, xmax=hpdu, color=sign_perc), height=0.0, linewidth=2) +
    scale_color_gradientn(colors = load_colors, limits=c(0, 1), name="p>0", na.value="grey75") + #this was where probability > 0
    scale_y_discrete(limits=rev, labels=trait_labels) +
    #scale_color_gradient2(low="orange", mid="white", high="purple", limits=c(-1, 1), name="L") +
    facet_grid(~ cat, scales="free_x", space="free_x") +
    #geom_tile() +
    #scale_fill_gradient2(low="orange", mid="white", high="purple", midpoint=0) +
    #scale_x_discrete(position = "top") +
    #scale_y_discrete() +
    labs(y="", x="Loadings Value") +
    theme_minimal() +
    theme(#axis.text.x = element_text(angle=0, hjust=1),
      text = element_text(family="Bahnschrift", face="bold", size=25, colour="black"),
      panel.border = element_rect(colour = "black", fill=NA),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.position = "right",
    ) +
    xlim(lims) +
    facet_grid(rows=vars(cat), #NULL if I want no labels to appear
               cols=vars(factor),
               labeller = labeller(factor=facet_labels),
               scales="free_y",
               space="free_y",
               switch = "y")
  if (length(cat_levs) == 1) {
    p <- p + theme(strip.text.y = element_blank())
  }
  
  # axis.title.y = element_text())
  n_traits <- length(trait_levs)
  ggsave(plot_name, width=width_scale * (k * 2 + 2), height= height_scale * (n_traits * 0.15 + 1), units="in", limitsize=FALSE)
  gc()
  return(p)
}

#example figures for paper
#set path to directory containing example files
my_dir = "C:/Users/shelb/OneDrive/Desktop/Chapter 2 Stuff/Simulation/script files/figures/example loadings matrices"

#good example
labels_path = file.path(my_dir, "loadingsStats_ex.csv")
stats_path = file.path(my_dir, "loadingsStats_ex.csv")
plot_loadings(stats_path, "test_ex.pdf", labels_path = labels_path)

#label switching example
labels_path = file.path(my_dir, "loadingsStats_ls.csv")
stats_path = file.path(my_dir, "loadingsStats_ls.csv")
plot_loadings(stats_path, "test_ls.pdf", labels_path = labels_path)

#bad example
labels_path = file.path(my_dir, "loadingsStats_bad.csv")
stats_path = file.path(my_dir, "loadingsStats_bad.csv")
plot_loadings(stats_path, "test_bad.pdf", labels_path = labels_path)
```

#### Loadings plots (code adapted from Hassler et al 2022) - simulation and emprirical loadings plots
```{r}
library(ggplot2)
library(wesanderson)
library(colorspace)
library(tidyr)
library(extrafont)
library(showtext)

#set up functions for plotting
f1 <- function(x) {
  return(x^2.5)
}

f2 <- function(x) {
  return(x^0.25)
}

custom_color_scale <- function(low, mid, high, f, n = 10) {
  lows = rev(half_scale(mid, low, f, n))
  highs = half_scale(mid, high, f, n)
  return(c(lows, hex(mid), highs))
}

half_scale <- function(c1, c2, f, n) {
  colors <- character(n)
  for (i in 1:n) {
    colors[i] <- hex(mixcolor(f(i / n), c1, c2))
  }
  return(colors)
}

##for blue & orange
#x <- 255
#c11 <- 100 #133
#c12 <- 190 #212
#c13 <- 200 #227
#c21 <- 250 #244
#c22 <- 181
#c23 <- 100 #189

##for green & purple
#x <- 255
#c11 <- 41 #133
#c12 <- 207 #212
#c13 <- 94 #227
#c21 <- 154 #244
#c22 <- 55
#c23 <- 184 #189

#for dark blue and light teal/green
x <- 255
c11 <- 21 #12 rgb(0% 21.9% 22.3%) new dark blue: 0, 52, 82; new blue: 18, 69, 89; new new blue: 14,52,68
c12 <- 71 #23
c13 <- 81 #94
c21 <- 65 #28 rgb(0% 86.1% 71.1%) new light green: 131, 234, 167; new light color: 158, 188, 159; new new light: 147,180,148
c22 <- 220 #153
c23 <- 145 #101

load_colors <- custom_color_scale(sRGB(c11/x, c12/x, c13/x), sRGB(0.75, 0.75, 0.75), sRGB(c21/x, c22/x, c23/x), f1)
fac_colors <- custom_color_scale(sRGB(1, 0, 0), sRGB(1, 1, 1), sRGB(0, 0, 1), f2)


# colors <- c("#0000FF", "#1515EA", "#2B2BD5", "#4040C0", "#5555AA", "#6B6B95", "#808080", "#956B6B", )

## Loadings plot

plot_loadings <- function(csv_path, plot_name, labels_path = NA, factors = NA,
                          height_scale=1.0, width_scale=1.0,
                          verbose = FALSE,
                          lims = NA){
  
  if (verbose) {print("Starting loadings plot")}
  if (verbose) {print("Reading data")}
  
  df  <- read.csv(csv_path, header=TRUE, encoding="UTF-8")
  if (!all(is.na(factors))) {
    if (verbose) {print("Subsetting factors")}
    df <- df[df$factor %in% factors,]
  }
  
  trait_levs <- unique(df$trait)
  trait_labels <- c()
  
  if (verbose) {print("X")}
  
  if (!is.na(labels_path)) {
    labels_df <- read.csv(labels_path, header=TRUE, fileEncoding="UTF-8-BOM")
    for (i in 1:length(trait_levs)[[1]]) {
      trait_labels[labels_df$trait[i]] = labels_df$pretty[i]
    }
    
    # trait_levs <- labels_df$pretty
    cat_levs <- unique(labels_df$cat)
    
    cat_dict <- labels_df$cat
    names(cat_dict) <- labels_df$trait
    df$cat <- cat_dict[df$trait]
    
    # names(trait_levs) <- labels_df$trait
    
    # df$trait <- trait_levs[df$trait]
  } else {
    for (i in 1:length(trait_levs)[[1]]) {
      trait_labels[trait_levs[i]] = trait_levs[i]
    }
    df$cat <- rep("NA", nrow(df))
    cat_levs <- c("NA")
  }
  
  
  wes_pal <- wes_palette("GrandBudapest2")
  pal <- c(wes_pal[1], "grey", wes_pal[4])
  
  df$trait <- factor(df$trait, levels=trait_levs)
  df$cat <- factor(df$cat, levels=cat_levs)
  df$L <- sapply(df$L, as.numeric)
  df$sign_perc <- df$perc
  for (i in 1:length(df$perc)) {
    if (df$L[i] < 0.0) {
      df$sign_perc[i] <- 1.0 - df$sign_perc[i]
    }
    if (df$hpdu[i] == 0 && df$hpdl[i] == 0) {
      df$sign_perc[i] = NA
    }
  }
  
  df$sign <- sapply(df$L, sign)
  n <- dim(df)[[1]]
  for (i in 1:n) {
    if (df$hpdl[i] <= 0 && df$hpdu[i] >= 0) {
      df$sign[i] = 0
    }
  }
  df$sign <- factor(df$sign, levels=c(1,0, -1))
  ymin = min(df$hpdl)
  ymax = max(df$hpdu)
  if (is.na(lims)) {
    absmax = max(abs(ymin), abs(ymax))
    #lims <- c(-absmax, absmax)
    lims <- c(-2, 2)
  }
  
  k <- max(df$factor)
  facet_labels <- character(k)
  facet_names <- integer(k)
  for (i in 1:k) {
    facet_labels[i] = paste("Factor", i)
    facet_names[i] = i
  }
  names(facet_labels) <- facet_names
  
  n_groups <- length(cat_levs)
  group_labels <- character(n_groups)
  if (length(n_groups) == 1) {
    names(group_labels) <- cat_levs
  }
  
  # ps <- c()
  # for (i in 1:max(df$row)) {
  #   ps <- c(ps, plot_single_row(df, 1, ymin, ymax))
  # }
  
  p <- ggplot(df) +
    geom_vline(xintercept=0, linetype="dashed") +
    geom_vline(xintercept=0.4, linetype="dashed", color='maroon') + #this is the cutoff for loadings that we want to show
    #geom_vline(xintercept=-0.4, linetype="dashed", color='maroon') +    
    geom_point(aes(y=trait, x=L, color=sign_perc), size=3) +
    geom_errorbarh(aes(y=trait, xmin=hpdl, xmax=hpdu, color=sign_perc), height=0.0, linewidth=2) +
    scale_color_gradientn(colors = load_colors, limits=c(0, 1), name="p > 0", na.value="grey75") +
    scale_y_discrete(limits=rev, labels=trait_labels) +
    # scale_color_gradient2(low="orange", mid="white", high="purple", limits=c(-1, 1), name="L") +
    #facet_grid(~ cat, scales="free_x", space="free_x") +
    #geom_tile() +
    #scale_fill_gradient2(low="orange", mid="white", high="purple", midpoint=0) +
    #scale_x_discrete(position = "top") +
    # scale_y_discrete() +
    labs(y="", x="loadings value") +
    theme_minimal() +
    theme(#axis.text.x = element_text(angle=0, hjust=1),
      text = element_text(family="Arial", face="bold", size=18, colour="black"),
      panel.border = element_rect(colour = "black", fill=NA),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.position = "right"
    ) +
    xlim(lims) +
    facet_grid(rows=vars(cat),
               cols=vars(factor),
               labeller = labeller(factor=facet_labels),
               scales="free_y",
               space="free_y",
               switch = "y")
  if (length(cat_levs) == 1) {
    p <- p + theme(strip.text.y = element_blank())
  }
  
  # axis.title.y = element_text())
  n_traits <- length(trait_levs)
  ggsave(plot_name, width=width_scale * (k * 2 + 2), height= height_scale * (n_traits * 0.15 + 1), units="in", limitsize=FALSE)
  gc()
  return(p)
}

###buildabat loadings summaries
#set path to directory containing all results files
my_dir = "C:/Users/shelb/OneDrive/Desktop/Chapter 2 Stuff/Simulation/script files/results_filesNov25"

#all continuous data
labels_path = file.path(my_dir, "buildabat_loadingsstats1.csv")
stats_path = file.path(my_dir, "buildabat_loadingsstats1.csv")
plot_loadings(stats_path, "buildabat_cont.pdf", labels_path = labels_path)

#some traits converted to binary
labels_path = file.path(my_dir, "buildabat_loadingsstats2.csv")
stats_path = file.path(my_dir, "buildabat_loadingsstats2.csv")
plot_loadings(stats_path, "buildabat_disc.pdf", labels_path = labels_path)

#all traits combined
labels_path = file.path(my_dir, "buildabat_loadingstatssum.csv")
stats_path = file.path(my_dir, "buildabat_loadingstatssum.csv")
plot_loadings(stats_path, "buildabat_all.pdf", labels_path = labels_path)

###error simulation loadings summaries
#scenario 1 - traits 11-20 with larger error
labels_path = file.path(my_dir, "error_loadingsStats1.csv")
stats_path = file.path(my_dir, "error_loadingsStats1.csv")
plot_loadings(stats_path, "error_scenario1.pdf", labels_path = labels_path) #550x600

#scenario 3 - traits 9/10 and 19/20 with even larger error
labels_path = file.path(my_dir, "error_loadingsStats3.csv")
stats_path = file.path(my_dir, "error_loadingsStats3.csv")
plot_loadings(stats_path, "error_scenario3.pdf", labels_path = labels_path) #550x600

#scenario 5 - traits 8/13/18 and 28/33/38 with random values
labels_path = file.path(my_dir, "error_loadingsStats5.csv")
stats_path = file.path(my_dir, "error_loadingsStats5.csv")
plot_loadings(stats_path, "error_scenario5.pdf", labels_path = labels_path) #550x900

###discrete simulation loadings summaries
#scenario 20C - 5 binary traits per factor
labels_path = file.path(my_dir, "discrete20_loadingsStats2.csv")
stats_path = file.path(my_dir, "discrete20_loadingsStats2.csv")
plot_loadings(stats_path, "discret20_scenarioC.pdf", labels_path = labels_path) #550x450

#scenario 40C - 5 binary traits per factor
labels_path = file.path(my_dir, "discrete40_loadingsStats2.csv")
stats_path = file.path(my_dir, "discrete40_loadingsStats2.csv")
plot_loadings(stats_path, "discrete40_scenarioC.pdf", labels_path = labels_path) #550x900

###empirical loadings summaries
#myotis (gunnell tree)
labels_path = file.path(my_dir, "myotisgunnell_loadingsstats.csv")
stats_path = file.path(my_dir, "myotisgunnell_loadingsstats.csv")
plot_loadings(stats_path, "myotis_gunnell.pdf", labels_path = labels_path)

#phyllostomidae (upham tree)
labels_path = file.path(my_dir, "phylloupham_loadingsstats.csv")
stats_path = file.path(my_dir, "phylloupham_loadingsstats.csv")
plot_loadings(stats_path, "phyllo_upham.pdf", labels_path = labels_path)
```

#### Loadings plots (code adapted from Hassler et al 2022) - buildabat combination plot
```{r}
library(ggplot2)
library(wesanderson)
library(colorspace)
library(tidyr)
library(extrafont)
library(showtext)

#set up functions for plotting
f1 <- function(x) {
  return(x^2.5)
}

f2 <- function(x) {
  return(x^0.25)
}

custom_color_scale <- function(low, mid, high, f, n = 10) {
  lows = rev(half_scale(mid, low, f, n))
  highs = half_scale(mid, high, f, n)
  return(c(lows, hex(mid), highs))
}

half_scale <- function(c1, c2, f, n) {
  colors <- character(n)
  for (i in 1:n) {
    colors[i] <- hex(mixcolor(f(i / n), c1, c2))
  }
  return(colors)
}

#for dark blue and light teal/green
x <- 255
c11 <- 21 #12 rgb(0% 21.9% 22.3%) new dark blue: 0, 52, 82; new blue: 18, 69, 89; new new blue: 14,52,68
c12 <- 71 #23
c13 <- 81 #94
c21 <- 65 #28 rgb(0% 86.1% 71.1%) new light green: 131, 234, 167; new light color: 158, 188, 159; new new light: 147,180,148
c22 <- 220 #153
c23 <- 145 #101

load_colors <- custom_color_scale(sRGB(c11/x, c12/x, c13/x), sRGB(0.75, 0.75, 0.75), sRGB(c21/x, c22/x, c23/x), f1)
fac_colors <- custom_color_scale(sRGB(1, 0, 0), sRGB(1, 1, 1), sRGB(0, 0, 1), f2)


# colors <- c("#0000FF", "#1515EA", "#2B2BD5", "#4040C0", "#5555AA", "#6B6B95", "#808080", "#956B6B", )

## Loadings plot

plot_loadings <- function(csv_path, plot_name, labels_path = NA, factors = NA,
                          height_scale=1.0, width_scale=1.0,
                          verbose = FALSE,
                          lims = NA){
  
  if (verbose) {print("Starting loadings plot")}
  if (verbose) {print("Reading data")}
  
  df  <- read.csv(csv_path, header=TRUE, encoding="UTF-8")
  if (!all(is.na(factors))) {
    if (verbose) {print("Subsetting factors")}
    df <- df[df$factor %in% factors,]
  }
  
  trait_levs <- unique(df$trait)
  trait_labels <- c()
  
  if (verbose) {print("X")}
  
  if (!is.na(labels_path)) {
    labels_df <- read.csv(labels_path, header=TRUE, fileEncoding="UTF-8-BOM")
    for (i in 1:length(trait_levs)[[1]]) {
      trait_labels[labels_df$trait[i]] = labels_df$pretty[i]
    }
    
    # trait_levs <- labels_df$pretty
    cat_levs <- unique(labels_df$cat)
    
    cat_dict <- labels_df$cat
    names(cat_dict) <- labels_df$trait
    df$cat <- cat_dict[df$trait]
    
    # names(trait_levs) <- labels_df$trait
    
    # df$trait <- trait_levs[df$trait]
  } else {
    for (i in 1:length(trait_levs)[[1]]) {
      trait_labels[trait_levs[i]] = trait_levs[i]
    }
    df$cat <- rep("NA", nrow(df))
    cat_levs <- c("NA")
  }
  
  
  wes_pal <- wes_palette("GrandBudapest2")
  pal <- c(wes_pal[1], "#f1f9f7", wes_pal[4])
  
  df$trait <- factor(df$trait, levels=trait_levs)
  df$cat <- factor(df$cat, levels=cat_levs)
  df$L <- sapply(df$L, as.numeric)
  df$sign_perc <- df$perc
  for (i in 1:length(df$perc)) {
    if (df$L[i] < 0.0) {
      df$sign_perc[i] <- 1.0 - df$sign_perc[i]
    }
    if (df$hpdu[i] == 0 && df$hpdl[i] == 0) {
      df$sign_perc[i] = NA
    }
  }
  
  df$sign <- sapply(df$L, sign)
  n <- dim(df)[[1]]
  for (i in 1:n) {
    if (df$hpdl[i] <= 0 && df$hpdu[i] >= 0) {
      df$sign[i] = 0
    }
  }
  df$sign <- factor(df$sign, levels=c(1,0, -1))
  ymin = min(df$hpdl)
  ymax = max(df$hpdu)
  if (is.na(lims)) {
    absmax = max(abs(ymin), abs(ymax))
    lims <- c(-absmax, absmax)
  }
  
  k <- max(df$factor)
  facet_labels <- character(k)
  facet_names <- integer(k)
  for (i in 1:k) {
    facet_labels[i] = paste("Factor", i)
    facet_names[i] = i
  }
  names(facet_labels) <- facet_names
  
  n_groups <- length(cat_levs)
  group_labels <- character(n_groups)
  if (length(n_groups) == 1) {
    names(group_labels) <- cat_levs
  }
  
  # ps <- c()
  # for (i in 1:max(df$row)) {
  #   ps <- c(ps, plot_single_row(df, 1, ymin, ymax))
  # }
  
  p <- ggplot(df) +
    geom_vline(xintercept=0, linetype="dashed", color="black") + #this is the 0 cutoff
    geom_vline(xintercept=0.4, linetype="dashed", color='#B62B5E') + #this is the cutoff for loadings that we want to show
    #geom_vline(xintercept=-0.4, linetype="dashed", color='maroon') +    
    geom_point(aes(y=trait, x=L, color=sign_perc), size=3) +
    geom_errorbarh(aes(y=trait, xmin=hpdl, xmax=hpdu, color=sign_perc), height=0.0, linewidth=2) +
    scale_color_gradientn(colors = load_colors, limits=c(0, 1), name="p > 0", na.value="#f1f9f7") +
    scale_y_discrete(limits=rev, labels=trait_labels) +
    # scale_color_gradient2(low="orange", mid="white", high="purple", limits=c(-1, 1), name="L") +
    #facet_grid(~ cat, scales="free_x", space="free_x") +
    #geom_tile() +
    #scale_fill_gradient2(low="orange", mid="white", high="purple", midpoint=0) +
    #scale_x_discrete(position = "top") +
    # scale_y_discrete() +
    labs(y="", x="loadings value") + #caption=c("All Continuous Traits","Including Discrete Traits")
    theme_minimal() +
    theme(#axis.text.x = element_text(angle=0, hjust=1),
      text = element_text(family="Arial", face="bold", size=18, colour="black"),
      panel.border = element_rect(colour = "black", fill=NA),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      #plot.caption=element_text(color=c("black","gray40"),hjust=0),
      legend.position = "right",
      axis.text.y=element_text(color=c("gray40","black"))
    ) +
    xlim(lims) +
    #annotate("text", x = 4, y=13000, label = "All Continuous Traits", color="black") + #add text to describe the two different colors
    #annotate("text", x = 4, y=13000, label = "Including Discrete Traits", color="gray40") +
    facet_grid(rows=vars(cat),
               cols=vars(factor),
               labeller = labeller(factor=facet_labels),
               scales="free_y",
               space="free_y",
               switch = "y")
  if (length(cat_levs) == 1) {
    p <- p + theme(strip.text.y = element_blank())
  }
  
  # axis.title.y = element_text())
  n_traits <- length(trait_levs)
  ggsave(plot_name, width=width_scale * (k * 2 + 2), height= height_scale * (n_traits * 0.15 + 1), units="in", limitsize=FALSE)
  gc()
  return(p)
}

#plot all traits together
#all traits combined
my_dir = "C:/Users/shelb/OneDrive/Desktop/Chapter 2 Stuff/Simulation/script files/compiled_results_May2025"
labels_path = file.path(my_dir, "buildabat_loadingstatssum.csv")
stats_path = file.path(my_dir, "buildabat_loadingstatssum.csv")
plot_loadings(stats_path, "buildabat_all.pdf", labels_path = labels_path)
#
```

#### Phylogeny/Factor Means plots (code adapted from Hassler et al 2022) - buildabat plot
```{r}
###Hassler factor-phylogeny plotting code
library(ggtree)
library(tidyr)
library(phytools)
library(tidytree)
library(ggplot2)
library(aplot)
library(RColorBrewer)
library(ggnewscale)
library(phyclust)
library(geiger)
require(treeio)
library(extrafont)
library(showtext)

#set up the functions for phylogeny plots

prep_trait <- function(tree, trait){
  fit <- phytools::fastAnc(tree, trait, vars=TRUE, CI=TRUE)
  td <- data.frame(node = nodeid(tree, names(trait)),
                   trait = trait)
  nd <- data.frame(node = names(fit$ace), trait = fit$ace)
  d <- rbind(td, nd)
  d$node <- as.numeric(d$node)
  tree <- full_join(tree, d, by = 'node')
  return(tree)
}

plot_tree <- function(tree, colors,
                      border=FALSE,
                      line_width=0.1,
                      tip_labels=FALSE,
                      layout="rectangular",
                      color_tree=TRUE,
                      fan.angle=15,
                      new_labels=NA,
                      limits=NA,
                      labels_offset=0) {
  p <- ggtree(tree, layout=layout, open.angle=fan.angle, size = 0)
  if (border){
    # p <- p + geom_tree(size=line_width)
  } else {
    #   p <- p + ggtree(tree)
  }
  if (color_tree) {
    p <- p + geom_tree(aes(color=trait),
                       continuous = TRUE, size=line_width)
    if (is.na(limits)) {
      #p <- p + scale_color_gradient2(midpoint = 0.0, low=colors[[1]], mid='grey', high=colors[[2]])
      p <- p + scale_color_gradient2(midpoint = 0.0, low=colors[[1]], mid='grey85', high=colors[[2]])
    } else {
      #p <- p + scale_color_gradient2(midpoint = 0.0, low=colors[[1]], mid='grey', high=colors[[2]], limits=limits)
      p <- p + scale_color_gradient2(midpoint = 0.0, low=colors[[1]], mid='grey85', high=colors[[2]], limits=limits)
      
    }
    
    p <- p + labs(color="Factor Value")
    
  } else {
    p <- p + geom_tree(size=line_width)
  }
  
  if (tip_labels) {
    if (!is.na(new_labels)) {
      p <- p %<+% new_labels + geom_tiplab(aes(label=classes, family="Bahnschrift", fontface="bold", size=18), offset=labels_offset)
    } else {
      p <- p + geom_tiplab(align=T, offset=labels_offset, geom="text", family="Bahnschrift", fontface="bold", size=4) #size=4, fontface=1
    }
  } else {
    # p <- p + geom_tiplab(aes(label=character(1)), align=TRUE)
  }
  
  # p <- ggtree(tree, size=2.8) +
  #   geom_tree(aes(color=trait),
  #             continuous = TRUE, size=2) +
  #   scale_color_gradient2(midpoint = 0, low='purple', mid='white', high='orange')
  return(p)
}

plot_factor_tree <- function(name, tree_path, factors_path, factors = NA,
                             class_path = NA,
                             height=NA,
                             width=NA,
                             border=FALSE,
                             line_width=1.0,
                             tip_labels=TRUE,
                             layout="rectangular",
                             class_palette=NA,
                             combined=TRUE,
                             fan.angle=30.0,
                             new_labels=NA,
                             common_scale=FALSE,
                             scale=TRUE,
                             heat_width=heat_width,
                             extra_offset=extra_offset,
                             legend_position=legend_position,
                             labels_offset=0.02,
                             fac_names=NA,
                             factor_fill = scale_fill_gradient2(midpoint = 0.0, low="#053439", mid='#F4F5F6', high="#440E33"),
                             factor_color = scale_color_gradient2(midpoint = 0.0, low="#053439", mid='#F4F5F6', high="#440E33"),
                             class_values=c("#729EA1","#B6D19E","#975E7C"),
                             legend_name="Example Class",
                             relabel = NA,
                             include_only = NA
) {
  
  x <- as.matrix(read.csv(factors_path, header=TRUE))
  k <- ncol(x) - 1
  if (all(is.na(factors))) {
    factors <- 1:k
  }
  
  base_tree <- read.tree(tree_path)
  taxa <- x[, 1]
  
  include_class = !is.na(class_path)
  if (include_class) {
    class_df = read.csv(class_path)
  }
  
  if (!all(is.na(include_only))) {
    drop_taxa = setdiff(taxa, include_only)
    keep_rows = match(include_only, taxa)
    x <- x[keep_rows, ]
    for (taxon in drop_taxa) {
      base_tree <- drop.tip(base_tree, taxon)
    }
    
    if (include_class) {
      keep_rows <- match(include_only, class_df$taxon)
      class_df <- class_df[keep_rows, ]
    }
  }
  
  taxa <- x[, 1]
  
  if (!all(is.na(relabel))) {
    matched_rows <- match(relabel$original, taxa)
    taxa[matched_rows] <- relabel$new
    
    matched_rows <- match(relabel$original, base_tree$tip.label)
    base_tree$tip.label[matched_rows] <- relabel$new
    
    if (include_class) {
      matched_rows <- match(relabel$original, class_df$taxon)
      class_df$taxon[matched_rows] <- relabel$new
    }
    
  }
  
  x <- x[, factors + 1]
  x <- apply(as.matrix(x), 2, as.numeric)
  rownames(x) <- taxa
  n <- length(taxa)
  if (scale) {
    x <- scale(x)
  }
  
  limits <- NA
  if (common_scale) {
    x_min = min(x)
    x_max = max(x)
    limits = c(x_min, x_max)
  }
  
  combined_layout = layout
  if (is.na(height)) {
    height <- 20
    if (layout == "rectangular") {
      height <- 0.2 * n
    }
  }
  
  if (is.na(width)) {
    width <- 20
  }
  if (layout=="circular") {
    combined_layout = "fan"
  }
  
    # include_class = !is.na(class_path)
  if (include_class) {
    
    # stop()
    class_df = read.csv(class_path)
    class_df[,2] <- as.factor(class_df[,2])
    classes <- data.frame(x = class_df[,2])
    colnames(classes) <- c(colnames(class_df)[2])
    row.names(classes) <- class_df$taxon
    
    class_fills = scale_fill_manual(values=class_values, name=legend_name)
    #class_fills = discrete_scale("fill", "manual", colorRampPalette(brewer.pal(min(nrow(unique(classes)), 12), "Set2")), name=colnames(classes)[1])
    
    class_colors = scale_color_manual(values=class_values, name=legend_name)
    #class_colors = discrete_scale("color", "manual", colorRampPalette(brewer.pal(min(nrow(unique(classes)), 12), "Set2")), name=colnames(classes)[1])
    
    if (!all(is.na(class_palette))) {
      #class_fills = scale_fill_manual(values=class_palette, name=colnames(classes)[1])
      #class_colors = scale_color_manual(values=class_palette, name=colnames(classes)[1])
    }
  }
  
  n_factors = length(factors)
  
  if (all(is.na(fac_names))) {
    fac_names <- character(n_factors)
    for (i in 1:n_factors) {
      fac_names[i] = paste("F", i)
    }
  }
  
  colnames(x) <- fac_names
  max_height <- max(node.depth.edgelength(base_tree))
  base_tree$edge.length <- base_tree$edge.length / max_height
  base_tree$edge.width <- base_tree$edge.width / (10*max_height)
  
  heat_width=heat_width

  
  if (combined) {
    legend_position=legend_position
    p1 <- plot_tree(base_tree, "", border=border, line_width=line_width, tip_labels = tip_labels, layout=combined_layout, color_tree=FALSE, new_labels=new_labels, labels_offset=labels_offset)
    pname <- paste(name, "_factors.svg", sep="")
    p2 <- p1 + theme(text = element_text(family="Arial", face="bold", size=10, colour="black"))
    if (include_class) {
      #p2 <- my_gheatmap(p1, classes, offset=extra_offset, width=heat_width, colnames_angle=90, colnames=TRUE, legend_title=colnames(classes)[1]) + class_fills #+ class_colors
      p2 <- my_gheatmap(p1, classes, offset=extra_offset, width=heat_width, colnames_angle=90, colnames=TRUE, legend_position=legend_position, legend_title="Ecotype") + 
        class_fills + class_colors +
        theme(text = element_text(family="Arial", face="bold", size=20, colour="black"))
    }
    
    p3 <- p2 + new_scale_fill() + new_scale_color()
    p4 <- my_gheatmap(p3, x, offset=extra_offset + heat_width * 2, width=heat_width * n_factors, colnames_angle=90, colnames_offset_y = 0, legend_position=legend_position, legend_title="Factor Value") +
      factor_fill + factor_color +
      labs(fill="Factor Value") +
      guides(color=FALSE) +
      theme(text = element_text(family="Arial", face="bold", size=12, colour="black"))
    
    svg(pname, height=height, width=width)
    print(p4)
    dev.off()
    gc()
    return(p4)
    
  }
  gc()
}


my_gheatmap <- function(p, data, legend_position=legend_position, offset=0, width=1, low="green", high="red", color="white",
                        colnames=T, colnames_position="bottom", colnames_angle=0, colnames_level=NULL, #colnames=TRUE for heatmap labeled
                        colnames_offset_x = 0, colnames_offset_y = -2, font.size=4, family="", hjust=0.5, legend_title = "value") {
  
  colnames_position %<>% match.arg(c("bottom", "top"))
  variable <- value <- lab <- y <- NULL
  
  ## convert width to width of each cell
  width <- width * (p$data$x %>% range(na.rm=TRUE) %>% diff) / ncol(data)
  
  isTip <- x <- y <- variable <- value <- from <- to <- NULL
  
  df <- p$data
  nodeCo <- intersect(df %>% filter(is.na(x)) %>%
                        select("parent", "node") %>% unlist(),
                      df %>% filter(!is.na(x)) %>%
                        select("parent", "node") %>% unlist())
  labCo <- df %>% filter("node" %in% nodeCo) %>%
    select("label") %>% unlist()
  selCo <- intersect(labCo, rownames(data))
  isSel <- "label" %in% selCo
  
  df <- df[df$isTip | isSel, ]
  start <- max(df$x, na.rm=TRUE) + offset
  
  dd <- as.data.frame(data)
  ## dd$lab <- rownames(dd)
  i <- order(df$y)
  
  ## handle collapsed tree
  ## https://github.com/GuangchuangYu/ggtree/issues/137
  i <- i[!is.na(df$y[i])]
  
  lab <- df$label[i]
  ## dd <- dd[lab, , drop=FALSE]
  ## https://github.com/GuangchuangYu/ggtree/issues/182
  dd <- dd[match(lab, rownames(dd)), , drop = FALSE]
  
  
  dd$y <- sort(df$y)
  dd$lab <- lab
  ## dd <- melt(dd, id=c("lab", "y"))
  dd <- gather(dd, variable, value, -c(lab, y))
  
  i <- which(dd$value == "")
  if (length(i) > 0) {
    dd$value[i] <- NA
  }
  if (is.null(colnames_level)) {
    dd$variable <- factor(dd$variable, levels=colnames(data))
  } else {
    dd$variable <- factor(dd$variable, levels=colnames_level)
  }
  V2 <- start + as.numeric(dd$variable) * width
  mapping <- data.frame(from=dd$variable, to=V2)
  mapping <- unique(mapping)
  
  dd$x <- V2
  dd$width <- width
  dd[[".panel"]] <- factor("Tree")
  height <- 1.00
  size <- 0.1
  if (is.null(color)) {
    p2 <- p + geom_tile(data=dd, aes(x, y, fill=value, color=value), color='white',lwd=1.5, width=width, height=height, size=size, inherit.aes=FALSE)
  } else {
    p2 <- p + geom_tile(data=dd, aes(x, y, fill=value, color=value), color='white',lwd=1.5, width=width, height=height, size=size, inherit.aes=FALSE)
  }
  if (is(dd$value,"numeric")) {
    p2 <- p2 + scale_fill_gradient(low=low, high=high, na.value=NA, name = legend_title) # "white")
    p2 <- p2 + scale_color_gradient(low=low, high=high, na.value=NA, name = legend_title) # "white")
  } else {
    p2 <- p2 + scale_fill_discrete(na.value=NA, name = legend_title) #"white")
    p2 <- p2 + scale_color_discrete(na.value=NA, name = legend_title) #"white")
  }
  
  if (colnames) {
    if (colnames_position == "bottom") {
      y <- 0
    } else {
      y <- max(p$data$y) + 1
    }
    mapping$y <- y
    mapping[[".panel"]] <- factor("Tree")
    p2 <- p2 + geom_text(data=mapping, aes(x=to, y = y, label=from), size=font.size, family=family, inherit.aes = FALSE, angle=colnames_angle, nudge_x=colnames_offset_x, nudge_y = colnames_offset_y, hjust=hjust)
  }
  
  legend_position=legend_position
  
  p2 <- p2 + theme(legend.position=legend_position)
  ## p2 <- p2 + guides(fill = guide_legend(override.aes = list(colour = NULL)))
  
  if (!colnames) {
    ## https://github.com/GuangchuangYu/ggtree/issues/204
    p2 <- p2 + scale_y_continuous(expand = c(0,0))
  }
  
  attr(p2, "mapping") <- mapping
  return(p2)
}

### buildabat tree plots
bat_dir = "C:/Users/shelb/OneDrive/Desktop/Chapter 2 Stuff/Simulation/script files/compiled_results_May2025"
newick = file.path(bat_dir, "buildabat_tree.txt")
facs = file.path(bat_dir, "buildabat_factorMeans2.csv")
class_path <- file.path(bat_dir, "buildabat_factorMeans2.csv")

plot_factor_tree("buildabat", newick, facs, tip_labels = T, layout="rectangular", extra_offset=0.4, heat_width=0.2, legend_position="right", height=40, width=40, class_path = NA, combined=T)

#width 400 and height 1200 for exporting
```

#### Phylogeny/Factor Means plots (code adapted from Hassler et al 2022) - Myotis plot
```{r}
###Hassler factor-phylogeny plotting code
library(ggtree)
library(tidyr)
library(phytools)
library(tidytree)
library(ggplot2)
library(aplot)
library(RColorBrewer)
library(ggnewscale)
library(phyclust)
library(geiger)
require(treeio)
library(extrafont)
library(showtext)

#set up the functions for phylogeny plots
prep_trait <- function(tree, trait){
  fit <- phytools::fastAnc(tree, trait, vars=TRUE, CI=TRUE)
  td <- data.frame(node = nodeid(tree, names(trait)),
                   trait = trait)
  nd <- data.frame(node = names(fit$ace), trait = fit$ace)
  d <- rbind(td, nd)
  d$node <- as.numeric(d$node)
  tree <- full_join(tree, d, by = 'node')
  return(tree)
}

plot_tree <- function(tree, colors,
                      border=FALSE,
                      line_width=0.1,
                      tip_labels=FALSE,
                      layout="rectangular",
                      color_tree=TRUE,
                      fan.angle=15,
                      new_labels=NA,
                      limits=NA,
                      labels_offset=0) {
  p <- ggtree(tree, layout=layout, open.angle=fan.angle, size = 0)
  if (border){
    # p <- p + geom_tree(size=line_width)
  } else {
    #   p <- p + ggtree(tree)
  }
  if (color_tree) {
    p <- p + geom_tree(aes(color=trait), continuous = TRUE, size=line_width)
    if (is.na(limits)) {
      #p <- p + scale_color_gradient2(midpoint = 0.0, low=colors[[1]], mid='grey', high=colors[[2]])
      p <- p + scale_color_gradient2(midpoint = 0.0, low=colors[[1]], mid='grey85', high=colors[[2]])
    } else {
      #p <- p + scale_color_gradient2(midpoint = 0.0, low=colors[[1]], mid='grey', high=colors[[2]], limits=limits)
      p <- p + scale_color_gradient2(midpoint = 0.0, low=colors[[1]], mid='grey85', high=colors[[2]], limits=limits)
      
    }
    
    p <- p + labs(color="Factor Value")
    
  } else {
    p <- p + geom_tree(size=line_width)
  }
  
  if (tip_labels) {
    if (!is.na(new_labels)) {
      p <- p %<+% new_labels + geom_tiplab(aes(label=classes, family="Bahnschrift", fontface="bold", size=18), offset=labels_offset)
    } else {
      p <- p + geom_tiplab(align=T, offset=labels_offset, geom="text", family="Bahnschrift", fontface="bold", size=4) #size=4, fontface=1
    }
  } else {
    # p <- p + geom_tiplab(aes(label=character(1)), align=TRUE)
  }
  
  # p <- ggtree(tree, size=2.8) +
  #   geom_tree(aes(color=trait),
  #             continuous = TRUE, size=2) +
  #   scale_color_gradient2(midpoint = 0, low='purple', mid='white', high='orange')
  return(p)
}

plot_factor_tree <- function(name, tree_path, factors_path, factors = NA,
                             class_path = NA,
                             height=NA,
                             width=NA,
                             border=FALSE,
                             line_width=1.0,
                             tip_labels=TRUE,
                             layout="rectangular",
                             class_palette=NA,
                             combined=TRUE,
                             fan.angle=30.0,
                             new_labels=NA,
                             common_scale=FALSE,
                             scale=TRUE,
                             heat_width=heat_width,
                             extra_offset=extra_offset,
                             legend_position=legend_position,
                             labels_offset=0.02,
                             fac_names=NA,
                             factor_fill = scale_fill_gradient2(midpoint = 0.0, low="#053439", mid='#F4F5F6', high="#440E33"), #low="#977294", mid='grey98', high="#F58551"
                             factor_color = scale_color_gradient2(midpoint = 0.0, low="#053439", mid='#F4F5F6', high="#440E33"), #low="#5D1E57", mid='#F4F5F6', high="#B33C05"
                             class_values=c("#729EA1","#B6D19E","#975E7C"), #69BADD #c("#6D554A","#06D6A0","#1B9AAA"), #c("#363537","#0CCE6B","#00A6FB"), #c("#519e8a","#7eb09b","#c5c9a4"),
                             legend_name="Ecotype",
                             relabel = NA,
                             include_only = NA
) {
  
  x <- as.matrix(read.csv(factors_path, header=TRUE))
  k <- ncol(x) - 1
  if (all(is.na(factors))) {
    factors <- 1:k
  }
  
  base_tree <- read.tree(tree_path)
  taxa <- x[, 1]
  
  include_class = !is.na(class_path)
  if (include_class) {
    class_df = read.csv(class_path)
  }
  
  if (!all(is.na(include_only))) {
    drop_taxa = setdiff(taxa, include_only)
    keep_rows = match(include_only, taxa)
    x <- x[keep_rows, ]
    for (taxon in drop_taxa) {
      base_tree <- drop.tip(base_tree, taxon)
    }
    
    if (include_class) {
      keep_rows <- match(include_only, class_df$taxon)
      class_df <- class_df[keep_rows, ]
    }
  }
  
  taxa <- x[, 1]
  
  if (!all(is.na(relabel))) {
    matched_rows <- match(relabel$original, taxa)
    taxa[matched_rows] <- relabel$new
    
    matched_rows <- match(relabel$original, base_tree$tip.label)
    base_tree$tip.label[matched_rows] <- relabel$new
    
    if (include_class) {
      matched_rows <- match(relabel$original, class_df$taxon)
      class_df$taxon[matched_rows] <- relabel$new
    }
    
  }
  
  x <- x[, factors + 1]
  x <- apply(as.matrix(x), 2, as.numeric)
  rownames(x) <- taxa
  n <- length(taxa)
  if (scale) {
    x <- scale(x)
  }
  
  limits <- NA
  if (common_scale) {
    x_min = min(x)
    x_max = max(x)
    limits = c(x_min, x_max)
  }
  
  combined_layout = layout
  if (is.na(height)) {
    height <- 20
    if (layout == "rectangular") {
      height <- 0.2 * n
    }
  }
  
  if (is.na(width)) {
    width <- 20
  }
  if (layout=="circular") {
    combined_layout = "fan"
  }
  
  # include_class = !is.na(class_path)
  if (include_class) {
    
    # stop()
    class_df = read.csv(class_path)
    class_df[,2] <- as.factor(class_df[,2])
    classes <- data.frame(x = class_df[,2])
    colnames(classes) <- c(colnames(class_df)[2])
    row.names(classes) <- class_df$taxon
    
    class_fills = scale_fill_manual(values=class_values, name=legend_name)
    #class_fills = scale_fill_manual(values=c("#519e8a","#7eb09b","#c5c9a4"), name="Ecotype")
    #class_fills = discrete_scale("fill", "manual", colorRampPalette(brewer.pal(min(nrow(unique(classes)), 12), "Set2")), name=colnames(classes)[1])
    
    class_colors = scale_color_manual(values=class_values, name=legend_name)
    #class_colors = scale_color_manual(values=c("#519e8a","#7eb09b","#c5c9a4"), name="Ecotype")
    #class_colors = discrete_scale("color", "manual", colorRampPalette(brewer.pal(min(nrow(unique(classes)), 12), "Set2")), name=colnames(classes)[1])
    
    if (!all(is.na(class_palette))) {
      #class_fills = scale_fill_manual(values=class_palette, name=colnames(classes)[1])
      #class_colors = scale_color_manual(values=class_palette, name=colnames(classes)[1])
    }
  }
  
  n_factors = length(factors)
  
  if (all(is.na(fac_names))) {
    fac_names <- character(n_factors)
    for (i in 1:n_factors) {
      fac_names[i] = paste("F", i)
    }
  }
  
  colnames(x) <- fac_names
  max_height <- max(node.depth.edgelength(base_tree))
  base_tree$edge.length <- base_tree$edge.length / max_height
  base_tree$edge.width <- base_tree$edge.width / (10*max_height)
  
  heat_width=heat_width

  
  if (combined) {
    legend_position=legend_position
    p1 <- plot_tree(base_tree, "", border=border, line_width=line_width, tip_labels = tip_labels, layout=combined_layout, color_tree=FALSE, new_labels=new_labels, labels_offset=labels_offset)
    pname <- paste(name, "_factors.svg", sep="")
    p2 <- p1 + theme(text = element_text(family="Arial", face="bold", size=18, colour="black"))
    if (include_class) {
      #p2 <- my_gheatmap(p1, classes, offset=extra_offset, width=heat_width, colnames_angle=90, colnames=TRUE, legend_title=colnames(classes)[1]) + class_fills #+ class_colors
      #p2 <- my_gheatmap(p1, classes, offset=extra_offset, width=heat_width, colnames_angle=90, colnames=TRUE, legend_position=legend_position, legend_title="Ecotype") + class_fills + class_colors
      p2 <- my_gheatmap(p1, classes, offset=extra_offset, width=heat_width, colnames_angle=90, colnames=F, legend_position=legend_position, legend_title="Ecotype") +
        class_fills + class_colors +
        theme(text = element_text(family="Arial", face="bold", size=18, colour="black"))
    }
    
    p3 <- p2 + new_scale_fill() + new_scale_color()
    p4 <- my_gheatmap(p3, x, offset=extra_offset + heat_width * 2, width=heat_width * n_factors, colnames_angle=90, colnames_offset_y = 0, legend_position=legend_position, legend_title="Factor Value") +
      factor_fill + factor_color +
      labs(fill="Factor Value") +
      guides(color=FALSE) +
      theme(text = element_text(family="Arial", face="bold", size=12, colour="black"))
    
    svg(pname, height=height, width=width)
    print(p4)
    dev.off()
    gc()
    return(p4)
    
  }
  gc()
}


my_gheatmap <- function(p, data, legend_position=legend_position, offset=0, width=2, low="green", high="red", color="white",
                        colnames=T, colnames_position="bottom", colnames_angle=0, colnames_level=NULL, #colnames=TRUE for heatmap labeled
                        colnames_offset_x = 0, colnames_offset_y = -2, font.size=4, family="", hjust=0.5, legend_title = "value", alpha=0.8) {
  
  colnames_position %<>% match.arg(c("bottom", "top"))
  variable <- value <- lab <- y <- NULL
  
  ## convert width to width of each cell
  width <- width * (p$data$x %>% range(na.rm=TRUE) %>% diff) / ncol(data)
  
  isTip <- x <- y <- variable <- value <- from <- to <- NULL
  
  df <- p$data
  nodeCo <- intersect(df %>% filter(is.na(x)) %>%
                        select("parent", "node") %>% unlist(),
                      df %>% filter(!is.na(x)) %>%
                        select("parent", "node") %>% unlist())
  labCo <- df %>% filter("node" %in% nodeCo) %>%
    select("label") %>% unlist()
  selCo <- intersect(labCo, rownames(data))
  isSel <- "label" %in% selCo
  
  df <- df[df$isTip | isSel, ]
  start <- max(df$x, na.rm=TRUE) + offset
  
  dd <- as.data.frame(data)
  ## dd$lab <- rownames(dd)
  i <- order(df$y)
  
  ## handle collapsed tree
  ## https://github.com/GuangchuangYu/ggtree/issues/137
  i <- i[!is.na(df$y[i])]
  
  lab <- df$label[i]
  ## dd <- dd[lab, , drop=FALSE]
  ## https://github.com/GuangchuangYu/ggtree/issues/182
  dd <- dd[match(lab, rownames(dd)), , drop = FALSE]
  
  
  dd$y <- sort(df$y)
  dd$lab <- lab
  ## dd <- melt(dd, id=c("lab", "y"))
  dd <- gather(dd, variable, value, -c(lab, y))
  
  i <- which(dd$value == "")
  if (length(i) > 0) {
    dd$value[i] <- NA
  }
  if (is.null(colnames_level)) {
    dd$variable <- factor(dd$variable, levels=colnames(data))
  } else {
    dd$variable <- factor(dd$variable, levels=colnames_level)
  }
  V2 <- start + as.numeric(dd$variable) * width
  mapping <- data.frame(from=dd$variable, to=V2)
  mapping <- unique(mapping)
  
  dd$x <- V2
  dd$width <- width
  dd[[".panel"]] <- factor("Tree")
  height <- 1.00
  size <- 0.1
  if (is.null(color)) {
    p2 <- p + geom_tile(data=dd, aes(x, y, fill=value, color=value), color='white',lwd=1.5, width=width, height=height, size=size, inherit.aes=FALSE, alpha=alpha)
  } else {
    p2 <- p + geom_tile(data=dd, aes(x, y, fill=value, color=value), color='white',lwd=1.5, width=width, height=height, size=size, inherit.aes=FALSE, alpha=alpha)
  }
  if (is(dd$value,"numeric")) {
    p2 <- p2 + scale_fill_gradient(low=low, high=high, na.value=NA, name = legend_title) # "white")
    p2 <- p2 + scale_color_gradient(low=low, high=high, na.value=NA, name = legend_title) # "white")
  } else {
    p2 <- p2 + scale_fill_discrete(na.value=NA, name = legend_title) #"white")
    p2 <- p2 + scale_color_discrete(na.value=NA, name = legend_title) #"white")
  }
  
  if (colnames) {
    if (colnames_position == "bottom") {
      y <- 0
    } else {
      y <- max(p$data$y) + 1
    }
    mapping$y <- y
    mapping[[".panel"]] <- factor("Tree")
    p2 <- p2 + geom_text(data=mapping, aes(x=to, y = y, label=from), size=font.size, family=family, inherit.aes = FALSE,
                         angle=colnames_angle, nudge_x=colnames_offset_x, nudge_y = colnames_offset_y, hjust=hjust)
  }
  
  legend_position=legend_position
  
  p2 <- p2 + theme(legend.position=legend_position)
  ## p2 <- p2 + guides(fill = guide_legend(override.aes = list(colour = NULL)))
  
  if (!colnames) {
    ## https://github.com/GuangchuangYu/ggtree/issues/204
    p2 <- p2 + scale_y_continuous(expand = c(0,0))
  }
  
  attr(p2, "mapping") <- mapping
  return(p2)
}

### myotis tree plots
bat_dir = "C:/Users/shelb/OneDrive/Desktop/Chapter 2 Stuff/Simulation/script files/compiled_results_May2025"
newick = file.path(bat_dir, "new_myotis_gunnell.txt")
facs = file.path(bat_dir, "myotisgunnell_factorMeans.csv")
class_path <- file.path(bat_dir, "myotis_classes.csv")

plot_factor_tree("myotis", newick, facs, tip_labels = T, layout="rectangular", extra_offset=1.5, heat_width=0.2, legend_position="right", height=120, width=40, class_path = class_path, combined=T)

#575 width & 1200 height for exporting the base figure image
```

#### Phylogeny/Factor Means plots (code adapted from Hassler et al 2022) - Phyllostomidae plot
```{r}
###Hassler factor-phylogeny plotting code
library(ggtree)
library(tidyr)
library(phytools)
library(tidytree)
library(ggplot2)
library(aplot)
library(RColorBrewer)
library(ggnewscale)
library(phyclust)
library(geiger)
require(treeio)
library(extrafont)
library(showtext)

#load in extra font
font_add(family = "Bahnschrift", regular = "Bahnschrift.ttf")
showtext_auto()

#set up functions for phylogeny plots
prep_trait <- function(tree, trait){
  fit <- phytools::fastAnc(tree, trait, vars=TRUE, CI=TRUE)
  td <- data.frame(node = nodeid(tree, names(trait)),
                   trait = trait)
  nd <- data.frame(node = names(fit$ace), trait = fit$ace)
  d <- rbind(td, nd)
  d$node <- as.numeric(d$node)
  tree <- full_join(tree, d, by = 'node')
  return(tree)
}

plot_tree <- function(tree, colors,
                      border=FALSE,
                      line_width=0.1,
                      tip_labels=FALSE,
                      layout="rectangular",
                      color_tree=TRUE,
                      fan.angle=15,
                      new_labels=NA,
                      limits=NA,
                      labels_offset=0) {
  p <- ggtree(tree, layout=layout, open.angle=fan.angle, size = 0)
  if (border){
    # p <- p + geom_tree(size=line_width)
  } else {
    #   p <- p + ggtree(tree)
  }
  if (color_tree) {
    p <- p + geom_tree(aes(color=trait),
                       continuous = TRUE, size=line_width)
    if (is.na(limits)) {
      #p <- p + scale_color_gradient2(midpoint = 0.0, low=colors[[1]], mid='grey', high=colors[[2]])
      p <- p + scale_color_gradient2(midpoint = 0.0, low=colors[[1]], mid='grey85', high=colors[[2]])
    } else {
      #p <- p + scale_color_gradient2(midpoint = 0.0, low=colors[[1]], mid='grey', high=colors[[2]], limits=limits)
      p <- p + scale_color_gradient2(midpoint = 0.0, low=colors[[1]], mid='grey85', high=colors[[2]], limits=limits)
      
    }
    
    p <- p + labs(color="Factor Value")
    
  } else {
    p <- p + geom_tree(size=line_width)
  }
  
  if (tip_labels) {
    if (!is.na(new_labels)) {
      p <- p %<+% new_labels + geom_tiplab(aes(label=classes, family="Bahnschrift", fontface="bold", size=18), offset=labels_offset)
    } else {
      p <- p + geom_tiplab(align=T, offset=labels_offset, geom="text", family="Bahnschrift", fontface="bold", size=4) #size=4, fontface=1
    }
  } else {
    # p <- p + geom_tiplab(aes(label=character(1)), align=TRUE)
  }
  
  # p <- ggtree(tree, size=2.8) +
  #   geom_tree(aes(color=trait),
  #             continuous = TRUE, size=2) +
  #   scale_color_gradient2(midpoint = 0, low='purple', mid='white', high='orange')
  return(p)
}

plot_factor_tree <- function(name, tree_path, factors_path, factors = NA,
                             class_path = NA,
                             height=NA,
                             width=NA,
                             border=FALSE,
                             line_width=1.0,
                             tip_labels=TRUE,
                             layout="rectangular",
                             class_palette=NA,
                             combined=TRUE,
                             fan.angle=30.0,
                             new_labels=NA,
                             common_scale=FALSE,
                             scale=TRUE,
                             heat_width=heat_width,
                             extra_offset=extra_offset,
                             legend_position=legend_position,
                             labels_offset=0.02,
                             fac_names=NA,
                             factor_fill = scale_fill_gradient2(midpoint = 0.0, low="#053439", mid='#F4F5F6', high="#440E33"), #low="#977294", mid='grey98', high="#F58551"
                             factor_color = scale_color_gradient2(midpoint = 0.0, low="#053439", mid='#F4F5F6', high="#440E33"),
                             class_values=c("#729EA1","#B6D19E","#975E7C"),
                             legend_name="Ecotype",
                             relabel = NA,
                             include_only = NA
) {
  
  x <- as.matrix(read.csv(factors_path, header=TRUE))
  k <- ncol(x) - 1
  if (all(is.na(factors))) {
    factors <- 1:k
  }
  
  base_tree <- read.tree(tree_path)
  taxa <- x[, 1]
  
  include_class = !is.na(class_path)
  if (include_class) {
    class_df = read.csv(class_path)
  }
  
  if (!all(is.na(include_only))) {
    drop_taxa = setdiff(taxa, include_only)
    keep_rows = match(include_only, taxa)
    x <- x[keep_rows, ]
    for (taxon in drop_taxa) {
      base_tree <- drop.tip(base_tree, taxon)
    }
    
    if (include_class) {
      keep_rows <- match(include_only, class_df$taxon)
      class_df <- class_df[keep_rows, ]
    }
  }
  
  taxa <- x[, 1]
  
  if (!all(is.na(relabel))) {
    matched_rows <- match(relabel$original, taxa)
    taxa[matched_rows] <- relabel$new
    
    matched_rows <- match(relabel$original, base_tree$tip.label)
    base_tree$tip.label[matched_rows] <- relabel$new
    
    if (include_class) {
      matched_rows <- match(relabel$original, class_df$taxon)
      class_df$taxon[matched_rows] <- relabel$new
    }
    
  }
  
  x <- x[, factors + 1]
  x <- apply(as.matrix(x), 2, as.numeric)
  rownames(x) <- taxa
  n <- length(taxa)
  if (scale) {
    x <- scale(x)
  }
  
  limits <- NA
  if (common_scale) {
    x_min = min(x)
    x_max = max(x)
    limits = c(x_min, x_max)
  }
  
  combined_layout = layout
  if (is.na(height)) {
    height <- 20
    if (layout == "rectangular") {
      height <- 0.2 * n
    }
  }
  
  if (is.na(width)) {
    width <- 20
  }
  if (layout=="circular") {
    combined_layout = "fan"
  }
  
  # include_class = !is.na(class_path)
  if (include_class) {
    
    # stop()
    class_df = read.csv(class_path)
    class_df[,2] <- as.factor(class_df[,2])
    classes <- data.frame(x = class_df[,2])
    colnames(classes) <- c(colnames(class_df)[2])
    row.names(classes) <- class_df$taxon
    
    class_fills = scale_fill_manual(values=class_values, name=legend_name)
    class_colors = scale_color_manual(values=class_values, name=legend_name)
    
    if (!all(is.na(class_palette))) {
      #class_fills = scale_fill_manual(values=class_palette, name=colnames(classes)[1])
      #class_colors = scale_color_manual(values=class_palette, name=colnames(classes)[1])
    }
  }
  
  n_factors = length(factors)
  
  if (all(is.na(fac_names))) {
    fac_names <- character(n_factors)
    for (i in 1:n_factors) {
      fac_names[i] = paste("F", i)
    }
  }
  
  colnames(x) <- fac_names
  max_height <- max(node.depth.edgelength(base_tree))
  base_tree$edge.length <- base_tree$edge.length / max_height
  base_tree$edge.width <- base_tree$edge.width / (10*max_height)
  
  heat_width=heat_width

  
  if (combined) {
    legend_position=legend_position
    p1 <- plot_tree(base_tree, "", border=border, line_width=line_width, tip_labels = tip_labels, layout=combined_layout, color_tree=FALSE, new_labels=new_labels, labels_offset=labels_offset)
    pname <- paste(name, "_factors.svg", sep="")
    p2 <- p1 + theme(text = element_text(family="Bahnschrift", face="bold", size=18, colour="black"))
    if (include_class) {
      #p2 <- my_gheatmap(p1, classes, offset=extra_offset, width=heat_width, colnames_angle=90, colnames=TRUE, legend_title=colnames(classes)[1]) + class_fills #+ class_colors
      p2 <- my_gheatmap(p1, classes, offset=extra_offset, width=heat_width, colnames_angle=90, colnames=TRUE, legend_position=legend_position, legend_title="Ecotype") + 
        class_fills + class_colors +
        theme(text = element_text(family="Bahnschrift", face="bold", size=18, colour="black"))
    }
    
    p3 <- p2 + new_scale_fill() + new_scale_color()
    p4 <- my_gheatmap(p3, x, offset=extra_offset + heat_width * 2, width=heat_width * n_factors, colnames_angle=90, colnames_offset_y = 0, legend_position=legend_position, legend_title="Factor Value") +
      factor_fill + factor_color +
      labs(fill="Factor Value") +
      guides(color=FALSE) +
      theme(text = element_text(family="Bahnschrift", face="bold", size=18, colour="black"))
    
    svg(pname, height=height, width=width)
    print(p4)
    dev.off()
    gc()
    return(p4)
    
  }
  gc()
}


my_gheatmap <- function(p, data, legend_position=legend_position, offset=0, width=1, low="green", high="red", color="white",
                        colnames=T, colnames_position="bottom", colnames_angle=0, colnames_level=NULL, #colnames=TRUE for heatmap labeled
                        colnames_offset_x = 0, colnames_offset_y = 0, font.size=4, family="", hjust=0.5, legend_title = "value", alpha=0.8) {
  
  colnames_position %<>% match.arg(c("bottom", "top"))
  variable <- value <- lab <- y <- NULL
  
  ## convert width to width of each cell
  width <- width * (p$data$x %>% range(na.rm=TRUE) %>% diff) / ncol(data)
  
  isTip <- x <- y <- variable <- value <- from <- to <- NULL
  
  df <- p$data
  nodeCo <- intersect(df %>% filter(is.na(x)) %>%
                        select("parent", "node") %>% unlist(),
                      df %>% filter(!is.na(x)) %>%
                        select("parent", "node") %>% unlist())
  labCo <- df %>% filter("node" %in% nodeCo) %>%
    select("label") %>% unlist()
  selCo <- intersect(labCo, rownames(data))
  isSel <- "label" %in% selCo
  
  df <- df[df$isTip | isSel, ]
  start <- max(df$x, na.rm=TRUE) + offset
  
  dd <- as.data.frame(data)
  ## dd$lab <- rownames(dd)
  i <- order(df$y)
  
  ## handle collapsed tree
  ## https://github.com/GuangchuangYu/ggtree/issues/137
  i <- i[!is.na(df$y[i])]
  
  lab <- df$label[i]
  ## dd <- dd[lab, , drop=FALSE]
  ## https://github.com/GuangchuangYu/ggtree/issues/182
  dd <- dd[match(lab, rownames(dd)), , drop = FALSE]
  
  
  dd$y <- sort(df$y)
  dd$lab <- lab
  ## dd <- melt(dd, id=c("lab", "y"))
  dd <- gather(dd, variable, value, -c(lab, y))
  
  i <- which(dd$value == "")
  if (length(i) > 0) {
    dd$value[i] <- NA
  }
  if (is.null(colnames_level)) {
    dd$variable <- factor(dd$variable, levels=colnames(data))
  } else {
    dd$variable <- factor(dd$variable, levels=colnames_level)
  }
  V2 <- start + as.numeric(dd$variable) * width
  mapping <- data.frame(from=dd$variable, to=V2)
  mapping <- unique(mapping)
  
  dd$x <- V2
  dd$width <- width
  dd[[".panel"]] <- factor("Tree")
  height <- 1.00
  size <- 0.1
  if (is.null(color)) {
    p2 <- p + geom_tile(data=dd, aes(x, y, fill=value, color=value), color='white',lwd=1.5, width=width, height=height, size=size, inherit.aes=FALSE, alpha=alpha)
  } else {
    p2 <- p + geom_tile(data=dd, aes(x, y, fill=value, color=value), color='white',lwd=1.5, width=width, height=height, size=size, inherit.aes=FALSE, alpha=alpha)
  }
  if (is(dd$value,"numeric")) {
    p2 <- p2 + scale_fill_gradient(low=low, high=high, na.value=NA, name = legend_title) # "white")
    p2 <- p2 + scale_color_gradient(low=low, high=high, na.value=NA, name = legend_title) # "white")
  } else {
    p2 <- p2 + scale_fill_discrete(na.value=NA, name = legend_title) #"white")
    p2 <- p2 + scale_color_discrete(na.value=NA, name = legend_title) #"white")
  }
  
  if (colnames) {
    if (colnames_position == "bottom") {
      y <- 0
    } else {
      y <- max(p$data$y) + 1
    }
    mapping$y <- y
    mapping[[".panel"]] <- factor("Tree")
    p2 <- p2 + geom_text(data=mapping, aes(x=to, y = y, label=from), size=font.size, family=family, inherit.aes = FALSE,
                         angle=colnames_angle, nudge_x=colnames_offset_x, nudge_y = colnames_offset_y, hjust=hjust)
  }
  
  legend_position=legend_position
  
  p2 <- p2 + theme(legend.position=legend_position)
  ## p2 <- p2 + guides(fill = guide_legend(override.aes = list(colour = NULL)))
  
  if (!colnames) {
    ## https://github.com/GuangchuangYu/ggtree/issues/204
    p2 <- p2 + scale_y_continuous(expand = c(0,0))
  }
  
  attr(p2, "mapping") <- mapping
  return(p2)
}

### phyllostomidae tree plots
bat_dir = "C:/Users/shelb/OneDrive/Desktop/Chapter 2 Stuff/Simulation/script files/compiled_results_May2025"
newick = file.path(bat_dir, "phyllo_tree_upham.txt")
facs = file.path(bat_dir, "phylloupham_factorMeans.csv")
class_path <- file.path(bat_dir, "phylloupham_classes.csv")

plot_factor_tree("phyllo", newick, facs, tip_labels = T, layout="rectangular", extra_offset=2.4, heat_width=0.3, legend_position="right",  height=40, width=40, class_path = NA, combined=T)

#width 550 height 1500 for export
```
