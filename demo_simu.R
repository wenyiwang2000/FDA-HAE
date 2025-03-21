rm(list=ls())
library(mgcv)
library(splines)
library(face)
library(MASS)
library(rTensor)


# Load functions
source('./R_code/GeneData.R')
source("./R_code/FDA_HAE_algorithm.R")


# Generate data 
load("initial_data_gen.RData")
set.seed(1)
data <- GeneData(n=713,subj=initial_data_gen$subj, time=initial_data_gen$time, group=initial_data_gen$group, 
                 position=initial_data_gen$position, Nij=initial_data_gen$Nij, R=2,L=3, X=initial_data_gen$X, 
                 lambda=initial_data_gen$lambda0, A=initial_data_gen$A,U1=initial_data_gen$U1,
                 U2=initial_data_gen$U2,U3=initial_data_gen$U3, phi=initial_data_gen$phi,
                 tBmatrix_theta=initial_data_gen$tBmatrix_theta,miss_rate=0)

#Get estimations
out <- FDA_HAE_algorithm(R=2,L=3,nm=7,p=3,iters=1000, data1=data,Miniter=20)









