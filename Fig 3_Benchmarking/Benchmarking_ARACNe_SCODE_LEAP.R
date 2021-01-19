install.packages("R.matlab")
library(R.matlab)
setwd("F:/Clinical Gene expression network Project/Reversion/Codes/Benchmarking_LPS_scRNAseq_data/")
X_stage <- readMat("scExp_Net_Ctime.mat")
X_stage=as.data.frame(X_stage)
Expression=X_stage[1:(dim(X_stage)[1]-1),]
Time_Stamp=X_stage[dim(X_stage)[1],]

Timepoint=readMat("DPP.mat")
Timepoint=as.numeric(unlist(Timepoint))

### ARACNe
# install.packages("minet")
# manually installed

library(minet)

mim <- build.mim(t(Expression),estimator="spearman")
net1 <- aracne(mim)
write.table(net1, file="F:/Clinical Gene expression network Project/Reversion/Codes/Benchmarking_LPS_scRNAseq_data/Results/net_ARACNe.csv",sep=',') 

net2 <- clr(mim)
write.table(net2, file="F:/Clinical Gene expression network Project/Reversion/Codes/Benchmarking_LPS_scRNAseq_data/Results/net_CLR.csv",sep=',') 

net3 <- mrnet(mim)
write.table(net3, file="F:/Clinical Gene expression network Project/Reversion/Codes/Benchmarking_LPS_scRNAseq_data/Results/net_mrnet.csv",sep=',') 

#### SCODE
library(MASS)

fdata <- Expression
ftime <- Timepoint
dir <- "F:/Clinical Gene expression network Project/Reversion/Codes/Benchmarking_LPS_scRNAseq_data/Results/"
tfnum <- dim(X_stage)[1]-1
pnum <- 4 # default
cnum <- dim(X_stage)[2]
maxite <- 100 # default

maxB <- 2.0
minB <- -10.0

system(paste("mkdir", dir, sep=" "))

# X <- as.matrix(read.table(fdata, sep="\t"))[1:tfnum,1:cnum]
X <- (Expression)
W <- matrix(rep(0,tfnum*pnum), nrow=tfnum, ncol=pnum)
Z <- matrix(rep(0,pnum*cnum), nrow=pnum, ncol=cnum)
WZ <- matrix(nrow=tfnum, ncol=cnum)

#read pseudo-time and normalize pseudo-time
# pseudotime <- read.table(ftime, sep="\t")[1:cnum,2]
pseudotime <- Timepoint
pseudotime <- pseudotime/max(pseudotime)

new_B <- rep(0, pnum)
old_B <- rep(0, pnum)

#initialization
RSS <- Inf
for(i in 1:pnum){
  new_B[i] <- runif(1, min=minB, max=maxB)
  old_B[i] <- new_B[i]
}

#function to sample Z
sample_Z <- function(){
  for(i in 1:pnum){
    for(j in 1:cnum){
      Z[i,j] <<- exp(new_B[i]*pseudotime[j]) + runif(1, min=-0.001, max=0.001)
    }
  }
}

#optimize W and B iteratively
for(ite in 1:maxite){
  #sampling B
  target <- floor(runif(1, min=1, max=pnum+1))
  new_B[target] <- runif(1, min=minB, max=maxB)
  
  #for last calculation
  if(ite == maxite){
    for(i in 1:pnum){
      new_B[i] <- old_B[i]
    }
  }
  
  #sample Z from new B
  sample_Z()
  
  #regression
  for(i in 1:tfnum){
    X.lm <- lm(as.numeric(X[i,]) ~ t(Z)-1)
    for(j in 1:pnum){
      W[i,j] <- X.lm$coefficients[j]
    }
    WZ[i,] <- W[i,] %*% Z
  }
  
  #RSS
  tmp_RSS <- sum((X-WZ)**2)
  if(tmp_RSS < RSS){
    RSS <- tmp_RSS
    old_B[target] <- new_B[target]
  }
  else{
    new_B[target] <- old_B[target]
  }
}

# #output RSS
# write.table(RSS, paste(dir,"/RSS.txt",sep=""), row.names=F, col.names=F, sep="\t")
# 
# #output W
# write.table(W, paste(dir,"/W.txt",sep=""), row.names=F, col.names=F, sep="\t")

#infer A
B <- matrix(rep(0,pnum*pnum), nrow=pnum, ncol=pnum)
for(i in 1:pnum){
  B[i,i] <- new_B[i]
}
invW <- ginv(W)
A <- W %*% B %*% invW

#output A and B
write.csv(A, paste(dir,"/net_SCODE.csv",sep=""))

#### LEAP
# install.packages('LEAP')
library(LEAP)

MAC_result=MAC_counter(Expression[,order(Timepoint)], max_lag_prop = 1/3, MAC_cutoff = 0.2, 
            file_name = F, lag_matrix = T, symmetric = T)


Netsize=dim(Expression)[1]  #48 

MAC_matrix=matrix(0,nrow=Netsize,ncol=Netsize)
for (i in 1:Netsize)
{
  MAC_matrix[MAC_result[i,3],MAC_result[i,4]]=MAC_result[i,1]
}

write.csv(MAC_matrix, paste(dir,"/net_LEAP.csv",sep=""))

### SINCERITIES
## SINCERITIES needs at least 5 time points so it can not be used for clinical data 