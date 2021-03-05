### UC SARC

setwd("F:/Clinical Gene expression network Project/Reversion/Codes/UC_SARC")  

AM_UC=read.csv('Results/AdjacentMatrix_UC.csv',header=F)
AM_SARC=read.csv('Results/AdjacentMatrix_SARC.csv',header=F)


NodeID=unique(unlist(EMT_reggene_names))  #EMT_reggene_names
###  Reform to txt for Cytoscape

## UC 
C=c()
AM=matrix(as.numeric(unlist(AM_UC)),dim(AM_UC))
C=matrix(NA,nrow=sum(sum(AM!=0)),ncol=3)
row=1
for (i in 1:length(NodeID))
{
  for (j in 1:length(NodeID))
  {
    if (AM[j,i]!=0)   # 注意这里是AM[j,i]不是AM[i,j]
    {
      C[row,1]=NodeID[i]
      C[row,3]=NodeID[j]
      C[row,2]=AM[j,i]  # 注意这里是AM[j,i]不是AM[i,j]
      row=row+1
      # C=rbind(C,c(NodeID[i],AM[j,i],NodeID[j]))
    }
  }
}

C_UC=C

## SARC 
AM=matrix(as.numeric(unlist(AM_SARC)),dim(AM_SARC))
C=c()
C=matrix(NA,nrow=sum(sum(AM!=0)),ncol=3)
row=1
for (i in 1:length(NodeID))
{
  for (j in 1:length(NodeID))
  {
    if (AM[j,i]!=0)   # 注意这里是AM[j,i]不是AM[i,j]
    {
      C[row,1]=NodeID[i]
      C[row,3]=NodeID[j]
      C[row,2]=AM[j,i]  # 注意这里是AM[j,i]不是AM[i,j]
      row=row+1
      # C=rbind(C,c(NodeID[i],AM[j,i],NodeID[j]))
    }
  }
}

C_SARC=C

write.table(C_UC, file="Results/Cytoscape_UC_Nerwork.txt",quote=F,row.names=F,col.names = F,sep = "\t")
write.table(C_SARC, file="Results/Cytoscape_SARC_Nerwork.txt",quote=F,row.names=F,col.names = F,sep = "\t")


## Differential network
AM1=matrix(as.numeric(unlist(AM_UC)),dim(AM_UC))
AM2=matrix(as.numeric(unlist(AM_SARC)),dim(AM_SARC))
 AM_D=AM2-AM1



AM=AM_D
C=c()
C=matrix(NA,nrow=sum(sum(AM!=0)),ncol=3)
row=1
for (i in 1:length(NodeID))
{
  for (j in 1:length(NodeID))
  {
    if (AM[j,i]!=0)   # 注意这里是AM[j,i]不是AM[i,j]
    {
      C[row,1]=NodeID[i]
      C[row,3]=NodeID[j]
      C[row,2]=AM[j,i]  # 注意这里是AM[j,i]不是AM[i,j]
      row=row+1
      # C=rbind(C,c(NodeID[i],AM[j,i],NodeID[j]))
    }
  }
}

C_Diff=C

write.table(C_Diff, file="Results/Cytoscape_Diff_Network.txt",quote=F,row.names=F,col.names = F,sep = "\t")

## Rewired network
AM_D=AM2*0
for (i in 1:length(NodeID))
{
  for (j in 1:length(NodeID))
  {
    if (AM2[i,j]!=0 & AM1[i,j]==0 )  # & abs(AM2[i,j])>=0.01
    {AM_D[i,j]=AM2[i,j]}
  }
}

AM=AM_D
C=c()
C=matrix(NA,nrow=sum(sum(AM!=0)),ncol=3)
row=1
for (i in 1:length(NodeID))
{
  for (j in 1:length(NodeID))
  {
    if (AM[j,i]!=0)   # 注意这里是AM[j,i]不是AM[i,j]
    {
      C[row,1]=NodeID[i]
      C[row,3]=NodeID[j]
      C[row,2]=AM[j,i]  # 注意这里是AM[j,i]不是AM[i,j]
      row=row+1
      # C=rbind(C,c(NodeID[i],AM[j,i],NodeID[j]))
    }
  }
}

C_gain=C
write.table(C_gain, file="Results/Cytoscape_SARC_Unique_Network.txt",quote=F,row.names=F,col.names = F,sep = "\t")

# lost edges
AM_D=AM2*0
for (i in 1:length(NodeID))
{
  for (j in 1:length(NodeID))
  {
    if (AM1[i,j]!=0 & AM2[i,j]==0 )  # & abs(AM2[i,j])>=0.01
    {AM_D[i,j]=AM1[i,j]}
  }
}

AM=AM_D
C=c()
C=matrix(NA,nrow=sum(sum(AM!=0)),ncol=3)
row=1
for (i in 1:length(NodeID))
{
  for (j in 1:length(NodeID))
  {
    if (AM[j,i]!=0)   # 注意这里是AM[j,i]不是AM[i,j]
    {
      C[row,1]=NodeID[i]
      C[row,3]=NodeID[j]
      C[row,2]=AM[j,i]  # 注意这里是AM[j,i]不是AM[i,j]
      row=row+1
      # C=rbind(C,c(NodeID[i],AM[j,i],NodeID[j]))
    }
  }
}

C_lost=C
write.table(C_lost, file="Results/Cytoscape_UC_Unique_Network.txt",quote=F,row.names=F,col.names = F,sep = "\t")


## UC_SARC whole network
AM_UC_SARC=read.csv('Results/AdjacentMatrix_UC_SARC.csv',header=F)
AM=matrix(as.numeric(unlist(AM_UC_SARC)),dim(AM_UC_SARC))
C=c()
C=matrix(NA,nrow=sum(sum(AM!=0)),ncol=3)
row=1
for (i in 1:length(NodeID))
{
  for (j in 1:length(NodeID))
  {
    if (AM[j,i]!=0)   # 注意这里是AM[j,i]不是AM[i,j]
    {
      C[row,1]=NodeID[i]
      C[row,3]=NodeID[j]
      C[row,2]=AM[j,i]  # 注意这里是AM[j,i]不是AM[i,j]
      row=row+1
      # C=rbind(C,c(NodeID[i],AM[j,i],NodeID[j]))
    }
  }
}

C_UC_SARC=C
write.table(C_UC_SARC, file="Results/Cytoscape_UC_SARC_Network.txt",quote=F,row.names=F,col.names = F,sep = "\t")

