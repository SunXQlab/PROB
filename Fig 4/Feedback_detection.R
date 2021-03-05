## Detect feedback loops in the GRN of EMT of UC to SARC


setwd("F:/胶质瘤诱导分化组学数据分析/Code-cluster_based/Results") 


############# Read sensitive network coefficient
B_S=read.table("Significant Coefficient_Sensitive Network.xls",header = F,sep="\t",fileEncoding = "UTF-8")
B_Net=B_S
rownames(B_Net)=B_Net[,1]
# colnames(B_Net)
B_Net=B_Net[,-1]   #53X53

# ############# Read resistant network coefficient
# A=read.table("Coefficient_Resistant Network.xls",header = T,sep="\t",fileEncoding = "UTF-8",row.names=NULL)
# 
# B=A[,-c(1,2)]
# B=as.matrix(B)
# dim(B)
# B[abs(B)<=0.01]=0
# sum(sum(B!=0))   # 100
# 
# rn=rownames(edge_R[rowSums(abs(edge_R))!=0, ])
# cn=colnames(edge_R[,colSums(abs(edge_R))!=0])
# Net_Gene=union(rn,cn)
# 
# rownames(B)=Net_Gene
# colnames(B)=Net_Gene
# 
# 
# 
# # B=read.table("Significant Coefficient_Resistant Network.xls",header = T,sep="\t")
# edge_R=read.table("edge matrix for Resistant genes.txt",header = T,sep="\t")
# rownames(edge_R)=edge_R[,1]
# edge_R=edge_R[,-1]
# 
# ### Select genes in the network
# a=rownames(edge_R[rowSums(abs(B))!=0, ])
# b=colnames(edge_R[,colSums(abs(B))!=0])
# DirectedNet_Gene=union(a,b)
# 
# B_Net=B[DirectedNet_Gene,DirectedNet_Gene]
# dim(B_Net)  # 44,44
# 
# write.table(B_Net, file="Coefficient_GeneNetwork_Resistant.txt",col.names = NA,sep = "\t")

############
B_R=read.table("Significant Coefficient_Resistant Network.xls",header = F,sep="\t",fileEncoding = "UTF-8")
B_Net=B_R
rownames(B_Net)=B_Net[,1]
# colnames(B_Net)
B_Net=B_Net[,-1]   #44X44
##############


###########  FB detection ################
N=dim(B_Net)[1]
NFB=NULL

for (i in 1:N)
{
  for (j in i+1:N)
  {
    if (j<=N)
    {
      print(i)
      print(j)
      if ((B_Net[i,j]<0 & B_Net[j,i]>0) || (B_Net[i,j]>0 & B_Net[j,i]<0))
      {
        NFB=rbind(NFB,cbind(rownames(B_Net)[i],rownames(B_Net)[j],B_Net[i,j],B_Net[j,i]))
      }
    }
    
  }
}

NFB


N=dim(B_Net)[1]

PFB=NULL

for (i in 1:N)
{
  for (j in i+1:N)
  {
    if (j<=N)
    {
      print(i)
      print(j)
      if ((B_Net[i,j]>0 & B_Net[j,i]>0))
      {
        PFB=rbind(PFB,cbind(rownames(B_Net)[i],rownames(B_Net)[j],B_Net[i,j],B_Net[j,i]))
      }
    }
    
  }
}

PFB


N=dim(B_Net)[1]

PFB=NULL

for (i in 1:N)
{
  for (j in i+1:N)
  {
    if (j<=N)
    {
      print(i)
      print(j)
      if ((B_Net[i,j]>0 & B_Net[j,i]>0) || (B_Net[i,j]<0 & B_Net[j,i]<0))
      {
        PFB=rbind(PFB,cbind(rownames(B_Net)[i],rownames(B_Net)[j],B_Net[i,j],B_Net[j,i]))
      }
    }
    
  }
}

PFB



PNFB=NULL  # coupled FB
M1=dim(PFB)[1]
M2=dim(NFB)[1]
for (k in 1:M1)
{
  for (m in 1:M2)
  {
    if (k<=M1 & m<=M2)
    {
      print(k)
      print(m)
      if (length(intersect(PFB[k,1:2],NFB[m,1:2]))!=0)
        PNFB=rbind(PNFB,c(PFB[k,],NFB[m,]))
    }
    
  }
}

PNFB

PPFB=NULL  # coupled FB
M1=dim(PFB)[1]
M2=dim(PFB)[1]
for (k in 1:M1)
{
  for (m in 1:M2)
  {
    if (k<=M1 & m<=M2)
    {
      print(k)
      print(m)
      if (length(intersect(PFB[k,1:2],PFB[m,1:2]))!=0)
        PPFB=rbind(PPFB,c(PFB[k,],PFB[m,]))
    }
    
  }
}

NNFB


NNFB=NULL  # coupled FB
M1=dim(NFB)[1]
M2=dim(NFB)[1]
for (k in 1:M1)
{
  for (m in 1:M2)
  {
    if (k<=M1 & m<=M2)
    {
      print(k)
      print(m)
      if (length(intersect(NFB[k,1:2],NFB[m,1:2]))!=0)
        NNFB=rbind(NNFB,c(NFB[k,],NFB[m,]))
    }
    
  }
}

NNFB


n1=dim(PFB)[1]
n2=dim(NFB)[1]
n3=dim(PNFB)[1]
n4=dim(PPFB)[1]
n5=dim(NNFB)[1]