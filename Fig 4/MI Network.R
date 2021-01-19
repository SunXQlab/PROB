##### Initial MI network

## Breast cancer
PPI_pair=read.csv("F:/Clinical Gene expression network Project/Code/TCGA Data/Breast cancer/PPI_String.csv",header = FALSE)
NetGene12=PPI_pair[-1,1:2]
NetGene=Gene_Selected_Name
Netsize=length(NetGene)  #72

PPI=matrix(0,nrow=Netsize,ncol=Netsize)
for (i in 1:Netsize)
{
  for (j in 1:Netsize)
  {
    for (k in 1:dim(NetGene12)[1])
    {
      #if (is.na(prod(which(as.matrix(PPI_pair)[k,]==c(names_SR[i],names_SR[j]))))==0)
      if (sum(c(NetGene[i],NetGene[j])==as.matrix(NetGene12)[k,])==2)
        PPI[i,j]=PPI[i,j]+1
      else
        PPI[i,j]=PPI[i,j]+0  
      
    }
    
  }  
}


for (i in 1:Netsize)
{
  for (j in 1:Netsize)
  {
    PPI[i,j]=max(PPI[i,j],PPI[j,i])
  }
}

rownames(PPI)=Gene_Selected_Name
colnames(PPI)=Gene_Selected_Name
sum(PPI)  #=  edges

setwd("F:/Clinical Gene expression network Project/Code/TCGA Data/Breast cancer/")  
write.table(PPI, file="PPI matrix.xls",sep='\t')

## Melanoma 

PPI_pair=read.csv("F:/Clinical Gene expression network Project/Code/TCGA Data/Melanoma/PPI_String.csv",header = FALSE)
NetGene12=PPI_pair[-1,1:2]
NetGene=Gene_Selected_Name
Netsize=length(NetGene)  #72

PPI=matrix(0,nrow=Netsize,ncol=Netsize)
for (i in 1:Netsize)
{
  for (j in 1:Netsize)
  {
    for (k in 1:dim(NetGene12)[1])
    {
      #if (is.na(prod(which(as.matrix(PPI_pair)[k,]==c(names_SR[i],names_SR[j]))))==0)
      if (sum(c(NetGene[i],NetGene[j])==as.matrix(NetGene12)[k,])==2)
        PPI[i,j]=PPI[i,j]+1
      else
        PPI[i,j]=PPI[i,j]+0  
      
    }
    
  }  
}


for (i in 1:Netsize)
{
  for (j in 1:Netsize)
  {
    PPI[i,j]=max(PPI[i,j],PPI[j,i])
  }
}

rownames(PPI)=Gene_Selected_Name
colnames(PPI)=Gene_Selected_Name
sum(PPI)  #=  edges

setwd("F:/Clinical Gene expression network Project/Code/TCGA Data/Melanoma/")  
write.table(PPI, file="PPI matrix.xls",sep='\t')



