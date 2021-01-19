### Breast cancer  GSE7390

setwd("F:/Clinical Gene expression network Project/Manuscript/Codes/More cancers/Breast cancer GEO/")  

AM=read.csv('Results/AdjacentMatrix.csv',header=F)

ind_selected=read.csv("Results/ind_selected.csv",header = FALSE)
ind_selected=as.numeric((ind_selected))

Gene_Expression=read.csv("Data/Gene_GSE7390.csv")
GeneSymbol=as.character(unlist(Gene_Expression[,1]))
Gene_selected=GeneSymbol[ind_selected]

NodeID=Gene_selected
###  Reform to txt for Cytoscape

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
    }
  }
}

write.table(C, file="Results/Cytoscape_BreastCancer_GSE7390_Nerwork_new.txt",quote=F,row.names=F,col.names = F,sep = "\t")
