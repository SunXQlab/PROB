
Cytoscape_Reformat <- function(AM,NodeID)
  # AM - adjacent matrix
  # NodeID - the ID or the symbol  of the genes in the AM
  
  {
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
  return(C)
}

# Cytoscape_Reformat(AM,NodeID)