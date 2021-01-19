###GEO data
setwd("F:/Clinical Gene expression network Project/Manuscript/Codes/Fig 4-5_Realistic applications/Breast cancer GEO_Fig 5/")  


## Breast cancer 
Gene=read.csv("Data/Validation/GSE40766_series_matrix.csv")
Gene_Expression=Gene[-c(1:65),-1]
Sample=Gene[30,-1]
Sample=(as.matrix(Sample))

Sample=as.character(Sample)
# Sample=gsub("mutation: ", "", Sample, fixed = TRUE)

rownames(Gene_Expression)=Gene[-c(1:65),1]
colnames(Gene_Expression)=Sample
Gene_Expression=Gene_Expression[-1,]


### ProbeID <--> GeneSymbol
probeid_geneid<-read.csv("Data/Validation/GPL10558-50081.csv")
probeid_geneid=probeid_geneid[-(1:31),c(1,6)]
probeid_geneid=as.matrix(probeid_geneid)
rownames(probeid_geneid)=probeid_geneid[,2]

### Selected gene expression
Gene_Selected_Name=c('FOXM1','ASPM',"CDCA8",'KIF2C','MCM10','MELK',"NCAPG",'SHCBP1',"STIL")
ProbeID=probeid_geneid[Gene_Selected_Name,]
Gene_Expression=Gene_Expression[ProbeID[,1],]

rownames(Gene_Expression)=ProbeID[,2]


# write.csv(Gene_Expression,file="Data/Validation/Gene_Selected_GSE2222.csv")

### Significant differential expression (BRCA1 mutation vs normal)
# Gene_Selected_siRNA_control=read.csv("Breast Cancer/Gene_Selected_siRNA_control.csv")
pvalue=matrix(0,length(Gene_Selected_Name),1)
Inc_Dec=matrix(0,length(Gene_Selected_Name),1)

for (i in 1:length(Gene_Selected_Name))
{
  G=as.numeric(as.matrix(Gene_Expression[i,seq(1,11,2)]))
  C=as.numeric(as.matrix(Gene_Expression[i,seq(2,12,2)]))
  Inc_Dec[i]=mean(C)-mean(G)
  if (Inc_Dec[i]>0)
  {
    T=wilcox.test(C,G,alternative = 'greater',paired = 0)  #wilcox.test(C-G)
  }
  
  if (Inc_Dec[i]<0)
  {
    T=wilcox.test(C,G,alternative = 'less',paired = 0)  #wilcox.test(C-G)
  }
  
  pvalue[i]=T$p.value
}

which( pvalue<0.05)
Gene_Selected_Name[pvalue<0.05]

## T test
pvalue=matrix(0,length(Gene_Selected_Name),1)
Inc_Dec=matrix(0,length(Gene_Selected_Name),1)
for (i in 1:length(Gene_Selected_Name))
{
  G=as.numeric(as.matrix(Gene_Expression[i,seq(1,11,2)]))
  C=as.numeric(as.matrix(Gene_Expression[i,seq(2,12,2)]))
  Inc_Dec[i]=mean(C)-mean(G)
  if (Inc_Dec[i]>0)
  {
    T=t.test(C,G,alternative = 'greater')  #wilcox.test(C-G)
  }
  
  if (Inc_Dec[i]<0)
  {
    T=t.test(C,G,alternative = 'less')  #wilcox.test(C-G)
  }
  
  pvalue[i]=T$p.value
}
which( pvalue<0.05)
Gene_Selected_Name[pvalue<0.05]

# write.csv(pvalue,file="F:/Clinical Gene expression network Project/Manuscript/Codes/Fig 4_5/Breast Cancer/pvalue_GeneExpression_siRNA_control.csv")
# write.csv(Inc_Dec,file="F:/Clinical Gene expression network Project/Manuscript/Codes/Fig 4_5/Breast Cancer/Increase_Decrease_GeneExpression_siRNA_control.csv")


### Visualization
Gene_Control=matrix(as.numeric(as.matrix(Gene_Expression[1:9,seq(2,12,2)])),9,6)
Gene_Case=matrix(as.numeric(as.matrix(Gene_Expression[1:9,seq(1,11,2)])),9,6)

Data=as.data.frame(cbind(Gene_Control,Gene_Case))
Data=t(na.omit(Data))

type=as.factor(c(rep('DMSO',6),rep('Thiostrepton',6)))
x=cbind(Data,(type))

# install.packages("vioplot")
library(vioplot)

D=x[x[,10]==1,1:9]
E=x[x[,10]==2,1:9]
## Boxplot
## brown and blue
dev.new()
boxplot(D[,1],D[,2],D[,3],D[,4],D[,5],D[,6],D[,7],D[,8],D[,9],at=seq(1,26,3), names=paste(as.character(Gene_Selected_Name), 'C',sep='_'), xlim=c(0,27), ylim=c(7,12),col="DarkKhaki", border="gray51", lty=1, lwd=1, rectCol="gray", 
        colMed="white", pchMed=19)
boxplot(E[,1],E[,2],E[,3],E[,4],E[,5],E[,6],E[,7],E[,8],E[,9],add=TRUE,at=seq(2,26,3), names=paste(as.character(Gene_Selected_Name), 'T',sep='_'), xlim=c(0,27), ylim=c(7,12),col="DeepSkyBlue3", border="gray51", lty=1, lwd=1, rectCol="gray", 
        colMed="white", pchMed=19)

## blue and red 
dev.new()
boxplot(D[,2],D[,3],D[,4],D[,5],D[,6],D[,7],D[,8],D[,9],at=seq(1,29,4), names=rep('DMSO',8), xlim=c(0,31), ylim=c(7,12),col="#00AFBB", border="gray51", lty=1, lwd=0.5, rectCol="gray", 
        colMed="white", pchMed=19)
boxplot(E[,2],E[,3],E[,4],E[,5],E[,6],E[,7],E[,8],E[,9],add=TRUE,at=seq(2,30,4), names=rep('Thiostrepton',8), xlim=c(0,31), ylim=c(7,12),col="#FC4E07", border="gray51", lty=1, lwd=0.5, rectCol="gray", 
        colMed="white", pchMed=19)

## violin
dev.new()
vioplot(D[,1],D[,2],D[,3],D[,4],D[,5],D[,6],D[,7],D[,8],D[,9],at=seq(1,26,3), names=paste(as.character(Gene_Selected_Name), 'C',sep='_'), xlim=c(0,26),col="DarkKhaki", border="gray51", lty=1, lwd=1, rectCol="gray", 
        colMed="white")
vioplot(E[,1],E[,2],E[,3],E[,4],E[,5],E[,6],E[,7],E[,8],E[,9],add=TRUE,at=seq(2,26,3), names=paste(as.character(Gene_Selected_Name), 'T',sep='_'), xlim=c(0,26),col="DeepSkyBlue3", border="gray51", lty=1, lwd=1, rectCol="gray", 
        colMed="white")


