setwd("F:/Clinical Gene expression network Project/Reversion/Codes/Benchmarking_LPS_scRNAseq_data/")

scExpression_time=readMat("input_data.mat")
scExpression_time=matrix(unlist(scExpression_time),dim(as.data.frame(scExpression_time )))
Capture_time=scExpression_time[dim(scExpression_time)[1],]
PT=readMat("DPT.mat")
PT=as.numeric(unlist(PT))

### Load the package or install if not present
# if (!require("RColorBrewer")) {
#   install.packages("RColorBrewer")
  library(RColorBrewer)
# }

library(vioplot)
dev.new()
vioplot(PT[Capture_time==1],PT[Capture_time==2],PT[Capture_time==4],PT[Capture_time==6],names=c("1h","2h","4h","6h"),col=brewer.pal(4,"Set3"),outer = FALSE,horizontal=FALSE,
        xlab="Capture time", ylab="Pseudotemporal progression")
# pdf(file="F:/Clinical Gene expression network Project/Reversion/Codes/Benchmarking_LPS_scRNAseq_data/Results/Pseudotemporal progression_vs_Capture time.pdf", width=6, height=6) #矢量图，pdf格式


mydata <- data.frame(y=PT, x1=Capture_time)
my_lm <- lm(y ~ x1, data = mydata)
my_lm
my_slm <- summary(my_lm)
names(my_slm)

rsquared <- my_slm$r.squared
rsquared
temp <- locator(1) # 在图表上，你喜欢的地方点击一下，文字就出来了
text(temp, paste(expression(R^2),"=",round(rsquared,3)))

## Heatmap
setwd("F:/Clinical Gene expression network Project/Reversion/Codes/Benchmarking_LPS_scRNAseq_data/")
X_stage <- readMat("scExp_Net_Ctime.mat")
X_stage=as.data.frame(X_stage)
Expression=X_stage[1:(dim(X_stage)[1]-1),]
rownames(Expression)=TFs


Data=as.matrix(Expression)
colnames(Data)=NULL
rownames(Data)=TFs


library(pheatmap)
## Normalized
M=apply(Data,1,max)
A=apply(Data,1,mean)
S=apply(Data,1,sd)
D=(Data-A)/S
dev.new()
pheatmap(D,cluster_row=T, cluster_cols=F, clustering_distance_rows='euclidean',clustering_method = "ward", color = colorRampPalette(c("CornflowerBlue", "white", "firebrick3"))(200), fontsize=12, fontsize_row=8,labRow=NA, show_colnames = T,cellwidth = 0.5, cellheight = 20) #自定义颜色



### Cytoscape visualization
setwd("F:/Clinical Gene expression network Project/Reversion/Codes/Benchmarking_LPS_scRNAseq_data/")

AM=read.csv('Results/AdjacentMatrix_LPS.csv',header=F)
AM=as.matrix(AM)
# AM=AM[abs(AM)>max(as.numeric(abs(AM))*0.5),]

AM_new=AM*0
for (i in 1:nrow(AM))
{
  for (j in 1:ncol(AM))
  {
    if (abs(AM[i,j])>max(as.numeric(abs(AM))*0.05))
    {
      AM_new[i,j] = as.numeric(AM[i,j])
    }
  }
}
AM=AM_new

NodeID=TFs
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

write.table(C, file="Results/Cytoscape_DC_LPS_Nerwork_new.txt",quote=F,row.names=F,col.names = F,sep = "\t")

