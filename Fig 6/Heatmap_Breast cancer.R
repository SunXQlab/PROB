
setwd("F:/Clinical Gene expression network Project/Manuscript/Codes/More cancers/Breast cancer GEO/")  

ind_selected=read.csv("Results/ind_selected.csv",header = FALSE)
ind_selected=as.numeric((ind_selected))

Gene_Expression=read.csv("Data/Gene_GSE7390.csv")
GeneSymbol=as.character(unlist(Gene_Expression[,1]))
Gene_selected=GeneSymbol[ind_selected]

Data_smooth_selected=read.csv("Results/Data_smooth_selected.csv", header = FALSE)
Data=as.matrix(Data_smooth_selected)
rownames(Data)=Gene_selected



# require(stats)
# library(pheatmap)
dev.off()

M=apply(Data,1,max)
A=apply(Data,1,mean)
S=apply(Data,1,sd)
D=(Data-A)/S
dev.new()
pheatmap(D,cluster_row=T, cluster_cols=F, clustering_distance_rows='euclidean',clustering_method = "ward", color = colorRampPalette(c("CornflowerBlue", "white", "firebrick3"))(200), fontsize=9, fontsize_row=6,labRow=NA, show_colnames = FALSE) #自定义颜色

kc = kmeans(dist(D) , 2)

GeneSet1 = Gene_selected[kc$cluster==1]
GeneSet2 = Gene_selected[kc$cluster==2]

G=list(1,2)
G[[1]]=as.character(GeneSet1)
G[[2]]=as.character(GeneSet2)
write.csv(rbind(GeneSet1,GeneSet2), file="Results/Clustered_GeneSet.csv")


#### Expression of FOXM1
Gene_Expression=as.matrix(Gene_Expression)
rownames(Gene_Expression)=GeneSymbol
Gene_Expression=Gene_Expression[,-1]
Grade=read.csv("Data/Grade_GSE7390.csv",header = TRUE)

FOXM1_G1=Gene_Expression['FOXM1',Grade[,2]==1]
FOXM1_G2=Gene_Expression['FOXM1',Grade[,2]==2]
FOXM1_G3=Gene_Expression['FOXM1',Grade[,2]==3]
# FOXM1_G4=Gene_Expression['FOXM1',Grade[,2]==4]

mean(as.numeric(FOXM1_G1))
mean(as.numeric(FOXM1_G2))
mean(as.numeric(FOXM1_G3))

# library(vioplot)
vioplot(as.numeric(FOXM1_G1),as.numeric(FOXM1_G2),as.numeric(FOXM1_G2), names=c("Stage 1","Stage 2","Stage 3"),col="gold",ylab="Expression level")
title("Stage-wise FOXM1 expression")

wilcox.test(as.numeric(FOXM1_G1),as.numeric(FOXM1_G2),alternative="less", paired=F)  # p-value = 0.04157

wilcox.test(as.numeric(FOXM1_G2),as.numeric(FOXM1_G3),alternative="less", paired=F)  # p-value = 1.644e-14

wilcox.test(as.numeric(FOXM1_G1),as.numeric(FOXM1_G3),alternative="less", paired=F)  # p-value = 1.644e-14

dataset <- data.frame(value = c(as.numeric(FOXM1_G1),as.numeric(FOXM1_G2),as.numeric(FOXM1_G3)), group = factor(rep(c("Stage 1", "Stage 2", "Stage 3"), times = c(length(FOXM1_G1), length(FOXM1_G2),length(FOXM1_G3)))))

dev.new()
boxplot( value ~ group,  notch = F, dataset, border = c( "#009E73","purple","red"),cex = 1,cex.axis=1,pars = list(boxwex = 0.5, staplewex = 0.5, outwex = 0.5),ylab="Expression level",ylim=c(3.8,13))  #,col.axis = "#009E73"

