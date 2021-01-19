
Data=read.csv("F:/Clinical Gene expression network Project/Code/TCGA Data/Breast cancer/GeneExpression_Selected_Smoothed.csv",header = F)  # Differential network genes
Data=read.csv("F:/Clinical Gene expression network Project/Code/TCGA Data/Breast cancer/GeneExpression_Selected_NonSmoothed.csv",header = F)  # Differential network genes


GeneName=read.csv("F:/Clinical Gene expression network Project/Code/TCGA Data/Breast cancer/SelectedGeneName_GEO_BreastCancer.txt",header = F)
rownames(Data)=as.character(unlist(GeneName))

require(stats)
library(pheatmap)
dev.off()

M=apply(Data,1,max)
A=apply(Data,1,mean)
S=apply(Data,1,sd)
D=(Data-A)/S
dev.new()
pheatmap(D,cluster_row=T, cluster_cols=F, clustering_distance_rows='euclidean',clustering_method = "ward", color = colorRampPalette(c("CornflowerBlue", "white", "firebrick3"))(200), fontsize=9, fontsize_row=6,labRow=NA, show_colnames = FALSE) #自定义颜色


