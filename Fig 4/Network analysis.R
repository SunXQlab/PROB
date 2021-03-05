### Network analysis
####### Calculating hUb score
library(igraph)
## UC hub score
Net_UC=as.data.frame(cbind(as.character(C_UC[,1]),as.character(C_UC[,3]),C_UC[,2]))  # Input data for R analysis, data.frame
net=graph_from_data_frame(Net_UC, directed = T, vertices = NULL) 
AM=as_adjacency_matrix(net)

AM=AM_UC[V(net),V(net)]
AM=abs(AM)
dim(AM)
EG1=eigen(AM*t(AM))
which(EG1$values==max(EG1$values))
hs_UC=EG1$vectors[,1]
names(hs_UC)= rownames(as.matrix(V(net)))
hs_UC=hs_UC[unique(EMT_reggene_names)]
sort(hs_UC,decreasing = T)


## SARC hub score
Net_SARC=as.data.frame(cbind(as.character(C_SARC[,1]),as.character(C_SARC[,3]),C_SARC[,2]))  # Input data for R analysis, data.frame
net=graph_from_data_frame(Net_SARC, directed = T, vertices = NULL) 
AM=as_adjacency_matrix(net)

AM=AM_SARC[V(net),V(net)]
AM=abs(AM)
dim(AM)
EG2=eigen(AM*t(AM))
which(EG2$values==max(EG2$values))
hs_SARC=EG2$vectors[,1]
names(hs_SARC)= rownames(as.matrix(V(net)))
hs_SARC=hs_SARC[unique(EMT_reggene_names)]
sort(hs_SARC,decreasing = T)

## Change of hub score
# hs_SARC=(hs_SARC-min(hs_SARC))/(max(hs_SARC)-min(hs_SARC))
# hs_UC=(hs_UC-min(hs_UC))/(max(hs_UC)-min(hs_UC))
hs_delta=hs_SARC[unique(EMT_reggene_names)]-hs_UC[unique(EMT_reggene_names)]
# hs_delta=rank(hs_SARC)-rank(hs_UC)
sort(hs_delta,decreasing = T)

hs=cbind(hs_UC,hs_SARC,hs_delta,abs(hs_delta))
setwd("F:/Clinical Gene expression network Project/Reversion/Codes/UC_SARC")  
write.csv(hs, "Results/Table of hub score.csv")


## Differential network
Net_Diff=as.data.frame(cbind(as.character(C_Diff[,1]),as.character(C_Diff[,3]),C_Diff[,2]))  # Input data for R analysis, data.frame
net=graph_from_data_frame(Net_Diff, directed = T, vertices = NULL) 
AM=as_adjacency_matrix(net)

AM=AM_D[V(net),V(net)]
AM=abs(AM)
dim(AM)
EG3=eigen(AM*t(AM))
which(EG3$values==max(EG3$values))
hs=EG1$vectors[,1]
names(hs)= rownames(as.matrix(V(net)))
sort(hs,decreasing = T)

## rewired network
Net_rewired=as.data.frame(cbind(as.character(C_rewired[,1]),as.character(C_rewired[,3]),C_rewired[,2]))  # Input data for R analysis, data.frame
net_rewired=graph_from_data_frame(Net_rewired, directed = T, vertices = NULL) 
Degree_rewired=degree(net_rewired,mode="out")
sort(Degree_rewired)

Net_lost=as.data.frame(cbind(as.character(C_lost[,1]),as.character(C_lost[,3]),C_lost[,2]))  # Input data for R analysis, data.frame
net_lost=graph_from_data_frame(Net_lost, directed = T, vertices = NULL) 
Degree_lost=degree(net_lost,mode="out")
sort(Degree_lost)

Degree=cbind(Degree_lost[unique(EMT_reggene_names)],Degree_rewired[unique(EMT_reggene_names)])

setwd("F:/Clinical Gene expression network Project/Reversion/Codes/UC_SARC")  
write.csv(Degree, "Results/Table of out degree.csv")

## Heatmap
GeneNames=unlist(EMT_reggene_names)  # code in Cytoscape visualization.R
# library(R.matlab)
PT1=readMat("DPT_sort1.mat")
PT1=as.numeric(unlist(PT1))

PT2=readMat("DPT_sort2.mat")
PT2=as.numeric(unlist(PT2))

NetGene1=readMat("NetGene1.mat")
NetGene2=readMat("NetGene2.mat")

NetGene1=as.data.frame(NetGene1)
NetGene2=as.data.frame(NetGene2)


# install.packages("gplots") #下载gplots程序包
library(gtools) #加载gplots程序
library(gplots) #加载gplots程序

Data=as.matrix(cbind(NetGene1,NetGene2))
colnames(Data)=NULL
rownames(Data)=unique(GeneNames)


library(pheatmap)
## Normalized
M=apply(Data,1,max)
A=apply(Data,1,mean)
S=apply(Data,1,sd)
D=(Data-A)/S
dev.new()
pheatmap(D,cluster_row=T, cluster_cols=F, clustering_distance_rows='euclidean',clustering_method = "ward", color = colorRampPalette(c("CornflowerBlue", "white", "firebrick3"))(200), fontsize=12, fontsize_row=8,labRow=NA, show_colnames = T,cellwidth = 1.5, cellheight = 8) #自定义颜色

# pheatmap(D,cluster_row=T, cluster_cols=F,  clustering_distance_rows='euclidean', cellwidth = 1.5, cellheight = 8,color = colorRampPalette(c("navy", "white", "firebrick3"))(100), fontsize=9, fontsize_row=6,labRow=NA) #自定义颜色

## Degree distribution
Net_UC=as.data.frame(cbind(as.character(C_UC[,1]),as.character(C_UC[,3]),C_UC[,2]))  # Input data for R analysis, data.frame
net_UC=graph_from_data_frame(Net_UC, directed = T, vertices = NULL) 
Net_SARC=as.data.frame(cbind(as.character(C_SARC[,1]),as.character(C_SARC[,3]),C_SARC[,2]))  # Input data for R analysis, data.frame
net_SARC=graph_from_data_frame(Net_SARC, directed = T, vertices = NULL) 
Deg1=degree(net_UC)
Deg2=degree(net_SARC)
F1=hist(Deg1,freq=F,breaks = length(Deg1))
F2=hist(Deg2,freq=F,breaks = length(Deg2))

WT=wilcox.test(Deg1,Deg2,alternative="greater", paired=F)


dev.new()
boxplot(Deg1,Deg2,names=c("UC network","SARC network"),col=c("#B3DE69","#80B1D3"),outline=F,boxwex = 0.4,boxwex=0.5,width=c(0.5,0.5),
        xlab="", ylab="Degree")

temp <- locator(1) # 在图表上，你喜欢的地方点击一下，文字就出来了
text(temp, paste("p-value =",WT$p.value))



### UC_SARC
Net_UC_SARC=as.data.frame(cbind(as.character(C_UC_SARC[,1]),as.character(C_UC_SARC[,3]),C_UC_SARC[,2]))  # Input data for R analysis, data.frame
net=graph_from_data_frame(Net_UC_SARC, directed = T, vertices = NULL) 
AM=as_adjacency_matrix(net)

AM=AM_UC_SARC[V(net),V(net)]
AM=abs(AM)
dim(AM)
EG1=eigen(AM*t(AM))
which(EG1$values==max(EG1$values))
hs=EG1$vectors[,1]
names(hs)= rownames(as.matrix(V(net)))
sort(hs,decreasing = T)

## EMT score 
output_data=readMat("output_data.mat")
GeneExpression_ordered=as.matrix(as.data.frame(output_data))


EMT_siggene_expression=GeneExpression_ordered[unique(EMT_siggene_ind),]
rownames(EMT_siggene_expression)=unique(EMT_siggene_names)

# CDH1
CDH1=EMT_siggene_expression['CDH1',]
weight=cor(t(EMT_siggene_expression),CDH1)
EMT_score1=t(EMT_siggene_expression)%*%weight
# 
# # SNAI2
# SNAI2=EMT_siggene_expression['CDH1',]
# weight=cor(t(EMT_siggene_expression),CDH1)
# EMT_score2=t(EMT_siggene_expression)%*%weight
# 
# # SNAI3
# SNAI3=EMT_siggene_expression['CDH1',]
# weight=cor(t(EMT_siggene_expression),CDH1)
# EMT_score3=t(EMT_siggene_expression)%*%weight
# # TWIST1
# TWIST1=EMT_siggene_expression['CDH1',]
# weight=cor(t(EMT_siggene_expression),CDH1)
# EMT_score4=t(EMT_siggene_expression)%*%weight
# 
# # TJP1
# TJP1=EMT_siggene_expression['CDH1',]
# weight=cor(t(EMT_siggene_expression),CDH1)
# EMT_score5=t(EMT_siggene_expression)%*%weight
# # CLDN1
# CLDN1=EMT_siggene_expression['CDH1',]
# weight=cor(t(EMT_siggene_expression),CDH1)
# EMT_score6=t(EMT_siggene_expression)%*%weight
# 
# EMT_score=EMT_score1+EMT_score2+EMT_score3+EMT_score4+EMT_score5+EMT_score6
EMT_score=EMT_score1
EMT_score=EMT_score-mean(EMT_score)
PT=c(PT1,PT2)
PT=PT/max(PT)

dev.new()
plot(EMT_score)
scatter.smooth(EMT_score)

EMT_score_UC=EMT_score[1:length(PT1)]
EMT_score_SARC=EMT_score[(length(PT1)+1):length(PT)]

WT=wilcox.test(EMT_score_UC,EMT_score_SARC,alternative="greater", paired=F)
dev.new()
boxplot(EMT_score_UC,EMT_score_SARC,names=c("UC","SARC"),col=c("#B3DE69","#80B1D3"),outline=F,boxwex = 0.4,boxwex=0.5,width=c(0.5,0.5),
        xlab="", ylab="EMT score")

temp <- locator(1) # 在图表上，你喜欢的地方点击一下，文字就出来了
text(temp, paste("p-value =",WT$p.value))

## Pseudotemporal Dynamics of the selected genes
PT_smoothed=readMat("DPP.mat")
PT_smoothed=as.numeric(unlist(PT_smoothed))

NetGene_smoothed=readMat("NetGene_smoothed.mat")
Data_smoothed=as.matrix(as.data.frame(NetGene_smoothed))
rownames(Data_smoothed)=unique(GeneNames)

# GNAI1=Data_smoothed['GNAI1',]
# WNT5B=Data_smoothed['WNT5B',]
# TP63=Data_smoothed['TP63',]
# CDH1=Data_smoothed['CDH1',]
# 
# GNAI1=GNAI1/GNAI1[1]
# WNT5B=WNT5B/WNT5B[1]
# TP63=TP63/TP63[1]
# CDH1=CDH1/CDH1[1]
# 
# library(RColorBrewer)
# 
# dev.new()
# plot(PT_smoothed,GNAI1, type = "l", col = "#4DAF4A", lwd = 5, xlim=c(0,1), ylim=c(0.7,1.3),xlab="Pseudotemporal progression",ylab="Expression relative to t=0")
# par(new=T)
# lines(rep(PT_smoothed[79],length(seq(0.65,1.25,0.01))),seq(0.65,1.25,0.01),lty=2 , col = "#808A87", lwd = 3)
# par(new=T)
# lines(PT_smoothed,WNT5B, type = "l", col = "#377EB8", lwd = 5)
# par(new=T)
# lines(PT_smoothed,TP63, type = "l", col = "#984EA3", lwd = 5)
# par(new=T)
# lines(PT_smoothed,CDH1, type = "l", col = "#E41A1C", lwd = 5)
# legend("topleft", legend = c('GNAI1','WNT5B','TP63','CDH1'),
#        col =  c("#4DAF4A","#377EB8","#984EA3","#E41A1C"), 
#        lty = 1,  lwd=6, bty="n",
#        inset = 0.01, xpd = TRUE, horiz = TRUE)
# ## test for expression changes between UC and SARC
# Pvalue=vector("double",4)
# WT=wilcox.test(GNAI1[1:79],GNAI1[80:102],alternative="greater", paired=F)
# Pvalue[1]=WT$p.value
# WT=wilcox.test(WNT5B[1:79],WNT5B[80:102],alternative="less", paired=F)
# Pvalue[2]=WT$p.value
# WT=wilcox.test(TP63[1:79],TP63[80:102],alternative="greater", paired=F)
# Pvalue[3]=WT$p.value
# WT=wilcox.test(CDH1[1:79],CDH1[80:102],alternative="greater", paired=F)
# Pvalue[4]=WT$p.value

### network rewiring

ACSS1=Data_smoothed['ACSS1',]
PTPN12=Data_smoothed['PTPN12',]
CDH1=Data_smoothed['CDH1',]

ACSS1=ACSS1/ACSS1[1]
PTPN12=PTPN12/PTPN12[1]
CDH1=CDH1/CDH1[1]

library(RColorBrewer)

dev.new()
plot(PT_smoothed,ACSS1, type = "l", col = "#377EB8", lwd = 3, xlim=c(0,1), ylim=c(0.85,1.3),xlab="Pseudotemporal progression",ylab="Expression relative to t=0")
par(new=T)
lines(rep(PT_smoothed[79],length(seq(0.65,1.25,0.01))),seq(0.65,1.25,0.01),lty=1 , col = "#808A87", lwd = 2)
par(new=T)
lines(PT_smoothed,PTPN12, type = "l", col = "#984EA3", lwd = 3)
par(new=T)
lines(PT_smoothed,CDH1, type = "l", col = "#E41A1C", lwd = 3)
legend("topleft", legend = c('ACSS1','PTPN12','CDH1'),
       col =  c("#377EB8","#984EA3","#E41A1C"), 
       lty = 1,  lwd=6, bty="n",
       inset = 0.01, xpd = TRUE, horiz = F)
## test for expression changes between UC and SARC
Pvalue=vector("double",3)
WT=wilcox.test(ACSS1[1:79],ACSS1[80:102],alternative="less", paired=F)
Pvalue[1]=WT$p.value
WT=wilcox.test(PTPN12[1:79],PTPN12[80:102],alternative="greater", paired=F)
Pvalue[2]=WT$p.value
WT=wilcox.test(CDH1[1:79],CDH1[80:102],alternative="greater", paired=F)
Pvalue[3]=WT$p.value



