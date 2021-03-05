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



## rewired network
Net_gain=as.data.frame(cbind(as.character(C_gain[,1]),as.character(C_gain[,3]),C_gain[,2]))  # Input data for R analysis, data.frame
net_gain=graph_from_data_frame(Net_gain, directed = T, vertices = NULL) 
Degree_gain=degree(net_gain,mode="out")
sort(Degree_gain)

Net_lost=as.data.frame(cbind(as.character(C_lost[,1]),as.character(C_lost[,3]),C_lost[,2]))  # Input data for R analysis, data.frame
net_lost=graph_from_data_frame(Net_lost, directed = T, vertices = NULL) 
Degree_lost=degree(net_lost,mode="out")
sort(Degree_lost)

Degree=cbind(Degree_lost[unique(EMT_reggene_names)],Degree_gain[unique(EMT_reggene_names)])

setwd("F:/Clinical Gene expression network Project/Reversion/Codes/UC_SARC")  
write.csv(Degree, "Results/Table of out degree.csv")

## Heatmap
GeneNames=unlist(EMT_reggene_names)  # code in Cytoscape visualization.R
# library(R.matlab)
PT1=readMat("DPT_sort1.mat")
PT1=as.numeric(unlist(PT1))

PT2=readMat("DPT_sort2.mat")
PT2=as.numeric(unlist(PT2))

install.packages("R.matlab")
library(R.matlab)
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


## Pseudotemporal Dynamics of the selected genes
PT_smoothed=readMat("DPP.mat")
PT_smoothed=as.numeric(unlist(PT_smoothed))

NetGene_smoothed=readMat("NetGene_smoothed.mat")
Data_smoothed=as.matrix(as.data.frame(NetGene_smoothed))
rownames(Data_smoothed)=unique(GeneNames)

ACSS1=Data_smoothed['ACSS1',]
PTPN12=Data_smoothed['PTPN12',]
CDH1=Data_smoothed['CDH1',]
PERP=Data_smoothed['PERP',]

ACSS1=ACSS1/ACSS1[1]
PTPN12=PTPN12/PTPN12[1]
CDH1=CDH1/CDH1[1]
PERP=PERP/PERP[1]

library(RColorBrewer)

dev.new()
plot(PT_smoothed,ACSS1, type = "l", col = "#377EB8", lwd = 3, xlim=c(0,1), ylim=c(0.85,1.3),xlab="Pseudotemporal progression",ylab="Expression relative to t=0")
par(new=T)
lines(rep(PT_smoothed[79],length(seq(0.65,1.25,0.01))),seq(0.65,1.25,0.01),lty=1 , col = "#808A87", lwd = 2)
par(new=T)
lines(PT_smoothed,PERP, type = "l", col = "#984EA3", lwd = 3)
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

## EMT score 
output_data=readMat("output_data.mat")
GeneExpression_ordered=as.matrix(as.data.frame(output_data))


EMT_siggene_expression=GeneExpression_ordered[unique(EMT_siggene_ind),]
rownames(EMT_siggene_expression)=unique(EMT_siggene_names)

# CDH1
CDH1=EMT_siggene_expression['CDH1',]
weight=cor(t(EMT_siggene_expression),CDH1)
EMT_score=t(EMT_siggene_expression)%*%weight

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

### Predicting phenotypes (epithelial vs mesenchymal)
# > ind_training
# [1]  69  52  15  34  27   6  73  50  44  82  24  30  29  63  54  75  42  17  32  40  65  56  46  22  58  51  68   7  70
# [30]  74  67  21  80  41  13  53  11  25  47  62  18   9 100  94 103  97 112  98 102 108  86  92  89 111 104 109
# > ind_test
# [1]   1   2   3   4   5   8  10  12  14  16  19  20  23  26  28  31  33  35  36  37  38  39  43  45  48  49  55  57  59
# [30]  60  61  64  66  71  72  76  77  78  79  81  83  84  85  87  88  90  91  93  95  96  99 101 105 106 107 110

ACSS1=Data['ACSS1',]
PTPN12=Data['PTPN12',]
PERP=Data['PERP',]
EMT_score=as.numeric(EMT_score<0)  # 0 for epithelial phenotype; 1 for mesenchymal phenotype
# ind_training=c(sample(c(1:84),42),sample(c(85:112),14))
# ind_test=setdiff(c(1:112),ind_training)
EMT_score_Training=EMT_score[ind_training]
EMT_score_Test=EMT_score[ind_test]
ACSS1_Training=ACSS1[ind_training]
PTPN12_Training=PTPN12[ind_training]
ACSS1_Test=ACSS1[ind_test]
PTPN12_Test=PTPN12[ind_test]
##
PERP_Training=PERP[ind_training]
PERP_Test=PERP[ind_test]
# library("stats")
EMT_model=glm(EMT_score_Training ~ PTPN12_Training + ACSS1_Training, family = binomial(link="logit"))
coefficients(EMT_model)
summary(EMT_model)
TestData=data.frame(PTPN12_Test,ACSS1_Test)
Prediction = predict(EMT_model,newx = TestData,type="response")

# library(ROCR)
pred <- prediction(Prediction,EMT_score_Test)

perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
performance(pred,"auc") # shows calculated AUC for model
dev.new()
plot(perf,colorize=FALSE, col="DarkCyan") # plot ROC curve
lines(c(0,1),c(0,1),col = "gray", lty = 4 )

# library(pROC)
roc0=roc(EMT_score_Test,as.numeric(Prediction))
AUC=as.numeric(auc(EMT_score_Test,as.numeric(Prediction)))
AUC

cor.test(ACSS1,EMT_score, method = "spearman")  # c("pearson", "kendall", "spearman")
cor.test(PERP,EMT_score, method = "spearman")

ACSS1_p=ACSS1[EMT_score==0]
ACSS1_m=ACSS1[EMT_score==1]
PTPN12_p=PTPN12[EMT_score==0]
PTPN12_m=PTPN12[EMT_score==1]
CDH1_p=CDH1[EMT_score==0]
CDH1_m=CDH1[EMT_score==1]
dev.new()
plot(CDH1_p,PTPN12_p,col="red") #B3DE69
par(new=T)
points(CDH1_m,PTPN12_m,col="blue")  #80B1D3

### Predicting histological subtypes (UC vs. SARC)
ACSS1_UC=ACSS1[1:84]
ACSS1_SARC=ACSS1[85:112]
PTPN12_UC=PTPN12[1:84]
PTPN12_SARC=PTPN12[85:112]
PERP_UC=PTPN12[1:84]
PERP_SARC=PTPN12[85:112]

dev.new()
plot(ACSS1_UC,PTPN12_UC,col="red") #B3DE69
par(new=T)
points(ACSS1_SARC,PTPN12_SARC,col="blue")  #80B1D3

# ind_training=c(sample(c(1:84),42),sample(c(85:112),14))
# ind_test=setdiff(c(1:112),ind_training)

Phenotype=c(rep(0,84),rep(1,28))
Phenotype_Training=Phenotype[ind_training]
Phenotype_Test=Phenotype[ind_test]
ACSS1_Training=ACSS1[ind_training]
PTPN12_Training=PTPN12[ind_training]
ACSS1_Test=ACSS1[ind_test]
PTPN12_Test=PTPN12[ind_test]

SARC_model=glm(Phenotype_Training ~ PTPN12_Training + ACSS1_Training , family = binomial(link="logit"))
coefficients(SARC_model)
summary(SARC_model)
TestData=data.frame(PTPN12_Test,ACSS1_Test)
Prediction = predict(SARC_model,newx = TestData,type="response")

# library(ROCR)
pred <- prediction(Prediction,Phenotype_Test)

perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
performance(pred,"auc") # shows calculated AUC for model
dev.new()
plot(perf,colorize=FALSE, col="DarkCyan") # plot ROC curve
lines(c(0,1),c(0,1),col = "gray", lty = 4 )

# library(pROC)
roc0=roc(EMT_score_Test,as.numeric(Prediction))
AUC=as.numeric(auc(Phenotype_Test,as.numeric(Prediction)))
AUC

### Differential expression analysis of ACSS1 in metastatic UCs
TCGAL_BLCA = read.csv("F:/Clinical Gene expression network Project/Reversion/Data/UC-SARC/TCGA_BLCA_logFPKM.csv")
#Metastasis
BLCA_phenotype = read.csv("F:/Clinical Gene expression network Project/Reversion/Data/UC-SARC/TCGA_BLCA_phenotype.csv")
Patient_ID=BLCA_phenotype[,1]
BLCA_event=BLCA_phenotype[,"new_neoplasm_event_type"]
BLCA_event=as.numeric(BLCA_event=='Metastatic')
names(BLCA_event)=Patient_ID

ACSS1=TCGAL_BLCA[TCGAL_BLCA[,1]=="ACSS1",]
ACSS1_metasis=ACSS1[,intersect(Patient_ID[BLCA_event==1],colnames(ACSS1))]
ACSS1_NOmetasis=ACSS1[,intersect(Patient_ID[BLCA_event==0],colnames(ACSS1))]
ACSS1_metasis=2^as.numeric(ACSS1_metasis)
ACSS1_NOmetasis=2^as.numeric(ACSS1_NOmetasis)

#Grade
BLCA_Grade=BLCA_phenotype[,"neoplasm_histologic_grade"]
names(BLCA_Grade)=Patient_ID
BLCA_Grade=BLCA_Grade[BLCA_Grade=="High Grade"|BLCA_Grade=="Low Grade"]
Patient_ID=names(BLCA_Grade)
BLCA_Grade=as.numeric(BLCA_Grade=='High Grade')
names(BLCA_Grade)=Patient_ID

# ACSS1
ACSS1=TCGAL_BLCA[TCGAL_BLCA[,1]=="ACSS1",]
ACSS1_HighGrade=ACSS1[,intersect(Patient_ID[BLCA_Grade==1],colnames(ACSS1))]
ACSS1_LowGrade=ACSS1[,intersect(Patient_ID[BLCA_Grade==0],colnames(ACSS1))]

ACSS1_HighGrade=2^as.numeric(ACSS1_HighGrade)
ACSS1_LowGrade=2^as.numeric(ACSS1_LowGrade)

WT=wilcox.test(ACSS1_HighGrade,ACSS1_LowGrade,alternative="less", paired=F)
dev.new()
boxplot(ACSS1_LowGrade,ACSS1_HighGrade,names=c("Low grade UC","High grade UC"),col=c("#B3DE69","#80B1D3"),outline=F,boxwex = 0.4,boxwex=0.5,width=c(0.5,0.5),
        xlab="", ylab="Gene expression")

#PTPN12
PTPN12=TCGAL_BLCA[TCGAL_BLCA[,1]=="PTPN12",]
PTPN12_HighGrade=PTPN12[,intersect(Patient_ID[BLCA_Grade==1],colnames(PTPN12))]
PTPN12_LowGrade=PTPN12[,intersect(Patient_ID[BLCA_Grade==0],colnames(PTPN12))]

PTPN12_HighGrade=2^as.numeric(PTPN12_HighGrade)
PTPN12_LowGrade=2^as.numeric(PTPN12_LowGrade)

WT=wilcox.test(PTPN12_HighGrade,PTPN12_LowGrade,alternative="greater", paired=F)
dev.new()
boxplot(PTPN12_LowGrade,PTPN12_HighGrade,names=c("Low grade UC","High grade UC"),col=c("#B3DE69","#80B1D3"),outline=F,boxwex = 0.4, width=c(0.5,0.5),
        xlab="", ylab="Gene expression")

CDH1=TCGAL_BLCA[TCGAL_BLCA[,1]=="CDH1",]
CDH1_HighGrade=CDH1[,intersect(Patient_ID[BLCA_Grade==1],colnames(CDH1))]
CDH1_LowGrade=CDH1[,intersect(Patient_ID[BLCA_Grade==0],colnames(CDH1))]

CDH1_HighGrade=2^as.numeric(CDH1_HighGrade)
CDH1_LowGrade=2^as.numeric(CDH1_LowGrade)

WT=wilcox.test(CDH1_HighGrade,CDH1_LowGrade,alternative="less", paired=F)
dev.new()
boxplot(CDH1_LowGrade,CDH1_HighGrade,names=c("Low grade UC","High grade UC"),col=c("#B3DE69","#80B1D3"),outline=F,boxwex = 0.4, width=c(0.5,0.5),
        xlab="", ylab="Gene expression")
