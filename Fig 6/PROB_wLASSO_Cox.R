
library(glmnet)  # manually load this package by downloading and installing from Toll. 
###GEO data
setwd("F:/Clinical Gene expression network Project/Manuscript/Codes/More cancers/Breast cancer GEO/")  


## Breast cancer (GSE7390) 
# Gene1=read.csv("Data/GSE7390_series_matrix.csv")
# Gene1_Expression=Gene1[-c(1:113),]
# Gene1_Expression=as.matrix(Gene1_Expression)
# rownames(Gene1_Expression)=Gene1_Expression[,1]
# colnames(Gene1_Expression)=Gene1_Expression[1,]
# Gene1_Expression=Gene1_Expression[-c(1,22284),-1]
# 
# 
# Clinic1=read.csv("Data/GSE7390_clinical information.csv")
# Grade1=Clinic1[,c(13)]
# names(Grade1)=Clinic1[,c(3)]
# 
# Grade1=na.omit(Grade1)
# 
# Gene1_Expression=Gene1_Expression[,intersect(colnames(Gene1_Expression),names(Grade1))]
# Grade1=Grade1[intersect(colnames(Gene1_Expression),names(Grade1))]
# 
# ProbeID=read.csv("Data/GPL96_ProbeIDannotation.csv")
# ProbeID=ProbeID[-c(1:16),c(1,11)]
# ProbeID=as.matrix(ProbeID)
# rownames(ProbeID)=ProbeID[,2]
# 
# Gene_Name=rownames(ProbeID[match(ProbeID[,1],rownames(Gene1_Expression)),])
# 
# rownames(Gene1_Expression)=Gene_Name

## Survival analysis
# Time=Clinic1[,c(15)]
# Status=Clinic1[,c(16)]
# names(Time)=Clinic1[,c(3)]
# names(Status)=Clinic1[,c(3)]
# 
# 
# Gene1_Expression=Gene1_Expression[,intersect(colnames(Gene1_Expression),names(Time))]
# Time=Time[intersect(colnames(Gene1_Expression),names(Time))]
# Status=Status[intersect(colnames(Gene1_Expression),names(Status))]

NodeScore=read.csv("Results/Node_Importance.csv",header = FALSE)
# Gene_marker=matrix(as.numeric(Gene1_Expression[NodeScore[,2],]),dim(Gene1_Expression[NodeScore[,2],]))
# Gene_marker=as.numeric(Gene1_Expression['FOXM1',])

ind_selected=read.csv("Results/ind_selected.csv",header = FALSE)
ind_selected=as.numeric((ind_selected))
Gene_marker=matrix(as.numeric(Gene1_Expression[ind_selected,]),dim(Gene1_Expression[ind_selected,]))

S0=as.numeric(NodeScore[,1])
S=S0
S=max(S)-S
# S=(S-min(S))/(max(S)-min(S))
# S=1/S
Gene_marker_w=Gene_marker  # Importance score-weighted expression
surv=Surv(as.double(Time),as.double(Status))
cvfit_wLASSO = cv.glmnet(t(Gene_marker_w),surv,family="cox",alpha=1,nfolds=3)  #,penalty.factor=S

Prediction = predict(cvfit_wLASSO,newx = t(Gene_marker_w),s="lambda.min")
pred <- prediction(Prediction,Status)
# perf_Training_wLASSO[i] <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
roc0=roc(Status,as.numeric(Prediction),quiet = TRUE)
as.numeric(auc(Status,as.numeric(Prediction)))

dev.new()
plot(cvfit_wLASSO)
A=coef(cvfit_wLASSO, s = "lambda.min")
A=coef(cvfit, s="lambda.1se")
A=as.numeric(A)
A[(A!=0)]
which(A!=0)
sum(A!=0)
Biomarker=rownames(Gene1_Expression)[which(A!=0)]


