source("http://bioc.ism.ac.jp/biocLite.R")
biocLite("AnnotationDbi")
install.packages("AnnotationDbi")
install.packages("bit")
library(bit)
library(AnnotationDbi)
library(org.Hs.eg.db)
###### Survival package
install.packages("lattice")
library(lattice)
install.packages("survival")
library(survival)

######## plot KM curves
install.packages("reshape2")
library(reshape2)
install.packages("data.table")
library(data.table)
install.packages("zoo")
library("zoo")
install.packages("survminer")
library("survminer")
library("survival")

#### ROC package
install.packages("gplots")
library(gplots)
install.packages("ROCR")
library(ROCR)
install.packages("dplyr")
library("dplyr")
install.packages("lubridate")
library(lubridate)
install.packages("robustbase")
install.packages("caret")
library(caret)
library(pROC)

install.packages("timeROC")
install.packages("prodlim")
library(prodlim)
library(quantreg)
install.packages("quantreg")
install.packages("polspline")
library("polspline")
library(timeROC)





###GEO data
# setwd("F:/Clinical Gene expression network Project/Manuscript/Codes/More cancers/Breast cancer GEO/")  


## GSE7390 
Gene1=read.csv("Data/GSE7390_series_matrix.csv")
Gene1_Expression=Gene1[-c(1:113),]
Gene1_Expression=as.matrix(Gene1_Expression)
rownames(Gene1_Expression)=Gene1_Expression[,1]
colnames(Gene1_Expression)=Gene1_Expression[1,]
Gene1_Expression=Gene1_Expression[-c(1,22284),-1]


Clinic1=read.csv("Data/GSE7390_clinical information.csv")
Grade1=Clinic1[,c(13)]
names(Grade1)=Clinic1[,c(3)]

# Grade1=na.omit(Grade1)

Gene1_Expression=Gene1_Expression[,intersect(colnames(Gene1_Expression),Clinic1[,c(3)])]
# Grade1=Grade1[intersect(colnames(Gene1_Expression),names(Grade1))]

ProbeID=read.csv("Data/GPL96_ProbeIDannotation.csv")
ProbeID=ProbeID[-c(1:16),c(1,11)]
ProbeID=as.matrix(ProbeID)
rownames(ProbeID)=ProbeID[,2]

Gene_Name=rownames(ProbeID[match(ProbeID[,1],rownames(Gene1_Expression)),])

rownames(Gene1_Expression)=Gene_Name

# write.csv(Gene1_Expression,file="Data/Gene_GSE7390.csv")
# write.csv(Grade1,file="Data/Grade_GSE7390.csv")

# ### GSE6532
# Data=load ('Data/LUMINAL_GSE6532 RData.RData')
# Gene2=data.untreated
# ProbeID2=annot.untreated
# ProbeID2=ProbeID2[,6]
# 
# Gene2=t(Gene2)
# 
# rownames(Gene2)=ProbeID2[intersect(rownames(Gene2),names(ProbeID2))]
# 
# Clinic2=demo.untreated
# Grade2=Clinic2[,c(7)]
# names(Grade2)=Clinic2[,1]
# 
# Grade2=na.omit(Grade2)
# 
# Gene2=Gene2[,intersect(colnames(Gene2),names(Grade2))]
# 
# 
# 
# ### GSE11121
# # The clinical information can not be found and download 
# 
# ### Data integration
# GeneID=(intersect(rownames(Gene1_Expression),t(rownames(Gene2))))
# GeneID=na.omit(GeneID)
# GeneID=as.vector(GeneID)
# Gene=cbind(Gene1_Expression[GeneID,],Gene2[GeneID,])
# 
# Gene_GSE7390_GSE6532=matrix(as.numeric(Gene),dim(Gene))
# rownames(Gene_GSE7390_GSE6532)=rownames(Gene)
# colnames(Gene_GSE7390_GSE6532)=colnames(Gene)
# 
# sample=union(colnames(Gene1_Expression),colnames(Gene2))
# 
# Grade=c(Grade1,Grade2)
# 
# colnames(Gene_GSE7390_GSE6532)==names(Grade)
# 
# write.csv(Gene_GSE7390_GSE6532,file="Data/Gene_GSE7390_GSE6532.csv")
# write.csv(Grade,file="Data/Grade_GSE7390_GSE6532.csv")

## Survival analysis
### RFS
Time=Clinic1[,c(15)]
Status=Clinic1[,c(16)]
names(Time)=Clinic1[,c(3)]
names(Status)=Clinic1[,c(3)]

### OS
Time=Clinic1[,c(17)]
Status=Clinic1[,c(18)]
names(Time)=Clinic1[,c(3)]
names(Status)=Clinic1[,c(3)]

### DMFS
Time=Clinic1[,c(19)]
Status=Clinic1[,c(20)]
names(Time)=Clinic1[,c(3)]
names(Status)=Clinic1[,c(3)]

Gene1_Expression=Gene1_Expression[,intersect(colnames(Gene1_Expression),names(Time))]
Time=Time[intersect(colnames(Gene1_Expression),names(Time))]
Status=Status[intersect(colnames(Gene1_Expression),names(Status))]

Targets=c('FOXM1','KIF2C', 'SHCBP1', 'CDCA8', 'NCAPG', 'ASPM', 'MELK','STIL', 'MCM10')
Gene_marker=t(matrix(as.numeric(Gene1_Expression[Targets,]),dim(Gene1_Expression[Targets,])))
# Gene_marker=as.numeric(Gene1_Expression['FOXM1',])

# status[time>10*365 & status==1]=0
# time[time>10*365]=10*365


surv=Surv(as.double(Time),as.double(Status))
fit <- coxph(surv ~ (Gene_marker))
A=coef(fit)  #0.05034427

### ROC 
# predicted=colSums (A*t(Gene_marker))  # for Targets
predicted=A*Gene_marker  # for FOXM1 single 
pred <- prediction(predicted,Status)
perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
performance(pred,"auc") # shows calculated AUC for model
dev.new()
plot(perf,colorize=FALSE, col="DarkCyan") # plot ROC curve
lines(c(0,1),c(0,1),col = "gray", lty = 4 )

roc0=roc((Status),predicted)
roc0
AUC=roc0$auc
AUC
## optimal combination
opt <- which.max(rowSums(cbind(roc0$sensitivities,roc0$specificities)))
## optimal cut-off point 
sort(predicted,F)[opt]

####### High and low risk groups

groups=matrix(0,1,length(Status))
groups[predicted<=sort(predicted,F)[opt]]=1
groups[predicted>sort(predicted,F)[opt]]=2
# groups[predicted<=median(predicted)]=1
# groups[predicted>median(predicted)]=2
groups=t(groups)
groups=as.numeric(groups)

# 
# status[time>10*365 & status==1]=0
# time[time>10*365]=10*365

fit <- coxph(surv ~ groups)
S=summary(fit)
waldtest_pvalue=S$waldtest


fit<- survfit(Surv(Time/365, Status) ~ groups, data = as.data.frame(Gene_marker))
# Drawing survival curves
dev.new()
ggsurvplot(fit, legend = c(0.2, 0.2))


# change line size --> 1
# Change line types by groups (i.e. "strata")
# and change color palette
ggsurvplot(fit, data = as.data.frame(Gene_marker), linetype = "strata", pval = TRUE, conf.int = FALSE,
           xlab='Time (Years)',
           risk.table = TRUE, risk.table.y.text.col = TRUE,
           legend = "none",
           legend.title = "Risk",
           legend.labs = c("Low","High"))

S$sctest  # show also Score (logrank) test p value
round(as.numeric(S$waldtest[3]),4)  # show also Wald test p value

###########################################################################
##### Bootstraping for significant test
waldtest_pvalue=matrix(0,1,10000)
HR=matrix(0,1,10000)
for (i in 1:10000)
{
  r=runif(1, min = 1, max = dim(Gene1_Expression)[1])
  r=floor(r)
  Gene_marker=as.numeric(Gene1_Expression[r,])
  
  surv=Surv(as.double(Time),as.double(Status))
  fit <- coxph(surv ~ (Gene_marker))
  A=coef(fit)
  HR[i]=exp(A)
  
  predicted=A*Gene_marker
  pred <- prediction(predicted,Status)
  perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
  # performance(pred,"auc") # shows calculated AUC for model
  # dev.new()
  # plot(perf,colorize=FALSE, col="DarkCyan") # plot ROC curve
  # lines(c(0,1),c(0,1),col = "gray", lty = 4 )
  
  roc0=roc((Status),predicted,quiet=TRUE)
  # roc0
  # AUC=roc0$auc
  # AUC
  ## optimal combination
  opt <- which.max(rowSums(cbind(roc0$sensitivities,roc0$specificities)))
  ## optimal cut-off point 
  # sort(predicted,F)[opt]
  
  ####### High and low risk groups
  
  groups=matrix(0,1,length(Status))
  groups[predicted<=sort(predicted,F)[opt]]=1
  groups[predicted>sort(predicted,F)[opt]]=2
  # groups[predicted<=median(predicted)]=1
  # groups[predicted>median(predicted)]=2
  groups=t(groups)
  groups=as.numeric(groups)
  fit <- coxph(surv ~ groups)
  S=summary(fit)
  waldtest_pvalue[i]=S$waldtest[3]
}

pvalue=na.omit(as.numeric(waldtest_pvalue))
sum(pvalue<0.001819515)/length(pvalue)  # 0.001819515 is pvalue for the biomarker identified

HR=na.omit(as.numeric(HR))
sum(HR>exp(0.1162123))/length(HR)  # 0.001819515 is pvalue for the biomarker identified

########### Plot histogram

# Kernel Density Plot

p_log=-log10(pvalue)
dev.new()
hist(p_log,breaks=100,prob=TRUE, col="#7F7FFF", xlab="-log10(pvalue)", ylab="Probability", main="",xlim = c(0, 4),ylim = c(0, 1),xaxs = "i", yaxs ="i")
lines(density(p_log), col="darkblue", lwd=2)
par(new=TRUE)
lines(rep(-log10(0.001819515),11),0:10,col="#FF0000", lwd=3)

p.value<-sum(pvalue<0.001819515)/length(pvalue)  # 0.001819515 is pvalue for the biomarker identified
p.value
