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





# ### multiple datasets
Gene1=read.csv("Data/Gene Expression_1089 patients.csv")
Gene1_Expression=as.matrix(Gene1)
rownames(Gene1_Expression)=Gene1_Expression[,1]
# colnames(Gene1_Expression)=Gene1_Expression[1,]
Gene1_Expression=Gene1_Expression[,-1]


Clinic0=read.csv("Data/Clinical Information_1089 patients.csv")

ProbeID=read.csv("Data/GPL96_ProbeIDannotation.csv")
ProbeID=ProbeID[-c(1:16),c(1,11)]
ProbeID=as.matrix(ProbeID)
rownames(ProbeID)=ProbeID[,2]

Gene_Name=rownames(ProbeID[!is.na(match(ProbeID[,1],rownames(Gene1_Expression))),])
rownames(Gene1_Expression)=Gene_Name

## Survival analysis
DatasetID=as.matrix(unique(Clinic0[,2]))
Sig_bootstraping=matrix(NA,1,length(DatasetID))
logrank_P=matrix(NA,1,length(DatasetID))

DatasetID_new=DatasetID[-c(5,8,12,13)]
for (GSEID in DatasetID_new)  # the 5-th dataset only has 1 sample
{
  Clinic1=Clinic0[Clinic0[,2]==GSEID,]
  

### RFS
Time=Clinic1[,c(4)]
Status=Clinic1[,c(3)]
names(Time)=Clinic1[,c(1)]
names(Status)=Clinic1[,c(1)]

# ### OS
# Time=Clinic1[,c(6)]
# Status=Clinic1[,c(5)]
# names(Time)=Clinic1[,c(1)]
# names(Status)=Clinic1[,c(1)]
# 
# ### DMFS
# Time=Clinic1[,c(8)]
# Status=Clinic1[,c(7)]
# names(Time)=Clinic1[,c(1)]
# names(Status)=Clinic1[,c(1)]

### Model 
Gene_Expression=Gene1_Expression[,intersect(colnames(Gene1_Expression),names(Time))]
rownames(Gene_Expression)=rownames(Gene1_Expression)
Time=Time[intersect(colnames(Gene1_Expression),names(Time))]
Status=Status[intersect(colnames(Gene1_Expression),names(Status))]

Targets=c('KIF2C', 'SHCBP1', 'CDCA8', 'NCAPG', 'ASPM', 'MELK','STIL', 'MCM10')
Gene_marker=t(matrix(as.numeric(Gene_Expression[Targets,]),dim(Gene_Expression[Targets,])))
FOXM1=as.numeric(Gene_Expression['FOXM1',])
aij=c(1.5525,3.5763,5.7171,1.6165,2.1345,3.7394,2.0968,3.787)
Gene_marker=(Gene_marker*FOXM1*aij)  #rowSums

# entropy score
# Int_force=Gene_marker*FOXM1*aij
# prob=abs(Int_force)/sum(abs(Int_force))
# S=rowSums(prob*log(prob))


# status[time>10*365 & status==1]=0
# time[time>10*365]=10*365

### COX model
surv=Surv(as.double(Time),as.double(Status))
fit <- coxph(surv ~ (Gene_marker))
A=coef(fit)  #0.05034427
# A0=A

### ROC 
predicted=colSums (A*t(Gene_marker))  # for Targets
# predicted=A*Gene_marker  # for FOXM1 single 
# predicted=S  # for entropy score

# pred <- prediction(predicted,Status)
# perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
# performance(pred,"auc") # shows calculated AUC for model
# dev.new()
# plot(perf,colorize=FALSE, col="DarkCyan") # plot ROC curve
# lines(c(0,1),c(0,1),col = "gray", lty = 4 )
# title(main=GSEID)
        
roc0=roc((Status),predicted)
roc0
AUC=roc0$auc
AUC
dev.new()
plot(roc0,col="DarkCyan")
text(0.4,0.2,paste0('AUC=',round(AUC,4)))
# dev.off()
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
surv=Surv(as.double(Time),as.double(Status))
fit <- coxph(surv ~ groups)
S=summary(fit)
# waldtest_pvalue=S$waldtest
logrank_P[which(DatasetID==GSEID)]=S$sctest[3]  # show also Score (logrank) test p value
# round(as.numeric(S$waldtest[3]),4)  # show also Wald test p value
p_value0=S$waldtest[3]

rm(fit)
fit<- survfit(Surv(Time, Status) ~ groups, data = as.data.frame(Gene_marker))
# Drawing survival curves
dev.new()
ggsurvplot(fit, legend = c(0.2, 0.2))

ggsurvplot(fit, data = as.data.frame(Gene_marker), linetype = "strata", pval = TRUE, conf.int = FALSE,
           xlab='Time (Years)',title=GSEID,
           risk.table = TRUE, risk.table.y.text.col = TRUE,
           legend = "none",
           legend.title = "Risk",
           legend.labs = c("Low","High"))

# dev.off()


###########################################################################
##### Bootstraping for significant test
waldtest_pvalue=matrix(0,1,10000)
HR=matrix(0,1,10000)
for (i in 1:10000)
{
  r=runif(8, min = 1, max = dim(Gene_Expression)[1])
  r=floor(r)
  Gene_marker=matrix(as.numeric(Gene_Expression[r,]),dim(Gene_Expression[r,]))
  # r=runif(1, min = 1, max = dim(Gene_Expression)[1])
  # r=floor(r)
  # TF=as.numeric(Gene_Expression[r,])
  # Gene_marker=(Gene_marker*TF)
  Gene_marker=(Gene_marker*FOXM1)
                      
  surv=Surv(as.double(Time),as.double(Status))
  fit <- coxph(surv ~ t(Gene_marker))
  A=coef(fit)  
  A[is.na(A)]=0
  # HR[i]=exp(A)
  

  
  predicted=colSums(A*Gene_marker)
  # pred <- prediction(predicted,Status)
  # perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
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
  waldtest_pvalue[i]=S$waldtest[3]  #S$sctest[3]
}

pvalue=na.omit(as.numeric(waldtest_pvalue))
Sig_bootstraping[which(DatasetID==GSEID)]=sum(pvalue<p_value0)/length(pvalue)  # 0.001819515 is pvalue for the biomarker identified
# 
# HR=na.omit(as.numeric(HR))
# sum(HR>exp(0.1162123))/length(HR)  # 0.001819515 is pvalue for the biomarker identified


} # this } correspond to ID   


########### Plot histogram

# Kernel Density Plot

p_log=-log10(pvalue)
dev.new()
hist(p_log,breaks=100,prob=TRUE, col="#7F7FFF", xlab="-log10(Wald p value)", ylab="Probability", main="",xlim = c(0, 5),ylim = c(0, 1),xaxs = "i", yaxs ="i")
lines(density(p_log), col="darkblue", lwd=2)
par(new=TRUE)
lines(rep(-log10(p_value0),11),0:10,col="#FF0000", lwd=3)

p.value<-sum(pvalue<p_value0)/length(pvalue)  # 0.001819515 is pvalue for the biomarker identified
p.value
text(4,0.8,paste0(paste0('Prob(p<'),format(p_value0,scientific=TRUE),')=',round(p.value,4)))
title(GSEID)
