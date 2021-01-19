######## Validation of DiffNet Genes' prognostic power
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

##################  ensemble ID  to gene symbol   #################################
ensembl2gene <- toTable(org.Hs.egENSEMBL2EG)
gene2symbol <- toTable(org.Hs.egSYMBOL)
ensemble2symbol <- merge(ensembl2gene, gene2symbol, by = 'gene_id')[2:3]
ensemble2symbol=as.matrix(ensemble2symbol)
rownames(ensemble2symbol)=ensemble2symbol[,2]

###TCGA data (BRCA)

setwd("F:/Clinical Gene expression network Project/Code/TCGA Data/")  


## Breast cancer 
Gene=read.csv("Breast Cancer/RNAseq_TCGA_BRCA.csv")
rownames(Gene)=Gene[,1]
colnames(Gene)=gsub(".", "-", colnames(Gene), fixed = TRUE)
Gene=Gene[,-1]
dim(Gene)  #20530  1218

Survival=read.csv("Breast Cancer/Phenotype_TCGA_BRCA.csv")
Survival=Survival[,c(1,158,159,169)]

Sample_intersect=intersect(strtrim(Survival[,1],15),colnames(Gene))

Survival=Survival[match(Sample_intersect,strtrim(Survival[,1],15)),]

rownames(Survival)=Sample_intersect
Gene=Gene[,Sample_intersect]

Survival[is.na(Survival[,2]),2]=0
Survival[is.na(Survival[,3]),3]=0

time=Survival[,2]+Survival[,3]
status=as.numeric(Survival[,4]=='dead')


######################## Select genes for signature   #####################################

Gene_marker=as.numeric(Gene['FOXM1',])


# status[time>10*365 & status==1]=0
# time[time>10*365]=10*365


surv=Surv(as.double(time),as.double(status))
fit <- coxph(surv ~ (Gene_marker))
A=coef(fit)

### ROC 
predicted=A*Gene_marker
pred <- prediction(predicted,status)
perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
performance(pred,"auc") # shows calculated AUC for model
dev.new()
plot(perf,colorize=FALSE, col="DarkCyan") # plot ROC curve
lines(c(0,1),c(0,1),col = "gray", lty = 4 )

roc0=roc((status),predicted)
roc0
AUC=roc0$auc
AUC
## optimal combination
opt <- which.max(rowSums(cbind(roc0$sensitivities,roc0$specificities)))
## optimal cut-off point 
sort(predicted,F)[opt]

####### High and low risk groups

groups=matrix(0,1,length(status))
groups[predicted<=sort(predicted,F)[opt]]=2
groups[predicted>sort(predicted,F)[opt]]=1
# groups[predicted<=median(predicted)]=1
# groups[predicted>median(predicted)]=2
groups=t(groups)
groups=as.numeric(groups)

# 
# status[time>10*365 & status==1]=0
# time[time>10*365]=10*365
fit<- survfit(Surv(time/365, status) ~ groups, data = as.data.frame(Gene_marker))
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


