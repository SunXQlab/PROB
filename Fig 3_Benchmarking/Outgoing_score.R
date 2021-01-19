
Shalek_data=read.csv('F:/Clinical Gene expression network Project/Reversion/Data/GSE48968-conquer database-0503/GSE48968_combined_new.csv')
LPS_data=Shalek_data[!is.na(Shalek_data[,2]),]
LPS_data=LPS_data[match(unique(LPS_data[,2]),LPS_data[,2]),]
rownames(LPS_data)=LPS_data[,2]
LPS_data=LPS_data[,-c(1,2)]
# LPS_data=LPS_data[,1:2412]

LPS_only= LPS_data[,grep(pattern="LPS",colnames(LPS_data))]
LPS_WT=LPS_only[,c(1:479)]
LPS_WT=na.omit(LPS_WT)

# Antiviral regulators and the targets
#select only the significant genes:
Regulators <- c("STAT2", "STAT1", "HHEX", "RBL1", "Timeless", "FUS")
#set the inferred regulators and the potential targets
# top_group <- c("STAT2", "IRF1", "IRF2")
Targets <- c("Ifnb1", "Cxcl9", "Ifit3", "Tsc22d1", "Tnfsf8", "Isg20", "Nmi", "Iigp1", "Irf7", "Lhx2", "Bbx", "Fus", "Tcf4", "Pml", "Usp12", "Irf7", "Ifit1", "Isg15", "Usp25", "Daxx", "Cd40", "Atm", "Lrf1", "Lgp2", "Mertk", "Cxcl11", "Trim12", "Trim21")

Regulators_Targets_ind=matrix(0,1,length(TFs))

Regulators_Targets=c(Regulators,Targets)

Regulators_Targets=intersect(toupper(Regulators_Targets),toupper(rownames(LPS_WT)))

Regulators_Targets_ind=matrix(0,1,length(Regulators_Targets))
for (i in 1:length(Regulators_Targets)){
  
  IDs=grep(toupper(Regulators_Targets[i]),toupper(rownames(LPS_WT)),value=F)
  for (ii in 1:length(IDs))
  {
    if (toupper(rownames(LPS_WT)[IDs[ii]])==toupper(Regulators_Targets[i]))
    {
      RT_ID=IDs[ii]
    }
    Regulators_Targets_ind[i]=RT_ID
  }
}

write.csv(Regulators_Targets_ind, "F:/Clinical Gene expression network Project/Reversion/Data/Regulators_Targets_ind.csv")


# S_PROB_RT=readMat("S_PROB_RT.mat")
# S_PROB_RT=unlist(S_PROB_RT)
# S_PROB_RT=matrix(S_PROB_RT,29,29)

AM=readMat("AM.mat")
AM=unlist(AM)
Aij=matrix(AM,29,29)

# Act_Inh_strength=readMat("Act_Inh_strength.mat")
# Act_Inh_strength=unlist(Act_Inh_strength)
# Aij=matrix(Act_Inh_strength,29,29)

## Outgoing causality score

# for (jj in 1:(dim(AM)[1]))
# {
#   OSC[jj]=sum(abs(AM[,jj]))
# }

setwd("F:/Clinical Gene expression network Project/Reversion/Codes/Benchmarking_LPS_scRNAseq_data/")
X_stage <- readMat("scExp_Net_Ctime.mat")
X_stage=as.data.frame(X_stage)
Expression=X_stage[1:(dim(X_stage)[1]-1),]
Expression[is.na(Expression),]=0

  ## Outgoing causality score (OCS)
# OCS=matrix(0,1,dim(Aij)[1])
OCS_top=list()
for (jj in 1:length(Regulators))
{
  OCS_top[[Regulators[jj]]]=c()
  for (ii in 1:(dim(Aij)[1]))
  {
    Ctarget=abs(Aij[ii,jj])*(as.numeric(Expression[ii,]))*(as.numeric(Expression[jj,]))  #*S_PROB_RT[ii,jj]
    if (sum(Ctarget)!=0)
    {
      OCS_top[[Regulators[jj]]]=c(OCS_top[[Regulators[jj]]],Ctarget)  #*Expression[ii,]*Expression[jj,]
    }
   
  }
}

OCS_bottom=list()
for (jj in 1:length(Targets))
{
  OCS_bottom[[Targets[jj]]]=c()
  for (ii in 1:(dim(Aij)[1]))
  {
    Ctarget=abs(Aij[ii,jj])*(as.numeric(Expression[ii,]))*(as.numeric(Expression[jj,]))  #*S_PROB_RT[ii,jj]
    if (sum(Ctarget)!=0)
    {
      OCS_bottom[[Targets[jj]]]=c(OCS_bottom[[Targets[jj]]],Ctarget)  #*Expression[ii,]*Expression[jj,]
    }
    
  }
}

# > mean(unlist(OCS_top))

# > 
#   > mean(unlist(OCS_bottom))


  # 
  # t.test(unlist(OCS_top),unlist(OCS_bottom),alternative = "greater")
  # ks.test(unlist(OCS_top),unlist(OCS_bottom),alternative ="less")
  # 

  WT=wilcox.test(unlist(OCS_top),unlist(OCS_bottom),alternative="greater", paired=F)  # Wilcoxon signed rank test with continuity correction;  p-value =0.01625
  # library(vioplot)
  # vioplot(unlist(OCS_top),unlist(OCS_bottom),names=c("Causality from Regulators","Causality from Targets"),col=c("green","brown"),outer = FALSE)
  # title("Outgoing Causality")
  dev.new()
  boxplot(log10(1+unlist(OCS_top)),log10(1+unlist(OCS_bottom)),names=c("Causality from Regulators","Causality from Targets"),col=c("#B3DE69","#80B1D3"),outline=F,boxwex = 0.4,boxwex=0.5,width=c(0.5,0.5))
  title("Outgoing Causality")
  temp <- locator(1) # 在图表上，你喜欢的地方点击一下，文字就出来了
  text(temp, paste("p-value =",WT$p.value))

  # dev.new()
  # vioplot(log10(1+unlist(OCS_top)),log10(1+unlist(OCS_bottom)),names=c("Causality from Regulators","Causality from Targets"),col=c("#B3DE69","#80B1D3"),outer = F,horizontal=FALSE,
  #         xlab="", ylab="Out-going causality score")
  # 
  