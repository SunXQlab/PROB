
Shalek_data=read.csv('F:/Clinical Gene expression network Project/Reversion/Data/GSE48968-conquer database-0503/GSE48968_combined_new.csv')
LPS_data=Shalek_data[!is.na(Shalek_data[,2]),]
LPS_data=LPS_data[match(unique(LPS_data[,2]),LPS_data[,2]),]
rownames(LPS_data)=LPS_data[,2]
LPS_data=LPS_data[,-c(1,2)]
# LPS_data=LPS_data[,1:2412]

LPS_only= LPS_data[,grep(pattern="LPS",colnames(LPS_data))]
LPS_WT=LPS_only[,c(1:479)]
LPS_WT=na.omit(LPS_WT)

# PAM_only= LPS_data[,grep(pattern="PAM",colnames(LPS_data))]
# PAM_WT=PAM_only[,c(1:341)]

Capture_time=as.numeric(substring(colnames(LPS_WT),5,5))

write.csv(LPS_WT, "F:/Clinical Gene expression network Project/Reversion/Data/LPS_WT_scExpression.csv")
write.csv(Capture_time, "F:/Clinical Gene expression network Project/Reversion/Data/LPS_Capture_time.csv")


### processed data
Shalek_data=read.csv('F:/Clinical Gene expression network Project/Reversion/Codes/Benchmarking_LPS_scRNAseq_data/trajectory inference/data/LPS_tpm.csv')
LPS_data=Shalek_data[!is.na(Shalek_data[,1]),]
LPS_data=LPS_data[match(unique(LPS_data[,1]),LPS_data[,1]),]
rownames(LPS_data)=LPS_data[,1]
LPS_data=LPS_data[,-c(1)]
LPS_WT=LPS_data

Capture_time=as.numeric(substring(colnames(LPS_WT),5,5))

write.csv(LPS_WT, "F:/Clinical Gene expression network Project/Reversion/Data/LPS_WT_scExpression.csv")
write.csv(Capture_time, "F:/Clinical Gene expression network Project/Reversion/Data/LPS_Capture_time.csv")

### 
TF_data=read.csv('F:/Clinical Gene expression network Project/Reversion/Data/TF-TF Coverage Score.csv')
TFs=colnames(TF_data)
TFs=TFs[-1]
TFs=toupper(TFs)
rownames(TF_data)=TF_data[,1]
TF_data=TF_data[,-1]
TF_network=matrix(as.numeric(TF_data>0.3),dim(TF_data))

write.csv(TF_network, "F:/Clinical Gene expression network Project/Reversion/Data/TF_network.csv")


# f<-function(index,x){
#   a<-c()
#   for (i in index){
#     a[i]<-grep(i,x,value=F)
#   }
#   return(a)
# }

TF_ind=matrix(0,1,length(TFs))

for (i in 1:length(TFs))
{
  IDs=grep(TFs[i],toupper(rownames(LPS_WT)),value=F)
  for (ii in 1:length(IDs))
  {
    if (toupper(rownames(LPS_WT)[IDs[ii]])==TFs[i])
    {
      TF_ID=IDs[ii]
    }
  TF_ind[i]=TF_ID
  }

 
  # TF_ind[i]=IDs
}

write.csv(TF_ind, "F:/Clinical Gene expression network Project/Reversion/Data/TF_network_ind.csv")



  