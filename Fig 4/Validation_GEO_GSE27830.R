###GEO data
setwd("F:/Clinical Gene expression network Project/Code/TCGA Data/")  


## Breast cancer 
Gene=read.csv("Breast Cancer/GSE27830_series_matrix.csv")
Gene_Expression=Gene[-c(1:61),-1]
Sample=Gene[41,-1]
Sample=(as.matrix(Sample))

Sample=as.character(Sample)
Sample=gsub("mutation status: ", "", Sample, fixed = TRUE)

rownames(Gene_Expression)=Gene[-c(1:61),1]
colnames(Gene_Expression)=Sample

Gene_Normal=Gene_Expression[,colnames(Gene_Expression)=='Non']
Gene_BRCA1_Mutated=Gene_Expression[,colnames(Gene_Expression)=='BRCA1']

# write.csv(Gene_Name,file="Breast Cancer/GeneName_GEO_BreastCancer.csv")

## Gene Index
Gene_TCGA=read.csv("Breast Cancer/Gene_TCGA_BRCA_processed.csv")
Gene_Name=read.csv("Breast Cancer/GeneName_GEO_BreastCancer.csv")
Gene_all=as.character(Gene_TCGA[,1])
Gene_selected=t(as.character(Gene_Name[,2]))
Gene_selected=t(unique(t(Gene_selected)))
Index=data.frame(1,length(Gene_selected))
for (i in 1:length(Gene_selected))
{
  Index[i]=which(Gene_all==Gene_selected[i])
}

Index=as.numeric(Index)
Index=na.omit(Index)

Gene_Selected_Name=Gene_all[Index]


# 
# write.csv(Gene_Selected_Name,file="Breast Cancer/SelectedGeneName_GEO_BreastCancer.csv")
# write.csv(Index,file="Breast Cancer/Gene_selected_Index.csv")
# write.csv(Gene_siRNA_control,file="Breast Cancer/Gene_siRNA_control_GEO_BreastCancer.csv")

### ProbeID <--> GeneSymbol
probeid_geneid<-read.csv("Breast Cancer/GPL570_ProbeID_GeneSymbol_annotation.csv") 
probeid_geneid=probeid_geneid[-(1:27),c(1,3)]
probeid_geneid=as.matrix(probeid_geneid)
rownames(probeid_geneid)=probeid_geneid[,2]

### Selected gene expression
ProbeID=probeid_geneid[Gene_Selected_Name,]

Gene_Normal=Gene_Normal[ProbeID[,1],]
rownames(Gene_Normal)=ProbeID[,2]

Gene_BRCA1_Mutated=Gene_BRCA1_Mutated[ProbeID[,1],]
rownames(Gene_BRCA1_Mutated)=ProbeID[,2]


Gene_Selected_GSE27830=cbind(Gene_Normal,Gene_BRCA1_Mutated)  #,Gene_non_BRCA_mutation
write.csv(Gene_Selected_GSE27830,file="Breast Cancer/Gene_Selected_GSE27830.csv")

### Validation
# 
# GATA3_normal=as.numeric(as.matrix(Gene_Normal[40,]))
# GATA3_BRCA1Mutated=as.numeric(as.matrix(Gene_BRCA1_Mutated[41,]))
# GATA3=data.frame(Samples=(c(rep("Normal",5),rep("BRCA1_Mutated",7))),Expression=(c(GATA3_normal,GATA3_BRCA1Mutated)))
# 
# CTNNB1_normal=as.numeric(as.matrix(Gene_Normal[30,]))
# GATA3_BRCA1Mutated=as.numeric(as.matrix(Gene_BRCA1_Mutated[41,]))
# GATA3=data.frame(Samples=(c(rep("Normal",5),rep("BRCA1_Mutated",7))),Expression=(c(GATA3_normal,GATA3_BRCA1Mutated)))


TOP2A_nomutation=as.numeric(as.matrix(Gene_Normal[71,]))
TOP2A_BRCA1Mutated=as.numeric(as.matrix(Gene_BRCA1_Mutated[71,]))
TOP2A=data.frame(Samples=(c(rep("non mutation",76),rep("BRCA1 mutated",47))),Expression=(c(TOP2A_nomutation,TOP2A_BRCA1Mutated)))

BRCA1_nomutation=as.numeric(as.matrix(Gene_Normal[7,]))
BRCA1_BRCA1Mutated=as.numeric(as.matrix(Gene_BRCA1_Mutated[7,]))
BRCA1=data.frame(Samples=(c(rep("non mutation",76),rep("BRCA1 mutated",47))),Expression=(c(BRCA1_nomutation,BRCA1_BRCA1Mutated)))

# install.packages("ggpubr")
library(ggpubr)
library(digest)


p<-ggboxplot(TOP2A, "Samples", "Expression",
             color = "Samples", palette =c("#00AFBB", "#FC4E07"),  #, "#E7B800"
             linetype = 1,
             title = FALSE, xlab = FALSE,
             font.label = list(size = 50, face = "plain"),
             order = c("non mutation","BRCA1 mutated"),
             # fill = "Samples", palette =c("#FC4E07", "#00AFBB"),
             add = "jitter")

dev.new()

p

compar<-list(c("non mutation","BRCA1 mutated"))


p+stat_compare_means(comparisons = compar, method = "wilcox.test",  method.args = list(alternative = "less"),aes(label = paste0("p = ", ..p.format..)))+stat_compare_means(label.y = 11, show.legend = FALSE)   #

dev.new()
boxplot(boxplot( Expression ~ Samples, TOP2A,  border = c("brown", "green")))


# ### Correlation of TOP2A and BRCA1
# TOP2A=as.numeric(Gene_Expression[probeid_geneid["TOP2A",1],])
# BRCA1=as.numeric(Gene_Expression[probeid_geneid["BRCA1",1],])
# 
# ## 
# TOP2A_nomutation=as.numeric(as.matrix(Gene_Normal[71,]))
# TOP2A_BRCA1Mutated=as.numeric(as.matrix(Gene_BRCA1_Mutated[71,]))
# TOP2A=c(TOP2A_nomutation,TOP2A_BRCA1Mutated)
# 
# BRCA1_nomutation=as.numeric(as.matrix(Gene_Normal[7,]))
# BRCA1_BRCA1Mutated=as.numeric(as.matrix(Gene_BRCA1_Mutated[7,]))
# BRCA1=c(BRCA1_nomutation,BRCA1_BRCA1Mutated)
# 
# CONDITION=c(rep(1,length(TOP2A_nomutation)),rep(0,length(TOP2A_BRCA1Mutated)))
# library("ggm")
# Data=data.frame(top2a=TOP2A,brca1=BRCA1,condition=CONDITION)
# r=pcor.test(pcor(c(1,2,3), var(Data)), 1, n=length(CONDITION))
# R=pcor(c('top2a','brca1','condition'), var(Data))
# 
# # 
# 
# cor.test(TOP2A,BRCA1)
# 
# cor.test(as.numeric(Gene_Selected_GSE27830[7,]),as.numeric(Gene_Selected_GSE27830[71,]),method = "kendall")
# 
# cor.test(BRCA1_nomutation,TOP2A_nomutation)
# 
# cor.test(BRCA1_BRCA1Mutated,TOP2A_BRCA1Mutated)
# 
# cor.test(BRCA1_BRCA1Mutated-mean(BRCA1_nomutation),TOP2A_BRCA1Mutated-mean(TOP2A_nomutation))
# 
# install.packages("entropy")
# library("entropy")
# mi.empirical(cbind(TOP2A,BRCA1))
