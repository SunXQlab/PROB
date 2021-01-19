
###TCGA data
setwd("F:/Clinical Gene expression network Project/Code/TCGA Data/")  


## Breast cancer 
Gene=read.csv("Breast Cancer/RNAseq_TCGA_BRCA.csv")
rownames(Gene)=Gene[,1]
colnames(Gene)=gsub(".", "-", colnames(Gene), fixed = TRUE)
Gene=Gene[,-1]
dim(Gene)  #20530  1218

Phenotype=read.csv("Breast Cancer/Phenotype_TCGA_BRCA.csv")
Phenotype=Phenotype[,c(1,168)]

rownames_Phenotype=Phenotype[,1]
rownames_Phenotype=strtrim(rownames_Phenotype,15)
Phenotype=Phenotype[,-1]
names(Phenotype)=rownames_Phenotype
length(Phenotype)  #1283
summary(Phenotype)



Phenotype=as.matrix(Phenotype)
Phenotype=gsub("stage", "", Phenotype, fixed = TRUE)
Phenotype=gsub("iv", "4", Phenotype, fixed = TRUE)
Phenotype=gsub("iii", "3", Phenotype, fixed = TRUE)
Phenotype=gsub("ii", "2", Phenotype, fixed = TRUE)
Phenotype=gsub("i", "1", Phenotype, fixed = TRUE)
Phenotype=strtrim(Phenotype,2)
Phenotype=as.numeric(Phenotype)
names(Phenotype)=rownames_Phenotype
# Phenotype=gsub(c("x","not reported"), "NA", Phenotype, fixed = TRUE)


# Phenotype[Phenotype==c('iv'),]=4
# Phenotype[Phenotype=='stage iii' || Phenotype=='stage iiia'|| Phenotype=='stage iiib'||Phenotype=='stage iiic',]=3
# Phenotype[Phenotype=='stage ii' || Phenotype=='stage iia'|| Phenotype=='stage iib',]=5
# Phenotype[Phenotype=='stage i' || Phenotype=='stage ia'|| Phenotype=='stage ib',]=1
# Phenotype[Phenotype=='stage x' || Phenotype=='not reported',]=NA


length(Phenotype)  #1283
summary(Phenotype)

Cases=intersect(colnames(Gene),names(na.omit(Phenotype)))
length(Cases)  #1217

Gene=Gene[,Cases]
Phenotype=Phenotype[Cases]


dim(Gene)
length(Phenotype)


write.csv(Gene,file="Breast Cancer/Gene_TCGA_BRCA_processed.csv")
write.csv(Phenotype,file="Breast Cancer/Phenotype_TCGA_BRCA_processed.csv")


#### Meneloma
# 
# source("http://bioc.ism.ac.jp/biocLite.R")
# biocLite("AnnotationDbi")
# install.packages("AnnotationDbi")
# install.packages("bit")
# biocLite("org.Hs.eg.db")
library(bit)
library(AnnotationDbi)
library(org.Hs.eg.db)


ensembl2gene <- toTable(org.Hs.egENSEMBL2EG)
gene2symbol <- toTable(org.Hs.egSYMBOL)
ensemble2symbol <- merge(ensembl2gene, gene2symbol, by = 'gene_id')[2:3]
ensemble2symbol=as.matrix(ensemble2symbol)
rownames(ensemble2symbol)=ensemble2symbol[,1]


Gene=read.csv("Melanoma/RNAseq_TCGA_SKCM.csv")
rownames(Gene)=Gene[,1]
colnames(Gene)=gsub(".", "-", colnames(Gene), fixed = TRUE)
Gene=Gene[,-1]
dim(Gene)  #60483   472


G=matrix(unlist(strsplit(rownames(Gene), ".", fixed = TRUE)),2,length(rownames(Gene)))
GeneID=G[1,]

G=ensemble2symbol[intersect(GeneID,ensemble2symbol[,1]),2]
Gene=Gene[names(G),]
Gene=as.matrix(Gene)
rownames(Gene)=G

Phenotype=read.csv("Melanoma/Phenotype_TCGA_SKCM.csv")
Phenotype=Phenotype[,c(1,112)]

rownames_Phenotype=Phenotype[,1]
# rownames_Phenotype=strtrim(rownames_Phenotype,15)
Phenotype=Phenotype[,-1]
names(Phenotype)=rownames_Phenotype
length(Phenotype)  #477
summary(Phenotype)



Phenotype=as.matrix(Phenotype)
Phenotype=gsub("stage ", "", Phenotype, fixed = TRUE)
Phenotype=gsub("iv", "4", Phenotype, fixed = TRUE)
Phenotype=gsub("iii", "3", Phenotype, fixed = TRUE)
Phenotype=gsub("ii", "2", Phenotype, fixed = TRUE)
Phenotype=gsub("i", "1", Phenotype, fixed = TRUE)
Phenotype=strtrim(Phenotype,1)
Phenotype=as.numeric(Phenotype)
names(Phenotype)=rownames_Phenotype


length(Phenotype)  #477
summary(Phenotype)

Cases=intersect(colnames(Gene),names(!is.na(Phenotype[Phenotype==1|Phenotype==2|Phenotype==3|Phenotype==4])))
length(Cases)  #426

Gene=Gene[,Cases]
Phenotype=Phenotype[Cases]


dim(Gene)
length(Phenotype)


write.csv(Gene,file="Melanoma/Gene_TCGA_SKCM_processed.csv")
write.csv(Phenotype,file="Melanoma/Phenotype_TCGA_SKCM_processed.csv")

Gene_Selected_Name=read.csv("Melanoma/GeneName_GEO_Melanoma.csv")
Gene_Selected_Name=as.character(Gene_Selected_Name[1:48,2])
Gene_Selected_Name=intersect(Gene_Selected_Name,rownames(Gene))


Gene_Selected_Index=as.matrix(0,1,length(Gene_Selected_Name))
for (i in 1:length(Gene_Selected_Name))
{
  Gene_Selected_Index[i]=which(rownames(Gene)==Gene_Selected_Name[i])
}

write.csv(Gene_Selected_Index,file="Melanoma/Gene_Selected_Index.csv")
write.csv(Gene_Selected_Name,file="Melanoma/Gene_Selected_Name.csv")  
  

##### COAD 

Gene=read.csv("COAD/RNAseq_COAD_TCGA.csv")
rownames(Gene)=Gene[,1]
colnames(Gene)=gsub(".", "-", colnames(Gene), fixed = TRUE)
Gene=Gene[,-1]
dim(Gene)  #20530  329

Phenotype=read.csv("COAD/COAD_clinicalMatrix.csv")
Phenotype=Phenotype[,c(1,89)]

rownames_Phenotype=Phenotype[,1]
# rownames_Phenotype=strtrim(rownames_Phenotype,15)
rownames_Phenotype=Phenotype[,1]
Phenotype=Phenotype[,-1]
names(Phenotype)=rownames_Phenotype
length(Phenotype)  #1283
summary(Phenotype)



Phenotype=as.matrix(Phenotype)
Phenotype=gsub("Stage ", "", Phenotype, fixed = TRUE)
Phenotype=gsub("IV", "4", Phenotype, fixed = TRUE)
Phenotype=gsub("III", "3", Phenotype, fixed = TRUE)
Phenotype=gsub("II", "2", Phenotype, fixed = TRUE)
Phenotype=gsub("I", "1", Phenotype, fixed = TRUE)
Phenotype=strtrim(Phenotype,1)
Phenotype=as.numeric(Phenotype)
names(Phenotype)=rownames_Phenotype
Phenotype=Phenotype[!is.na(Phenotype)]


length(Phenotype)  #535
summary(Phenotype)

Cases=intersect(colnames(Gene),names(na.omit(Phenotype)))
length(Cases)  #316

Gene=Gene[,Cases]
Phenotype=Phenotype[Cases]


dim(Gene)
length(Phenotype)


write.csv(Gene,file="COAD/Gene_TCGA_COAD_processed.csv")
write.csv(Phenotype,file="COAD/Phenotype_TCGA_COAD_processed.csv")
