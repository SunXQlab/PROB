
library(limma)

library(DESeq2)
## Breast cancer 
Gene=read.csv("Data/GSE58626_genes.htseq.csv")
Gene_Expression=as.matrix(Gene)
rownames(Gene_Expression)=Gene[,1]
Gene_Expression=Gene_Expression[,-1]

Gene_Expression=matrix(as.numeric(Gene_Expression),dim(Gene_Expression))
rownames(Gene_Expression)=Gene[,1]

colnames(Gene_Expression)=c(1:12)

batch = c(rep("Control",3),rep("Treatment",9))

Gene_Normalized <- removeBatchEffect(Gene_Expression, batch)

Gene_selected = Gene_Normalized[c('KIF2C', 'SHCBP1', 'CDCA8', 'NCAPG', 'ASPM', 'MELK','STIL', 'MCM10'),]

Gene_selected_0 = Gene_Expression[c('KIF2C', 'SHCBP1', 'CDCA8', 'NCAPG', 'ASPM', 'MELK','STIL', 'MCM10'),]

Gene_Selected_Name = c('KIF2C', 'SHCBP1', 'CDCA8', 'NCAPG', 'ASPM', 'MELK','STIL', 'MCM10')

## DEseq2
countData=Gene_Expression
condition <- factor(c(rep("Control",3),rep("Treatment",9))) # 定义condition
colData <- data.frame(row.names=colnames(countData), condition) # 样品信息矩阵

dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), design= ~ condition ) # 构建dds矩阵
head(dds) # 查看dds矩阵的前6行

dds <- DESeq(dds) # 对dds进行Normalize
resultsNames(dds) # 查看结果的名称
res <- results(dds) # 使用results()函数获取结果，并赋值给res
head (res, n=5) # 查看res矩阵的前5行

mcols(res,use.names= TRUE) # 查看res矩阵每一列的含义
summary(res) # 对res矩阵进行总结

# table(res$padj<0.05) # 取padj小于0.05的数据
res <- res[order(res$padj),] # 按照padj的大小将res重新排列
# diff_gene_deseq2 <- subset(res,padj < 0.05 & (log2FoldChange >0 | log2FoldChange < 0)) # 获取padj小于0.05，表达倍数取以2为对数后绝对值大于1的差异表达基因，赋值给diff_gene_deseq2

diff_gene_deseq2 <- subset(res,padj < 0.05)
intersect(rownames(diff_gene_deseq2),Gene_Selected_Name)



### Wilcox Test for Differential expression analysis
Gene_Expression=Gene_selected_0
pvalue=matrix(0,length(Gene_Selected_Name),1)
Inc_Dec=matrix(0,length(Gene_Selected_Name),1)

for (i in 1:length(Gene_Selected_Name))
{
  G=as.numeric(as.matrix(Gene_Expression[i,4:12]))
  C=as.numeric(as.matrix(Gene_Expression[i,1:3]))
  Inc_Dec[i]=mean(C)-mean(G)
  if (Inc_Dec[i]>0)
  {
    T=wilcox.test(C,G,alternative = 'greater',paired = 0)  #wilcox.test(C-G)
  }
  
  if (Inc_Dec[i]<0)
  {
    T=wilcox.test(C,G,alternative = 'less',paired = 0)  #wilcox.test(C-G)
  }
  
  pvalue[i]=T$p.value
}

which( pvalue<0.05)
Gene_Selected_Name[pvalue<0.05]

## T test
for (i in 1:length(Gene_Selected_Name))
{
  G=as.numeric(as.matrix(Gene_Expression[i,4:12]))
  C=as.numeric(as.matrix(Gene_Expression[i,1:3]))
  Inc_Dec[i]=mean(C)-mean(G)
  if (Inc_Dec[i]>0)
  {
    T=t.test(C,G,alternative = 'greater')  #wilcox.test(C-G)
  }
  
  if (Inc_Dec[i]<0)
  {
    T=t.test(C,G,alternative = 'less')  #wilcox.test(C-G)
  }
  
  pvalue[i]=T$p.value
}
which( pvalue<0.05)
Gene_Selected_Name[pvalue<0.05]

### Validation
KIF2C_C=as.numeric(as.matrix(Gene_Expression[1,1:3]))
KIF2C_K=as.numeric(as.matrix(Gene_Expression[1,4:12]))
KIF2C=data.frame(Samples=(c(rep("Control",3),rep("FOXM1 inhibition",9))),Expression=(c(KIF2C_C,KIF2C_K))/mean(KIF2C_C),Gene=rep("KIF2C",12))
# 

CDCA8_C=as.numeric(as.matrix(Gene_Expression[3,1:3]))
CDCA8_K=as.numeric(as.matrix(Gene_Expression[3,4:12]))
CDCA8=data.frame(Samples=(c(rep("Control",3),rep("FOXM1 inhibition",9))),Expression=(c(CDCA8_C,CDCA8_K))/mean(CDCA8_C),Gene=rep("CDCA8",12))
# 

ASPM_C=as.numeric(as.matrix(Gene_Expression[5,1:3]))
ASPM_K=as.numeric(as.matrix(Gene_Expression[5,4:12]))
ASPM=data.frame(Samples=(c(rep("Control",3),rep("FOXM1 inhibition",9))),Expression=(c(ASPM_C,ASPM_K))/mean(ASPM_C),Gene=rep("ASPM",12))
# 
Data=rbind(KIF2C,CDCA8,ASPM)

### 箱线图
# install.packages("ggpubr")
library(ggpubr)
library(digest)


p<-ggboxplot(ASPM, "Samples", "Expression",
             color = "Samples", palette =c("#00AFBB", "#FC4E07"),  #, "#E7B800"
             linetype = 1,
             title = "ASPM", xlab = FALSE,
             ylab="Normalized expression",
             #error.plot="errorbar",
             font.label = list(size = 50, face = "plain"),
             order = c("Control","FOXM1 inhibition"),
             add="boxplot") #，fill = "Samples",)

dev.new()

p

compar<-list(c("Control","FOXM1 inhibition"))
pp=p+stat_compare_means(comparisons = compar, method = "wilcox.test",  method.args = list(alternative = "greater"))+stat_compare_means(label.y = 11, show.legend = FALSE)   #

ggpar(pp,ylim=c(0,1.6))


# ### Multi-group comparison
# p<-ggboxplot(Data, "Gene", "Expression",
#              color = "Samples", palette =c("#00AFBB", "#FC4E07"),  #, "#E7B800"
#              linetype = 1,
#              title = "KIF2C", xlab = FALSE,
#              ylab="Normalized expression",
#              #error.plot="errorbar",
#              font.label = list(size = 50, face = "plain"),
#              # order = c("Control","FOXM1 inhibition"),
#              add="boxplot") #，fill = "Samples",)
# 
# dev.new()
# 
# p
# 
# compare_means( Expression ~ Samples, data = Data, 
#                group.by = "Gene")
# 
# p + stat_compare_means(label = "p.format")
# 
# 
# p + stat_compare_means(aes(group = Samples), label = "p.format")
# 
# 
# compar<-list(c("Control","FOXM1 inhibition"))
# pp=p+stat_compare_means( method = "wilcox.test",  method.args = list(alternative = "greater"))+stat_compare_means(label.y = 11, show.legend = FALSE)   #
# 
# ggpar(pp,ylim=c(0,1.6))
# 


####### GSE55204
Gene=read.csv("Data/GSE55204.csv",header=FALSE)
Gene=as.matrix(Gene)
Gene_Expression=Gene[-c(1:66),]

Gene_Expression[1:5,1:5]
RowName=Gene_Expression[,1]
Gene_Expression=matrix(as.numeric(Gene_Expression),dim(Gene_Expression))
rownames(Gene_Expression)=RowName
Gene_Expression=Gene_Expression[,-1]


# probeid_geneid<-read.csv("Data/Validation/GPL96-57554.csv") 
# probeid_geneid=probeid_geneid[-(1:16),c(1,11)]
# probeid_geneid=as.matrix(probeid_geneid)
# rownames(probeid_geneid)=probeid_geneid[,2]

### Selected gene expression
Gene_Selected_Name=c('KIF2C', 'SHCBP1', 'CDCA8', 'NCAPG', 'ASPM', 'MELK','STIL', 'MCM10')
ProbeID=probeid_geneid[Gene_Selected_Name,]
Gene_Expression=Gene_Expression[ProbeID[,1],]
# Gene_Expression=Gene_Expression[,-1]
rownames(Gene_Expression)=ProbeID[,2]


####### Wilcox Test
pvalue=matrix(0,length(Gene_Selected_Name),1)
Inc_Dec=matrix(0,length(Gene_Selected_Name),1)

for (i in 1:length(Gene_Selected_Name))
{
  G=as.numeric(as.matrix(Gene_Expression[i,5:7]))
  C=as.numeric(as.matrix(Gene_Expression[i,1:4]))
  Inc_Dec[i]=mean(C)-mean(G)
  if (Inc_Dec[i]>0)
  {
    T=wilcox.test(C,G,alternative = 'greater',paired = 0)  #wilcox.test(C-G)
  }
  
  if (Inc_Dec[i]<0)
  {
    T=wilcox.test(C,G,alternative = 'less',paired = 0)  #wilcox.test(C-G)
  }
  
  pvalue[i]=T$p.value
}

which( pvalue<0.05)
Gene_Selected_Name[pvalue<0.05]


## T test
for (i in 1:length(Gene_Selected_Name))
{
  G=as.numeric(as.matrix(Gene_Expression[i,5:7]))
  C=as.numeric(as.matrix(Gene_Expression[i,1:4]))
  Inc_Dec[i]=mean(C)-mean(G)
  if (Inc_Dec[i]>0)
  {
    T=t.test(C,G,alternative = 'greater')  #wilcox.test(C-G)
  }
  
  if (Inc_Dec[i]<0)
  {
    T=t.test(C,G,alternative = 'less')  #wilcox.test(C-G)
  }
  
  pvalue[i]=T$p.value
}
which( pvalue<0.05)
Gene_Selected_Name[pvalue<0.05]


