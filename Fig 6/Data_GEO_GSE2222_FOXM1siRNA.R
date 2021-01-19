###GEO data
setwd("F:/Clinical Gene expression network Project/Manuscript/Codes/Fig 4-5_Realistic applications/Breast cancer GEO_Fig 5/")  


## Breast cancer 
Gene=read.csv("Data/Validation/GSE2222_series_matrix.csv")
Gene_Expression=Gene[-c(1:49),]
Sample=Gene[24,-1]
Sample=(as.matrix(Sample))

Sample=as.character(Sample)
# Sample=gsub("mutation: ", "", Sample, fixed = TRUE)

rownames(Gene_Expression)=Gene[-c(1:49),1]
colnames(Gene_Expression)=Sample



### ProbeID <--> GeneSymbol
probeid_geneid<-read.csv("Data/Validation/GPL96-57554.csv") 
probeid_geneid=probeid_geneid[-(1:16),c(1,11)]
probeid_geneid=as.matrix(probeid_geneid)
rownames(probeid_geneid)=probeid_geneid[,2]

### Selected gene expression
Gene_Selected_Name=c('SHCBP1',"NCAPG","CDCA8",'MELK','ASPM','MCM10','KIF2C','FOXM1')
ProbeID=probeid_geneid[Gene_Selected_Name,]
Gene_Expression=Gene_Expression[ProbeID[,1],]
Gene_Expression=Gene_Expression[,-1]
rownames(Gene_Expression)=ProbeID[,2]


# write.csv(Gene_Expression,file="Data/Validation/Gene_Selected_GSE2222.csv")

### Significant differential expression (BRCA1 mutation vs normal)
# Gene_Selected_siRNA_control=read.csv("Breast Cancer/Gene_Selected_siRNA_control.csv")
pvalue=matrix(0,length(Gene_Selected_Name),1)
Inc_Dec=matrix(0,length(Gene_Selected_Name),1)

for (i in 1:length(Gene_Selected_Name))
{
  G=as.numeric(as.matrix(Gene_Expression[i,7:9]))
  C=as.numeric(as.matrix(Gene_Expression[i,1:6]))
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

## T test
for (i in 1:length(Gene_Selected_Name))
{
  G=as.numeric(as.matrix(Gene_Expression[i,7:9]))
  C=as.numeric(as.matrix(Gene_Expression[i,1:6]))
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

write.csv(pvalue,file="F:/Clinical Gene expression network Project/Manuscript/Codes/Fig 4_5/Breast Cancer/pvalue_GeneExpression_siRNA_control.csv")
write.csv(Inc_Dec,file="F:/Clinical Gene expression network Project/Manuscript/Codes/Fig 4_5/Breast Cancer/Increase_Decrease_GeneExpression_siRNA_control.csv")


### Validation

CDCA8_C=as.numeric(as.matrix(Gene_Expression[3,1:6]))
CDCA8_K=as.numeric(as.matrix(Gene_Expression[3,7:9]))
CDCA8=data.frame(Samples=(c(rep("Control",6),rep("FOXM1 siRNA",3))),Expression=(c(CDCA8_C,CDCA8_K)))
# 

ASPM_C=as.numeric(as.matrix(Gene_Expression[5,1:6]))
ASPM_K=as.numeric(as.matrix(Gene_Expression[5,7:9]))
ASPM=data.frame(Samples=(c(rep("Control",6),rep("FOXM1 siRNA",3))),Expression=(c(ASPM_C,ASPM_K)))
# 

KIF2C_C=as.numeric(as.matrix(Gene_Expression[7,1:6]))
KIF2C_K=as.numeric(as.matrix(Gene_Expression[7,7:9]))
KIF2C=data.frame(Samples=(c(rep("Control",6),rep("FOXM1 siRNA",3))),Expression=(c(KIF2C_C,KIF2C_K)))
# 
### 箱线图
# install.packages("ggpubr")
library(ggpubr)
library(digest)


p<-ggboxplot(ASPM, "Samples", "Expression",
             color = "Samples", palette =c("#00AFBB", "#FC4E07"),  #, "#E7B800"
             linetype = 1,
             title = "ASPM", xlab = FALSE,
             #error.plot="errorbar",
             font.label = list(size = 50, face = "plain"),
             order = c("Control","FOXM1 siRNA"),
             add="boxplot") #，fill = "Samples",)

dev.new()

p

compar<-list(c("Control","FOXM1 siRNA"))
p+stat_compare_means(comparisons = compar, method = "wilcox.test",  method.args = list(alternative = "greater"))+stat_compare_means(label.y = 11, show.legend = FALSE)   #

# ####  ErrorBar plot  +  p  value
# Control <- ASPM_C
# Case <- ASPM_K
# i=5
# m <- c(mean(Control), mean(Case))
# s <- c(sd(Control), sd(Case))
# d <- data.frame(V = c("Control", "FOXM1 siRNA"), mean = m, sd = s)
# d$V <- factor(d$V, levels = c("Control", "FOXM1 siRNA"))
# 
# # library(ggplot2)
# dev.new()
# mm=max(mean(Control)+sd(Control),mean(Case)+sd(Case))
# mm=mm*(1+0.1)
# p <- ggplot(d, aes(V, mean, fill = V, width = 0.2))+
#     opts(panel.grid.major = none, panel.grid.minor = none) + opts(panel.border = none)
# p <- p + geom_errorbar(aes(ymin = mean, ymax = mean + sd, width = 0.2),
#                        position = position_dodge(width = 0.8)) +
#   geom_bar(stat = "identity", position = position_dodge(width = 0.8),
#            colour = "black") +
#   scale_fill_manual(values = c("#00AFBB", "#FC4E07")) +
#   labs(x = "",y = "Esxpression levels", title = Gene_Selected_Name[i])+
#   
#   theme_bw() + theme(plot.title = element_text(size = 20,face = "bold", vjust = 0.5, hjust = 0.5),
#                      legend.position = "none") + xlab("") + # ylab("") +
#   
#   theme(axis.text.x = element_text(face = "bold", size = 12),
#         axis.text.y = element_text(face = "bold", size = 12),
#         panel.grid = element_blank()) +
#   
#   geom_segment(aes(x = 1, y = mm*0.93, xend = 1, yend = mm*0.95)) +
#   geom_segment(aes(x = 2, y = mm*0.93, xend = 2, yend = mm*0.95)) +
#   geom_segment(aes(x = 1, y = mm*0.95, xend = 1.45, yend = mm*0.95)) + #与源代码有差异，使得“*”在连线中间
#   geom_segment(aes(x = 1.55, y = mm*0.95, xend = 2, yend = mm*0.95)) +
#   annotate("text", x = 1.5, y = mm*0.95, label = "*")+
#   annotate("text", x = 1.5, y = mm*1, label = paste0("p = ",round(as.numeric(pvalue[i]),4)))+
#   ylim(0, mm*1.05)+
#   scale_y_continuous(expand = c(0,0))
# 
# plot(p)
# 
# 
# ### 小提琴图
# # install.packages("ggsignif")
# library(ggsignif)
# # devtools::install_github('cttobin/ggthemr')
# library(ggthemr)
# 
# ggthemr("flat")
# compared = list(c("Control", "FOXM1 siRNA"))
# p <- ggplot(d, aes(V, mean, fill = V, width = 0.2)) +
#   geom_violin() +
#   ylim(0, mm*1.05)+
#   geom_signif(comparisons = compared,
#               step_increase = 0.5,
#               map_signif_level = F,
#               test = wilcox.test)
# p
