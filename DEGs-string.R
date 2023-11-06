
################################################
rm(list=ls())
library(limma)
load("Pancreatic_string_NES.rda")
load("cluster-2-train.rda")
label<-read.table('pred-2-train.txt')
colnames(label)<-'cluster'
rownames(label)<-rownames(label_random)
label$cluster[which(label$cluster==1)]<-2
label$cluster[which(label$cluster==0)]<-1

clu1<-rownames(label)[which(label$cluster==2)]
clu2<-rownames(label)[which(label$cluster==1)]

c1<-string_NES[,clu1]
c2<-string_NES[,clu2]

data<-cbind(c2,c1)

sampleType<-c(rep("Cluster 2",ncol(c1)),rep("Cluster 1",ncol(c2)))
sampleType<-factor(sampleType,levels=c("Cluster 2","Cluster 1"))
design<-model.matrix(~sampleType)


fit<-lmFit(data,design)
fit<-eBayes(fit)
DEGs_result=topTable(fit,coef=2,n=Inf)
save(DEGs_result,file="DEGs_result.rda")

sig_DEGs<-rownames(DEGs_result[which(abs(DEGs_result$logFC)>0.8&DEGs_result$P.Value<0.05),])
save(sig_DEGs,file = "sig_DEGs.rda")

#####################################################
sig_DEGs1<-rownames(DEGs_result[which((DEGs_result$logFC)>0.8&DEGs_result$P.Value<0.05),])
sig_DEGs2<-rownames(DEGs_result[which((DEGs_result$logFC)<(-0.8)&DEGs_result$P.Value<0.05),])



#############################################??
rm(list= ls()) 
library(ggpubr)
library(ggrepel)
library(ggplot2)
library(ggrepel)
load('DEGs_result.rda')

##########################

deg_data<-DEGs_result
head(deg_data)
deg_data$logP<- -log10(deg_data$P.Value)
deg_data$group <- "not-significant"
deg_data$group[which((deg_data$P.Value < 0.05) & (deg_data$logFC > 0.5))] <-"up-regulated"
deg_data$group[which((deg_data$P.Value < 0.05) & (deg_data$logFC < -0.5))] <-"down-regulated"
table(deg_data$group)

deg_data$label<-""
deg_data <- deg_data[order(deg_data$P.Value),]
deg_data$Symbol <- rownames(deg_data)
up_genes <- head(deg_data$Symbol[which(deg_data$group == "up-regulated")],10)
down_genes <- head(deg_data$Symbol[which(deg_data$group == "down-regulated")],10)
deg_top10_genes<-c(as.character(up_genes), as.character(down_genes))
deg_data$label[match(deg_top10_genes,deg_data$Symbol)] <- deg_top10_genes
colnames(deg_data)[8] <- c("Group")


#################################
p1<-ggplot(deg_data,aes(logFC,logP))+
  ggtitle("DEGs for between two subtypes")+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+ # 横向水平参考线
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "#999999")+ # 纵向垂直参考线
  geom_point(aes(size=logP, color=logP))+ # 
  scale_color_gradientn(values = seq(0,1,0.2), colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+ # 指定颜色渐变模式
  scale_size_continuous(range = c(0.5,4))+ # 
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25))+ #
  theme_bw()+ # 
  theme(panel.grid = element_blank(),
        legend.position ="none",
        #legend.position = c(0.81,0.67),
        legend.justification = c(0,1))+ #
  theme(legend.background=element_rect(fill=rgb(1,1,1,alpha=0.001),colour=NA))+ #
  guides(col = guide_colourbar(title = "-Log10_P-value"))+ #
  geom_text_repel(data=deg_data,family="serif",
                  aes(label=label),size=2.5, box.padding=unit(0.25, "lines"),
                  point.padding=unit(0.30,"lines"),
                  max.overlaps = getOption("ggrepel.max.overlaps", default =30))+
  xlab("Log2(Fold change)")+
  ylab("-Log10(P-value)") #

p1
p2<-p1+theme(plot.title = element_text(size=15,family="serif",color="black",face="bold",hjust=0.5),
             axis.title=element_text(size=13,family="serif",color="black",face= "bold"),
             axis.text.y=element_text(size=12,family="serif",color="black",face="bold"),
             axis.text.x=element_text(size=12,family="serif",color="black",face="bold"),
             legend.title=element_text(size=12,family="serif",color="black",face="bold"),
             legend.text=element_text(size=12,family="serif",color="black",face="bold")
)
p2
Pancreatic_edgeR_valcano<-p2
save(Pancreatic_edgeR_valcano,file="Pancreatic_edgeR_valcano.rda")
ggsave("Pancreatic_edgeR_valcano.pdf",width =6.5,height=6)

##############################################################


