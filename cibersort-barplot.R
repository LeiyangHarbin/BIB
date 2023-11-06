
######################################
rm(list = ls())
library(reshape2)
library(ggpmisc)
library(ggpubr)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
load("cluster-2-train.rda")
label<-read.table('pred-2-train.txt')
load("cibersort_raw.rda")
colnames(label)<-'cluster'
rownames(label)<-rownames(label_random)
label$cluster[which(label$cluster==1)]<-2
label$cluster[which(label$cluster==0)]<-1


rownames(cibersort_raw)<-cibersort_raw[,1]

cibersort_raw<-cibersort_raw[rownames(label),]

label<-label[cibersort_raw[,1],]

data<-cbind(cibersort_raw,label)
colnames(data)[24]<-"subtype"
data$subtype[which(data$subtype==1)]<-"Cluster 1"
data$subtype[which(data$subtype==2)]<-"Cluster 2"
data<-data[,-1]

plotA=cbind(rownames(data),data[,c(23,1:22)])
plotAA=melt(plotA,id=c("subtype","rownames(data)"))
colnames(plotAA)=c("Subtype","Patient","Cell_type","Proportion")
plotAA<-data.frame(plotAA)
dataB<-plotAA
save(dataB,file="dataB.rda")

#

mypalette<-colorRampPalette(brewer.pal(8,"Set1"))

#mypalette<-colorRampPalette(brewer.pal(12,"Set3"))
col<-c("cadetblue1", "palevioletred1")
col<-c("cadetblue1", "orange")

#Checkpoint<-c(cluster1="Cluster 1",cluster2="Cluster 2")
#,ncol=2
p<-ggplot(dataB,aes(Patient,Proportion,fill = Cell_type))+
  facet_grid(~Subtype,scales="free",space="free"
  )+ 
  geom_bar(position = "stack",stat = "identity") +
  labs(fill = "Cell Type",x = "Sample",y = "Estiamted Proportion")+theme_bw()+
  theme(axis.text.x = element_blank())+theme(axis.ticks.x=element_blank())+
  scale_y_continuous(expand=c(0.01,0)) +
  scale_fill_manual(values=mypalette(23))
#
p1<-p+ theme(legend.position = "bottom")
p2<-p1+theme(axis.title=element_text(size=15,family="serif",color="black",face= "bold"),
             axis.text.y=element_text(size=11,family="serif",color="black",face="bold"),#y轴刻度内容调整
             legend.title=element_text(size=11,family="serif",color="black",face="bold"),
             legend.text=element_text(size=9,family="serif",color="black",face="bold"),
)
p3<-p2+ggtitle("")+theme(plot.title=element_text(size=18,family="serif",color="black",face="bold",hjust=0.5))
p4<-p3+theme(strip.text.x=element_text(size=15,colour="black",face="bold",family="serif"))
p4

CIB_barplot<-p4

save(CIB_barplot,file="CIB_barplot_string.rda")
ggsave("CIB_barplot_string.pdf",width=12,height=9,units="in")

##########################################################
