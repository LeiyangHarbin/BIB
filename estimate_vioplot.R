
#####################################
rm(list=ls())
load("estimateScore.rda")
load("cluster-2-train.rda")
label<-read.table('pred-2-train.txt')
colnames(label)<-'cluster'
rownames(label)<-rownames(label_random)
label$cluster[which(label$cluster==1)]<-2
label$cluster[which(label$cluster==0)]<-1
esti<-estimateScore[[1]]
sample<-rownames(esti)
bbb1<-cbind(sample,esti)

bb=data.frame(label)
sample<-rownames(bb)
bbb2<-cbind(sample,bb)

MerMat=merge(bbb1,bbb2,by="sample")
data1<-MerMat[,c(2:7)]
colnames(data1)<-c("CYT","Stromal.score",
                   "Immune.score","ESTIMATE.score","Tumor.purity","cluster")

#-------------------------------------------------------------------------
library(reshape2)
library(ggpmisc)
library(ggpubr)
library(cowplot)
library(ggplot2)
library(dplyr)
require(magrittr)
library(forcats)
library(Hmisc)
library(scales)

data1$cluster[which(data1$cluster == 1)]<-"Cluster 1"
data1$cluster[which(data1$cluster == 2)]<-"Cluster 2"
rownames(data1)<-MerMat$sample

plotA=data1[,c(1,6)]
plotAA=melt(plotA,id="cluster")
colnames(plotAA)=c("cluster","Immune.Checkpoint","CYT");
plotAA<-data.frame(plotAA)
dataA=plotAA
dataA[,3]=as.numeric(as.vector(dataA[,3]))

col<-c("#00AFBB", "#FC4E07")
Checkpoint<-c(CYT="CYT")
mode(dataA)
pB<-ggplot(dataA,aes(x=cluster,y=CYT,fill=cluster))+
  geom_violin()+
  geom_boxplot(width=0.2,notch=T)+
  facet_grid(~Immune.Checkpoint,space="free",scale="free_x",labeller=labeller(Immune.Checkpoint=Checkpoint))+
  scale_fill_manual(values=col,name="")+guides(fill=guide_legend(title="cluster"))+scale_y_log10()
p1<-pB+theme(axis.title=element_text(size=15,family="serif",color="black",face= "bold"),
             axis.text.y=element_text(size=12,family="serif",color="black",face="bold"),
             axis.text.x=element_text(size=12,family="serif",color="black",face="bold"),
             legend.title=element_text(size=13,family="serif",color="black",face="bold"),
             legend.text=element_text(size=10,family="serif",color="black",face="bold"))
#p2<-p1+ggtitle("Tumor purity")+theme(plot.title=element_text(size=15,family="serif",color="black",face="bold",hjust=0.5))
p3<-p1+labs(x="Subtype",y="CYT")+theme(strip.text.x=element_text(size=15,colour="black",face="bold",family="serif"))
p4<-p3+stat_compare_means(method="wilcox.test",label.y.npc=0.99,label.x = 1.3)#label.y,label.x调整显示值的位置
p5A<-p4+theme(panel.background = element_rect(fill="snow",colour="snow",size=0.5,linetype="solid"),
              panel.grid.major=element_line(size=0.5,linetype='solid',colour="white"),
              panel.grid.minor=element_line(size=0.25,linetype='solid',colour="white"),legend.position="none")
p5A


plotA=data1[,c(2,6)]
plotAA=melt(plotA,id="cluster")
colnames(plotAA)=c("cluster","Immune.Checkpoint","Stromal.score");
plotAA<-data.frame(plotAA)
dataA=plotAA
dataA[,3]=as.numeric(as.vector(dataA[,3]))

#col<-c("cadetblue1", "palevioletred1")
Checkpoint<-c(Stromal.score="Stromal score")
mode(dataA)
pB<-ggplot(dataA,aes(x=cluster,y=Stromal.score,fill=cluster))+
  geom_violin()+
  geom_boxplot(width=0.2,notch = T)+
  facet_grid(~Immune.Checkpoint,space="free",scale="free_x",labeller=labeller(Immune.Checkpoint=Checkpoint))+
  scale_fill_manual(values=col,name="")+guides(fill=guide_legend(title="cluster"))+ylim(-2000,2000)
p1<-pB+theme(axis.title=element_text(size=15,family="serif",color="black",face= "bold"),
             axis.text.y=element_text(size=12,family="serif",color="black",face="bold"),
             axis.text.x=element_text(size=12,family="serif",color="black",face="bold"),
             legend.title=element_text(size=13,family="serif",color="black",face="bold"),
             legend.text=element_text(size=10,family="serif",color="black",face="bold"))
#p2<-p1+ggtitle("Tumor purity")+theme(plot.title=element_text(size=15,family="serif",color="black",face="bold",hjust=0.5))
p3<-p1+labs(x="Subtype",y="Stromal score")+theme(strip.text.x=element_text(size=15,colour="black",face="bold",family="serif"))
p4<-p3+stat_compare_means(method="wilcox.test",label.y.npc=0.99,label.x = 1.3)#label.y,label.x调整显示值的位置
p5B<-p4+theme(panel.background = element_rect(fill="snow",colour="snow",size=0.5,linetype="solid"),
              panel.grid.major=element_line(size=0.5,linetype='solid',colour="white"),
              panel.grid.minor=element_line(size=0.25,linetype='solid',colour="white"),legend.position="none")
p5B

plotA=data1[,c(3,6)]
plotAA=melt(plotA,id="cluster")
colnames(plotAA)=c("cluster","Immune.Checkpoint","Immune.score");
plotAA<-data.frame(plotAA)
dataA=plotAA
dataA[,3]=as.numeric(as.vector(dataA[,3]))

#col<-c("cadetblue1", "palevioletred1")
Checkpoint<-c(Immune.score="Immune score")
mode(dataA)
pB<-ggplot(dataA,aes(x=cluster,y=Immune.score,fill=cluster))+
  geom_violin()+
  geom_boxplot(width=0.2,notch = T)+
  facet_grid(~Immune.Checkpoint,space="free",scale="free_x",labeller=labeller(Immune.Checkpoint=Checkpoint))+
  scale_fill_manual(values=col,name="")+guides(fill=guide_legend(title="cluster"))+ylim(-1000,4000)
p1<-pB+theme(axis.title=element_text(size=15,family="serif",color="black",face= "bold"),
             axis.text.y=element_text(size=12,family="serif",color="black",face="bold"),
             axis.text.x=element_text(size=12,family="serif",color="black",face="bold"),
             legend.title=element_text(size=13,family="serif",color="black",face="bold"),
             legend.text=element_text(size=10,family="serif",color="black",face="bold"))
#p2<-p1+ggtitle("Tumor purity")+theme(plot.title=element_text(size=15,family="serif",color="black",face="bold",hjust=0.5))
p3<-p1+labs(x="Subtype",y="Immune score")+theme(strip.text.x=element_text(size=15,colour="black",face="bold",family="serif"))
p4<-p3+stat_compare_means(method="wilcox.test",label.y.npc=0.96,label.x = 1.3)#label.y,label.x调整显示值的位置
p5C<-p4+theme(panel.background = element_rect(fill="snow",colour="snow",size=0.5,linetype="solid"),
              panel.grid.major=element_line(size=0.5,linetype='solid',colour="white"),
              panel.grid.minor=element_line(size=0.25,linetype='solid',colour="white"),legend.position="none")
p5C

plotA=data1[,c(4,6)]
plotAA=melt(plotA,id="cluster")
colnames(plotAA)=c("cluster","Immune.Checkpoint","ESTIMATE.score");
plotAA<-data.frame(plotAA)
dataA=plotAA
dataA[,3]=as.numeric(as.vector(dataA[,3]))

#col<-c("cadetblue1", "palevioletred1")
Checkpoint<-c(ESTIMATE.score="ESTIMATE score")
mode(dataA)
pB<-ggplot(dataA,aes(x=cluster,y=ESTIMATE.score,fill=cluster))+
  geom_violin()+
  geom_boxplot(width=0.2,notch = T)+
  facet_grid(~Immune.Checkpoint,space="free",scale="free_x",labeller=labeller(Immune.Checkpoint=Checkpoint))+
  scale_fill_manual(values=col,name="")+guides(fill=guide_legend(title="cluster"))+ylim(-2000,6000)
p1<-pB+theme(axis.title=element_text(size=15,family="serif",color="black",face= "bold"),
             axis.text.y=element_text(size=12,family="serif",color="black",face="bold"),
             axis.text.x=element_text(size=12,family="serif",color="black",face="bold"),
             legend.title=element_text(size=13,family="serif",color="black",face="bold"),
             legend.text=element_text(size=10,family="serif",color="black",face="bold"))
#p2<-p1+ggtitle("Tumor purity")+theme(plot.title=element_text(size=15,family="serif",color="black",face="bold",hjust=0.5))
p3<-p1+labs(x="Subtype",y="ESTIMATE score")+theme(strip.text.x=element_text(size=15,colour="black",face="bold",family="serif"))
p4<-p3+stat_compare_means(method="wilcox.test",label.y.npc=0.99,label.x = 1.3)#label.y,label.x调整显示值的位置
p5D<-p4+theme(panel.background = element_rect(fill="snow",colour="snow",size=0.5,linetype="solid"),
              panel.grid.major=element_line(size=0.5,linetype='solid',colour="white"),
              panel.grid.minor=element_line(size=0.25,linetype='solid',colour="white"),legend.position="none")
p5D

plotA=data1[,c(5,6)]
plotAA=melt(plotA,id="cluster")
colnames(plotAA)=c("cluster","Immune.Checkpoint","Tumor.purity");
plotAA<-data.frame(plotAA)
dataA=plotAA
dataA[,3]=as.numeric(as.vector(dataA[,3]))

#col<-c("cadetblue1", "palevioletred1")
Checkpoint<-c(Tumor.purity="Tumor purity")
mode(dataA)
pB<-ggplot(dataA,aes(x=cluster,y=Tumor.purity,fill=cluster))+
  geom_violin()+
  geom_boxplot(width=0.2,notch = T)+
  facet_grid(~Immune.Checkpoint,space="free",scale="free_x",labeller=labeller(Immune.Checkpoint=Checkpoint))+
  ylim(0,1)+
  scale_fill_manual(values=col,name="")+guides(fill=guide_legend(title="cluster"))
p1<-pB+theme(axis.title=element_text(size=15,family="serif",color="black",face= "bold"),
             axis.text.y=element_text(size=12,family="serif",color="black",face="bold"),
             axis.text.x=element_text(size=12,family="serif",color="black",face="bold"),
             legend.title=element_text(size=13,family="serif",color="black",face="bold"),
             legend.text=element_text(size=10,family="serif",color="black",face="bold"))
#p2<-p1+ggtitle("Tumor purity")+theme(plot.title=element_text(size=15,family="serif",color="black",face="bold",hjust=0.5))
p3<-p1+labs(x="Subtype",y="Tumor purity")+theme(strip.text.x=element_text(size=15,colour="black",face="bold",family="serif"))
p4<-p3+stat_compare_means(method="wilcox.test",label.y=1.0,label.x = 1.3)#label.y,label.x调整显示值的位置
p5E<-p4+theme(panel.background = element_rect(fill="snow",colour="snow",size=0.5,linetype="solid"),
              panel.grid.major=element_line(size=0.5,linetype='solid',colour="white"),
              panel.grid.minor=element_line(size=0.25,linetype='solid',colour="white"),legend.position="none")
p5E

figure1<-plot_grid(p5B,p5C,p5D,p5E,ncol=2)
figure2=plot_grid(figure1+ggtitle("")+
                    theme(plot.title=element_text(size=18,family="serif",color="black",face="bold",hjust=0.5)))
Estimate_vio_train<-figure1
Estimate_vio_train
save(Estimate_vio_train,file="Estimate_vio_train.rda")
ggsave("Estimate_vio_train.pdf",width=8,height=8,units="in")


###########################
CYT_vio_train<-p5A
CYT_vio_train
save(CYT_vio_train,file="CYT_vio_train.rda")
ggsave("CYT_vio_train.pdf",width=5,height=5,units="in")

