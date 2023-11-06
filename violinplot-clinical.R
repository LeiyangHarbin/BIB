
###########################################
rm(list=ls())
library(reshape2)
load('cluster-2.rda')
label<-read.table('pred-2.txt')
colnames(label)<-'cluster'
rownames(label)<-rownames(label_random)
label$cluster[which(label$cluster==1)]<-2
label$cluster[which(label$cluster==0)]<-1
load("pancreatic cancer.rda")


Pancreatic_cli<-Pancreatic_cancer[rownames(label),]
dat<-Pancreatic_cli[,c("tumor_purity","met_count","met_site_count","tmb","fga","msi_score")]
dat1<-cbind(label,dat)
dat1$cluster[which(dat1$cluster==1)]<-"Cluster 1"
dat1$cluster[which(dat1$cluster==2)]<-"Cluster 2"

plotA<-dat1
plotAA=melt(plotA,id="cluster")
colnames(plotAA)=c("Subtype","Clinical","Expression")
plotAA<-data.frame(plotAA)

#----------------------------------------------------------
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(ggpubr)
#####################################1
geneMat=colnames(plotA)[2:ncol(plotA)]
dataA=plotAA[which(plotAA[,2]==geneMat[1]),]
col<-c("lawngreen", "orange")
Checkpoint<-c(tumor_purity="Tumor purity")

pA<-ggplot(dataA,aes(x=Subtype,y=Expression,fill=Subtype))+
  geom_violin()+
  geom_boxplot(width=0.2,notch=TRUE)+#scale_y_continuous(trans="log10")+
  #geom_jitter(position=position_jitter(width=0.2,height=0.2),
  #            alpha=0.7,shape=21, size =1.0)+
  facet_grid(~Clinical,space="free",scale="free_x",labeller=labeller(Clinical=Checkpoint))+
  scale_fill_manual(values=col,name="")+
  #scale_fill_manual(values=c(brewer.pal(6,"Set1")[c(3,5)]))+
  guides(fill=guide_legend(title="Subtype"))

p1A<-pA+theme(axis.title=element_text(size=15,family="serif",color="black",face= "bold"),
              axis.text.y=element_text(size=11,family="serif",color="black",face="bold"), 
              axis.text.x=element_text(size=11,family="serif",color="black",face="bold"),
              legend.title=element_text(size=11,family="serif",color="black",face="bold"),
              legend.text=element_text(size=12,family="serif",color="black",face="bold"),
)
p2A<-p1A#+ggtitle("")+theme(plot.title=element_text(size=18,family="serif",color="black",face="bold",hjust=0.5))
p3A<-p2A+labs(x="Subtype",y="Tumor purity")+theme(strip.text.x=element_text(size=15,colour="black",face="bold",family="serif"))
p4A<-p3A+stat_compare_means(method="wilcox.test",label.y=95,label.x = 1.3)

p7A<-p4A+theme(panel.background=element_rect(fill="honeydew",colour="honeydew"
                                             ,size=0.5,linetype="solid"),
               panel.grid.major=element_line(size=0.5,linetype='solid',colour="white")
               ,panel.grid.minor=element_line(size=0.25,linetype='solid',colour="white"))
p7A

###################################2
j=2
geneMat=colnames(plotA)[2:ncol(plotA)]
dataA=plotAA[which(plotAA[,2]==geneMat[j]),]
col<-c("lawngreen", "orange")
Checkpoint<-c(met_count="Metastatic count")

pB<-ggplot(dataA,aes(x=Subtype,y=Expression,fill=Subtype))+
  geom_violin()+
  geom_boxplot(width=0.2,notch=TRUE)+
  scale_y_log10()+
  #geom_boxplot(notch=TRUE)+#scale_y_continuous(trans="log10")+
  #geom_jitter(position=position_jitter(width=0.2,height=0.2),
  #           alpha=0.7,shape=21,size =1.0)+
  facet_grid(~Clinical,space="free",scale="free_x",labeller=labeller(Clinical=Checkpoint))+
  scale_fill_manual(values=col,name="")+guides(fill=guide_legend(title="Subtype"))

p1B<-pB+theme(axis.title=element_text(size=15,family="serif",color="black",face= "bold"),
              axis.text.y=element_text(size=11,family="serif",color="black",face="bold"), 
              axis.text.x=element_text(size=11,family="serif",color="black",face="bold"),
              legend.title=element_text(size=11,family="serif",color="black",face="bold"),
              legend.text=element_text(size=9,family="serif",color="black",face="bold"),
)
p2B<-p1B#+ggtitle("")+theme(plot.title=element_text(size=18,family="serif",color="black",face="bold",hjust=0.5))
p3B<-p2B+labs(x="Subtype",y="Met count")+theme(strip.text.x=element_text(size=15,colour="black",face="bold",family="serif"))
p4B<-p3B+stat_compare_means(method="wilcox.test",label.x=1.3,label.y=1.7)

p7B<-p4B+theme(panel.background=element_rect(fill="honeydew",colour="honeydew"
                                             ,size=0.5,linetype="solid"),
               panel.grid.major=element_line(size=0.5,linetype='solid',colour="white")
               ,panel.grid.minor=element_line(size=0.25,linetype='solid',colour="white"))
p7B
####################################3
j=3
geneMat=colnames(plotA)[2:ncol(plotA)]
dataA=plotAA[which(plotAA[,2]==geneMat[j]),]
col<-c("lawngreen", "orange")
Checkpoint<-c(met_site_count="Metastatic site count")

pC<-ggplot(dataA,aes(x=Subtype,y=Expression,fill=Subtype))+
  geom_violin()+
  scale_y_log10()+
  geom_boxplot(width=0.2,notch=TRUE)+#scale_y_continuous(trans="log10")+
  #geom_jitter(position=position_jitter(width=0.2,height=0.2),
  #           alpha=0.7,shape=21, size =1.0)+
  facet_grid(~Clinical,space="free",scale="free_x",labeller=labeller(Clinical=Checkpoint))+
  scale_fill_manual(values=col,name="")+guides(fill=guide_legend(title="IRPS"))

p1C<-pC+theme(axis.title=element_text(size=15,family="serif",color="black",face= "bold"),
              axis.text.y=element_text(size=11,family="serif",color="black",face="bold"), 
              axis.text.x=element_text(size=11,family="serif",color="black",face="bold"),
              legend.title=element_text(size=11,family="serif",color="black",face="bold"),
              legend.text=element_text(size=9,family="serif",color="black",face="bold"),
)
p2C<-p1C#+ggtitle("")+theme(plot.title=element_text(size=18,family="serif",color="black",face="bold",hjust=0.5))
p3C<-p2C+labs(x="Subtype",y="Met site count")+theme(strip.text.x=element_text(size=15,colour="black",face="bold",family="serif"))
p4C<-p3C+stat_compare_means(method="wilcox.test",label.x=1.3,label.y=1.3)

p7C<-p4C+theme(panel.background=element_rect(fill="honeydew",colour="honeydew"
                                             ,size=0.5,linetype="solid"),
               panel.grid.major=element_line(size=0.5,linetype='solid',colour="white")
               ,panel.grid.minor=element_line(size=0.25,linetype='solid',colour="white"))
p7C
####################################4
j=4
geneMat=colnames(plotA)[2:ncol(plotA)]
dataA=plotAA[which(plotAA[,2]==geneMat[j]),]
col<-c("lawngreen", "orange")
Checkpoint<-c(tmb="TMB")

pD<-ggplot(dataA,aes(x=Subtype,y=Expression,fill=Subtype))+
  geom_violin()+
  scale_y_log10()+
  geom_boxplot(width=0.2,notch=TRUE)+#scale_y_continuous(trans="log10")+
  #geom_jitter(position=position_jitter(width=0.2,height=0.2),
  #            alpha=0.7,shape=21, size =1.0)+
  facet_grid(~Clinical,space="free",scale="free_x",labeller=labeller(Clinical=Checkpoint))+
  scale_fill_manual(values=col,name="")+guides(fill=guide_legend(title="Subtype"))

p1D<-pD+theme(axis.title=element_text(size=15,family="serif",color="black",face= "bold"),
              axis.text.y=element_text(size=11,family="serif",color="black",face="bold"), 
              axis.text.x=element_text(size=11,family="serif",color="black",face="bold"),
              legend.title=element_text(size=11,family="serif",color="black",face="bold"),
              legend.text=element_text(size=9,family="serif",color="black",face="bold"),
)
p2D<-p1D#+ggtitle("")+theme(plot.title=element_text(size=18,family="serif",color="black",face="bold",hjust=0.5))
p3D<-p2D+labs(x="Subtype",y="TMB")+theme(strip.text.x=element_text(size=15,colour="black",face="bold",family="serif"))
p4D<-p3D+stat_compare_means(method="wilcox.test",label.x=1.3,label.y=2.7)

p7D<-p4D+theme(panel.background=element_rect(fill="honeydew",colour="honeydew"
                                             ,size=0.5,linetype="solid"),
               panel.grid.major=element_line(size=0.5,linetype='solid',colour="white")
               ,panel.grid.minor=element_line(size=0.25,linetype='solid',colour="white"))
p7D

####################################4
j=5
geneMat=colnames(plotA)[2:ncol(plotA)]
dataA=plotAA[which(plotAA[,2]==geneMat[j]),]
col<-c("lawngreen", "orange")
Checkpoint<-c(fga="FGA")

pE<-ggplot(dataA,aes(x=Subtype,y=Expression,fill=Subtype))+
  geom_violin()+
  scale_y_log10()+
  geom_boxplot(width=0.2,notch=TRUE)+#scale_y_continuous(trans="log10")+
  #geom_jitter(position=position_jitter(width=0.2,height=0.2),
  #            alpha=0.7,shape=21, size =1.0)+
  facet_grid(~Clinical,space="free",scale="free_x",labeller=labeller(Clinical=Checkpoint))+
  scale_fill_manual(values=col,name="")+guides(fill=guide_legend(title="Subtype"))

p1E<-pE+theme(axis.title=element_text(size=15,family="serif",color="black",face= "bold"),
              axis.text.y=element_text(size=11,family="serif",color="black",face="bold"), 
              axis.text.x=element_text(size=11,family="serif",color="black",face="bold"),
              legend.title=element_text(size=11,family="serif",color="black",face="bold"),
              legend.text=element_text(size=9,family="serif",color="black",face="bold"),
)
p2E<-p1E#+ggtitle("")+theme(plot.title=element_text(size=18,family="serif",color="black",face="bold",hjust=0.5))
p3E<-p2E+labs(x="Subtype",y="FGA")+theme(strip.text.x=element_text(size=15,colour="black",face="bold",family="serif"))
p4E<-p3E+stat_compare_means(method="wilcox.test",label.x=1.3,label.y=0.5)

p7E<-p4E+theme(panel.background=element_rect(fill="honeydew",colour="honeydew"
                                             ,size=0.5,linetype="solid"),
               panel.grid.major=element_line(size=0.5,linetype='solid',colour="white")
               ,panel.grid.minor=element_line(size=0.25,linetype='solid',colour="white"))
p7E

####################################4
j=6
geneMat=colnames(plotA)[2:ncol(plotA)]
dataA=plotAA[which(plotAA[,2]==geneMat[j]),]
col<-c("lawngreen", "orange")
Checkpoint<-c(msi_score="MSIsensor score")

pF<-ggplot(dataA,aes(x=Subtype,y=Expression,fill=Subtype))+
  geom_violin()+
  scale_y_log10()+
  geom_boxplot(width=0.2,notch=TRUE)+#scale_y_continuous(trans="log10")+
  #geom_jitter(position=position_jitter(width=0.2,height=0.2),
  #           alpha=0.7,shape=21, size =1.0)+
  facet_grid(~Clinical,space="free",scale="free_x",labeller=labeller(Clinical=Checkpoint))+
  scale_fill_manual(values=col,name="")+guides(fill=guide_legend(title="Subtype"))

p1F<-pF+theme(axis.title=element_text(size=15,family="serif",color="black",face= "bold"),
              axis.text.y=element_text(size=11,family="serif",color="black",face="bold"), 
              axis.text.x=element_text(size=11,family="serif",color="black",face="bold"),
              legend.title=element_text(size=11,family="serif",color="black",face="bold"),
              legend.text=element_text(size=9,family="serif",color="black",face="bold"),
)
p2F<-p1F#+ggtitle("")+theme(plot.title=element_text(size=18,family="serif",color="black",face="bold",hjust=0.5))
p3F<-p2F+labs(x="Subtype",y="Msi score")+theme(strip.text.x=element_text(size=15,colour="black",face="bold",family="serif"))
p4F<-p3F+stat_compare_means(method="wilcox.test",label.x=1.3,label.y=1.5)

p7F<-p4F+theme(panel.background=element_rect(fill="honeydew",colour="honeydew"
                                             ,size=0.5,linetype="solid"),
               panel.grid.major=element_line(size=0.5,linetype='solid',colour="white")
               ,panel.grid.minor=element_line(size=0.25,linetype='solid',colour="white"))
p7F



figure2<-plot_grid(
  p7A+theme(legend.position="none"),
  p7B+theme(legend.position="none"),
  p7C+theme(legend.position="none"),
  p7D+theme(legend.position="none"),
  p7E+theme(legend.position="none"),
  p7F+theme(legend.position="none"),
  ncol=3)
figure2

cli_violinplot<-figure2
ggsave("Pancreatic-violinplot-clinical-train.pdf",width=10,height=7,units="in")

save(cli_violinplot,file="Pancreatic_cli_violinplot-train.rda")

###########################################

