
#######################################################
rm(list=ls())
library(survival)
library(ggplot2)
library(survminer)
library(RColorBrewer)
load('cluster-2-train.rda')
load('Sta_tim_Pancreatic_cancer.rda')
label<-read.table('pred-2-train.txt')
colnames(label)<-'cluster'
rownames(label)<-rownames(label_random)
label$cluster[which(label$cluster==1)]<-2
label$cluster[which(label$cluster==0)]<-1


dat<-cbind(Sta_tim_PAN[rownames(label),c(2,1)],label)
colnames(dat)[c(1,2)]<-c("times","status")
dat1<-data.frame(apply(dat, 2, as.numeric))
dat1[,1]<-dat1[,1]/30

fit2 <- survfit(Surv(times, status)~cluster, data = dat1) 
Bcox<-coxph(Surv(times, status)~cluster, data=dat1)
summcph<-summary(Bcox)
summcph
summcph$sctest[3]
mypalette<-rainbow(10)[c(1,4)]
p1<-ggsurvplot(fit2,legend.title="Subtype",
              main="Overall survival",
              legend.labs=c("Cluster 1","Cluster 2"),
              pval=F, conf.int = T,#
              risk.table =T, # Add risk table
              #risk.table.col = "strata", # Change risk table color by groups
              linetype = "strata", # Change line type by groups
              surv.median.line = "hv", # Specify median survival
              ggtheme=theme_survminer(base_size=12,base_family="serif",
                                      font.main=c(14, "bold","black"),
                                      font.submain= c(13, "bold", "black"),
                                      font.caption = c(14, "bold", "black"),
                                      font.x=c(14, "bold", "black"),font.legend = c(14,"bold"),
                                      font.y=c(14, "bold", "black"),
                                      font.tickslab=c(13, "bold")),
              palette =mypalette )

p2<-p1+labs(x="Month")
p3<-p2$plot+ggtitle("Kaplan-Meier Curve for training cohort")+
  theme(plot.title=element_text(size=15,family="serif",color="black",face="bold",hjust=0.1))
p4<-p3+annotate("text", x=7, y=0.25, label="P=2.26e-08",color="black",cex=5)
p5<-p2$table+theme(legend.position="none")
p6<-ggarrange(p4,p5,heights=c(5,1.8),ncol=1,nrow=2, align = "v")
p6
metabolic_cluster_survival_sdcn2<-p6
ggsave("metabolic_survival_sdcn2.pdf",width =7,height = 7)
save(metabolic_cluster_survival_sdcn2,file = "metabolic_survival_sdcn2.rda")

####################################################
