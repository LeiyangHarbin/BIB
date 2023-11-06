
###############################################
rm(list=ls())
library(factoextra)
library(survival)
library(survminer)
library(cowplot)
load('cluster-2-train.rda')
label<-read.table('pred-2-train.txt')
colnames(label)<-'cluster'
rownames(label)<-rownames(label_random)
label$cluster[which(label$cluster==1)]<-2
label$cluster[which(label$cluster==0)]<-1
load("pancreatic cancer.rda")

Sta_tim_msk<-Pancreatic_cancer[,c(1,27,28)]
Sta_tim_msk$os_status[which(Sta_tim_msk$os_status=="alive")]<-0
Sta_tim_msk$os_status[which(Sta_tim_msk$os_status=="dead")]<-1

dat<-as.data.frame(label)
surdat<-Sta_tim_msk[rownames(dat),c(2,3)]
surdat1<-as.data.frame(apply(surdat,2,as.numeric))
dat1<-cbind(surdat1,dat)
colnames(dat1)<-c("times","status","Cluster")

colon<-dat1
Pancreatic_cancer<-Pancreatic_cancer[rownames(colon),]

colon$Sex<-Pancreatic_cancer$sex
colon$Age<-as.numeric(Pancreatic_cancer$last_contact_age)
colon$Sample_type<-Pancreatic_cancer$sample_type
colon$TMB<-Pancreatic_cancer$tmb
colon$FGA<-Pancreatic_cancer$fga
colon$Msi_score<-Pancreatic_cancer$msi_score
colon$Msi_type<-Pancreatic_cancer$msi_type
colon$Tumor_purity<-Pancreatic_cancer$tumor_purity

colon<-na.omit(colon)

colon[colon$Cluster == 1,3] = "Cluster 1"
colon[colon$Cluster == 2,3] = "Cluster 2"

a<-median(colon$Age)
colon[colon$Age >= a,5] = "Old"
colon[colon$Age < a,5] = "Young"
table(colon$Age)
A1<-which(colon$Age=="8741")
A2<-which(colon$Age=="9807")
colon[c(A1,A2),5]<-"Young"

a<-median(colon$TMB)
colon[colon$TMB >= a,7] = "TMB(high)"
colon[colon$TMB < a,7] = "TMB(low)"

a<-median(colon$FGA)
colon[colon$FGA >= a,8] = "FGA(high)"
colon[colon$FGA < a,8] = "FGA(low)"

a<-median(colon$Msi_score)
colon[colon$Msi_score >= a,9] = "MSIsensor (high)"
colon[colon$Msi_score < a,9] = "MSIsensor (low)"

a<-median(colon$Tumor_purity)
colon[colon$Tumor_purity >= a,11] = "Tumor purity(high)"
colon[colon$Tumor_purity < a,11] = "Tumor purity(low)"
table(colon$Tumor_purity)
colon[colon$Tumor_purity=="5",11]="Tumor purity(low)"

colnames(colon)<-c("times","status","Cluster","Sex","Age","Sample_type",
                   "TMB","FGA","MSIsensor_score","MSIsensor_type","Tumor_purity")

colon1<- within(colon, {
  Sex <- factor(Sex, labels = c("female", "male"))
  Cluster <- factor(Cluster, labels = c("Cluster 1", "Cluster 2"))
  Age <- factor(Age, labels = c("Old", "Young"))
  Sample_type <- factor(Sample_type,labels = c("Primary","Metastasis"))
  TMB <- factor(TMB,labels = c("TMB(high)","TMB(low)"))
  FGA <- factor(FGA,labels = c("FGA(high)","FGA(low)"))
  MSIsensor_score <- factor(MSIsensor_score,labels = c("MSIsensor (high)","MSIsensor (low)"))
  MSIsensor_type <- factor(MSIsensor_type,labels = c("Stable","Indeterminate","Instable","Do not report"))
  Tumor_purity <- factor(Tumor_purity,labels = c("Tumor purity(high)","Tumor purity(low)"))
})

bigmodel <-
  coxph(Surv(times, status) ~ Cluster+Sex+Age+Sample_type+TMB+FGA+MSIsensor_score+MSIsensor_type+Tumor_purity,
        data = colon )

summary(bigmodel)

ggforest(bigmodel,fontsize = 1.0)
ggsave("pancreatic_muticox_forest-train.pdf",width=9,height=8)

#####################################

