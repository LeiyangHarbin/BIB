
################################################
rm(list=ls())
library(rms)
library(survival)
library(survminer)
library(survcomp)
library(glmnet)
library(GEOquery)
load('RS_train-test_all.rda')
load('pancreatic_GSEA.rda')
load('pancreatic_data.rda')

train<-pancreatic$train
rownames(RS_train)<-rownames(train)

train_rs<-RS_train[,c(1,2,79)]
data<-cbind(train_rs,Pan_cli[rownames(train_rs),])
#colnames(data)<-c('times','status','Risk_score','Sex','Sample_type','Tumor_purity','Age','TMB','FGA','MSI_score','MSI_type')
colnames(data)<-c('times','status','Risk_score','Sex','Sample_type','Tumor_purity','Age','TMB','FGA','MSI_score')
data1<-na.omit(data)
data1$Age<-data1$Age/365
data1$times<-data1$times/30

save(data1,file="train.rda")

dd=datadist(data1)#data駅倬頁data.frame
options(datadist="dd")
fit<-cph(Surv(times,status)~Risk_score+Sex+Sample_type+Tumor_purity+Age+TMB+FGA+MSI_score,data=data1,x=TRUE,y=TRUE,surv=TRUE)
mode(data1)
B<-rcorrcens(Surv(times,status) ~ predict(fit), data =data1)
B
C_Index=(1-B[1])
C_Index

survival<-Survival(fit)
survival1<-function(x)survival(12,x)
survival3<-function(x)survival(36,x)
survival5<-function(x)survival(60,x)#

nom <- nomogram(fit,fun=list(survival1,survival3,survival5),
                funlabel = c("1 Year Survival","3 Year Survival","5 Year Survival"))
save(nom,file = "train_nomplot.rda")


pdf("train_nomplot.pdf",width=12,height=10)
par(family="serif",ps="14")
plot(nom,xfrac=.15)
title(main = "Training cohort")
dev.off()


#！！！！！！！！！！！！！！！！！！！！！！！！！！！！calibration curve！！！！！！！！！！！！！！！！！！！！！！！！！！
fit1<- cph(Surv(times,status)~Risk_score+Sex+Sample_type+Tumor_purity+Age+TMB+FGA+MSI_score,
           data=data1,x=TRUE,y=TRUE,surv=TRUE,time.inc=12)
fit2<- cph(Surv(times,status)~Risk_score+Sex+Sample_type+Tumor_purity+Age+TMB+FGA+MSI_score,
           data=data1,x=TRUE,y=TRUE,surv=TRUE,time.inc=36)
fit3<- cph(Surv(times,status)~Risk_score+Sex+Sample_type+Tumor_purity+Age+TMB+FGA+MSI_score,
           data=data1,x=TRUE,y=TRUE,surv=TRUE,time.inc=60)

cal1 <- calibrate(fit1, cmethod="KM", u=12, m=100, B=1000)#

cal2<- calibrate(fit2, cmethod="KM", u=36, m=100, B=1000)#

cal3<- calibrate(fit3, cmethod="KM", u=60, m=100, B=1000)#




################
pdf("Cal-train-1.pdf",width=7.0,height=5.0)
par(family="serif",ps="10")

plot(cal1,xlim = c(0,1),ylim= c(0,1),errbar.col="orange",col="orange",lwd=3,pch=0,
     xlab="Nomogram-Predicted OS",ylab="Observed OS(%)")
lines(cal1[,c("mean.predicted","KM")],lwd=3,col="orange",type= "b",cex=1,pch=15)#
legend("topleft", inset=0.05,c("1-year"), col=c("orange"), lty =1,lwd =2, bty = "n")

dev.off()



pdf("Cal-train-2.pdf",width=7.0,height=5.0)
par(family="serif",ps="10")
plot(cal2,xlim = c(0,1),ylim= c(0,1),errbar.col="turquoise3",col="turquoise3",lwd=3,pch=1,
     xlab="Nomogram-Predicted OS",ylab="Observed OS(%)")
lines(cal2[,c("mean.predicted","KM")],lwd=3,col="turquoise3",type= "b",cex=1.5,pch=16)#
legend("topleft", inset=0.05,c("3-year"), col=c("turquoise3"), lty =1,lwd =2, bty = "n")

dev.off()



pdf("Cal-train-3.pdf",width=7.0,height=5.0)
par(family="serif",ps="10")
plot(cal3,xlim = c(0,1),ylim= c(0,1),errbar.col="seagreen3",col="seagreen3",lwd=3,#
     xlab="Nomogram-Predicted OS",ylab="Observed OS(%)")  #
lines(cal3[,c("mean.predicted","KM")],lwd=3,col="seagreen3",type= "b",cex=1.5,pch=17)#

abline(0,1,lty=3,lwd=2,col=c(rgb(0,118,192,maxColorValue=255)))

legend("topleft", inset=0.05,c("5-year"), col=c("seagreen3"), lty =1,lwd =2, bty = "n")

dev.off()





################
pdf("Cal-train.pdf",width=6.0,height=5.0)
par(family="serif",ps="10")

plot(cal1,xlim = c(0,1),ylim= c(0,1),errbar.col="orange",col="orange",lwd=3,pch=0,
     xlab="Nomogram-Predicted OS",ylab="Observed OS(%)")
lines(cal1[,c("mean.predicted","KM")],lwd=3,col="orange",type= "b",cex=1,pch=15)#

par(new=TRUE,family="serif",ps="10")
plot(cal2,xlim = c(0,1),ylim= c(0,1),errbar.col="turquoise3",col="turquoise3",lwd=3,pch=1,
     xlab="Nomogram-Predicted OS",ylab="Observed OS(%)")
lines(cal2[,c("mean.predicted","KM")],lwd=3,col="turquoise3",type= "b",cex=1.5,pch=16)#

par(new=TRUE,family="serif",ps="10")
plot(cal3,xlim = c(0,1),ylim= c(0,1),errbar.col="seagreen3",col="seagreen3",lwd=3,#
     xlab="Nomogram-Predicted OS",ylab="Observed OS(%)")  #
lines(cal3[,c("mean.predicted","KM")],lwd=3,col="seagreen3",type= "b",cex=1.5,pch=17)#

abline(0,1,lty=3,lwd=2,col=c(rgb(0,118,192,maxColorValue=255)))

legend("topleft", inset=0.05,c("1-year", "3-year","5-year"), col=c("orange", "turquoise3","seagreen3"), lty =1,lwd =2, bty = "n")

dev.off()
##################################################





