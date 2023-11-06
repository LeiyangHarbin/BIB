
################################################
rm(list=ls())
library(survminer)
library(survival)
library(survcomp)
library(prodlim)
library(ggplot2)
library(RColorBrewer)
library(MetBrewer)
load('RS_train-test_all.rda')
data<-cbind(RS_train$OS.time,RS_train$OS,RS_train$`SetpCox [backward] + RSF`)
colnames(data)<-c("times","status","risk_score")
data<-as.data.frame(data)
data$risk_score<-as.numeric(as.character(data$risk_score))
data$times<-data$times/30


sur.cut <- surv_cutpoint(data,time="times",event="status",variables="risk_score",minprop=0.4)
summary(sur.cut)
sur.cat <- surv_categorize(sur.cut)
head(sur.cat)
fit1 <- survfit(Surv(times, status) ~risk_score, data = sur.cat)  
Bcox<-coxph(Surv(times, status)~risk_score, data = sur.cat)
summcph<-summary(Bcox)
summcph
summcph$sctest[3]
#mypalette<-brewer.pal(8,"Dark2")[c(5,6)]
#mypalette<-met.brewer("VanGogh2")[c(1,5)]
mypalette<-c("lawngreen","red3")
p1<-ggsurvplot(fit1,legend.title="Subtype",
              main="Overall survival",
              legend.labs=c("High","Low"),
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
p3<-p2$plot+ggtitle("Kaplan-Meier Curve for training cohort")+#设置plot图形
  theme(plot.title=element_text(size=15,family="serif",color="black",face="bold",hjust=0.1))
p4<-p3+annotate("text", x=8, y=0.15, label="P=2.48e-139 ",color="black",cex=5)
p4
p5<-p2$table+theme(legend.position="none")#设置table图形
p6<-ggarrange(p4,p5,heights=c(5,1.8),ncol=1,nrow=2, align = "v")
p6
Pancreatic_train_ov<-p6
p6
save(Pancreatic_train_ov,file="Pancreatic_train_ov.rda")
ggsave("Pancreatic_train_ov.pdf",width=6,height=6)

dev.off()

#################################


