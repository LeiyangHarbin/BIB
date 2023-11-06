
######################################
rm(list=ls())
library(GSVA)
library(clusterProfiler)
library(enrichplot)
library(tidyverse)
library(ggrepel)
library(cowplot)
load("DEGs_result.rda")

immport<-read.csv("Geneappend3.csv")
immport<-immport[,c(6,1)]

ge = DEGs_result$logFC
names(ge) = rownames(DEGs_result)
ge = sort(ge,decreasing = T)
head(ge)

GSEA_immport<-GSEA(ge, TERM2GENE = immport,
                   nPerm=1000,
                   verbose=FALSE,by="fgsea",
                   pAdjustMethod="BH",
                   pvalueCutoff=1)

save(GSEA_immport,file = "GSEA_immport.rda")

################################



######################################
rm(list=ls())
library(GSVA)
library(clusterProfiler)
library(enrichplot)
library(tidyverse)
library(ggrepel)
library(cowplot)
setwd("D:\\OneDrive - hrbmu.edu.cn\\yanglei004\\Pancreatic cancer-metabolic\\DEG\\GSEA\\Immport\\")
load("GSEA_immport.rda")

a<-order(GSEA_immport$NES,decreasing = T)[1:8]

gseaplot<-list()
for (i in 1:length(a)) {
  gseaplot[[i]]<-gseaplot2(GSEA_immport, geneSetID = a[i], title =GSEA_immport$Description[a[i]])
}


immport_GSEA_BRCA1<-plot_grid(gseaplot[[1]],gseaplot[[2]],
                              gseaplot[[3]],gseaplot[[4]],nrow=2)
ggsave("immport_GSEA_1.pdf",width =8,height =8)

save(immport_GSEA_BRCA1,file="immport_GSEA_1.rda")

###########################
b<-order(GSEA_immport$NES,decreasing =F)[1:8]

gseaplot<-list()
for (i in 1:length(b)) {
  gseaplot[[i]]<-gseaplot2(GSEA_immport, geneSetID = b[i], title =GSEA_immport$Description[b[i]])
}


immport_GSEA_BRCA2<-plot_grid(gseaplot[[1]],gseaplot[[2]],
                              gseaplot[[3]],gseaplot[[4]],nrow=2)
ggsave("immport_GSEA_2.pdf",width =8,height =8)

save(immport_GSEA_BRCA2,file="immport_GSEA_2.rda")


immport_GSEA<-plot_grid(immport_GSEA_BRCA1,
                        immport_GSEA_BRCA2,nrow=2)
ggsave("immport_GSEA.pdf",width =8,height=16)

##############################################
B<-GSEA_immport@result
kk<-GSEA_immport
sortkk<-kk[order(kk$NES,decreasing=T)]

p1<-gseaplot2(kk,row.names(sortkk)[c(1,2,3,4)],color="green",pvalue_table =F,rel_heights=c(1.5,0.5,0),subplots = 1:2)
p1
Immport_GSEA_3<-p1
ggsave("Immport_GSEA_3.pdf",width=8.0,height=6.0,units="in")
save(Immport_GSEA_3,file="Immport_GSEA_3.rda")

##################################################

xx<-GSEA_immport
res<-xx@result
res<-res[order(res[,2],decreasing=F),]
res$Significant<-ifelse(res$pvalue<0.05,"S","NS")
res$Enriched<-ifelse(res$NES>0,"E","D")
res<-tidyr::unite(res, "Sig_Enriched",Significant,Enriched)

p<-ggplot(res,aes(x=NES, y=-log10(pvalue)))+scale_x_continuous(limits=c(-5,5))+
  geom_point(aes(color=Sig_Enriched),size=5.0)+scale_color_manual(values=c("green","darkmagenta","red","orange"))+
  theme_bw(base_size=12)+theme(legend.position="None")+
  geom_text_repel(data=subset(res,NES>-10&pvalue<1),family="serif",
                  aes(label=Description),size=2.5, box.padding=unit(0.25, "lines"),
                  point.padding=unit(0.30,"lines"),
                  max.overlaps = getOption("ggrepel.max.overlaps", default =30),
  )#only p<0.05
p
p1<-p+geom_vline(xintercept=0)+geom_hline(yintercept=1.30103)
p2<-p1+theme(axis.title=element_text(size=15,family="serif",color="black",face= "bold"),
             axis.text.y=element_text(size=14,family="serif",color="black",face="bold"),
             axis.text.x=element_text(size=14,family="serif",color="black",face="bold"),
             legend.title=element_text(size=14,family="serif",color="black",face="bold"),
             legend.text=element_text(size=14,family="serif",color="black",face="bold")
)
#p3<-p2+ggtitle("TCGA-risk")+theme(plot.title=element_text(size=20,family="serif",color="black",face="bold",hjust=0.5))
p4A<-p2+labs(x="NES",y="-log10(P-value)")
p4A

GSEA_immport_vol<-p4A

ggsave("GSEA_immport_vol.pdf",width=6,height=6,units="in")
save(GSEA_immport_vol,file="GSEA_immport_vol.rda")

#########################################