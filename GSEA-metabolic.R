######################################
rm(list=ls())
library(GSVA)
library(clusterProfiler)
library(enrichplot)
library(tidyverse)
library(ggrepel)
library(cowplot)
load("DEGs_result.rda")

metabolic<-data.table::fread("metabolic_gene-114.txt",header = F)
metabolic<-metabolic[,c(2,1)]

ge = DEGs_result$logFC
names(ge) = rownames(DEGs_result)
ge = sort(ge,decreasing = T)
head(ge)

GSEA_metabolic_type<-GSEA(ge, TERM2GENE = metabolic,
                   nPerm=1000,
                   verbose=FALSE,by="fgsea",
                   pAdjustMethod="BH",
                   pvalueCutoff=1)

save(GSEA_metabolic_type,file="GSEA_metabolic_type.rda")

################################



######################################
rm(list=ls())
library(GSVA)
library(clusterProfiler)
library(enrichplot)
library(tidyverse)
library(ggrepel)
library(cowplot)
setwd("D:\\OneDrive - hrbmu.edu.cn\\yanglei004\\Pancreatic cancer-metabolic\\DEG\\GSEA\\metabolic\\")
load("GSEA_metabolic_type.rda")

a<-order(GSEA_metabolic_type$NES,decreasing = T)[1:8]

gseaplot<-list()
for (i in 1:length(a)) {
  gseaplot[[i]]<-gseaplot2(GSEA_metabolic_type,geneSetID=a[i],title=GSEA_metabolic_type$Description[a[i]])
}


metabolic_GSEA_BRCA1<-plot_grid(gseaplot[[1]],gseaplot[[2]],
                              gseaplot[[3]],gseaplot[[4]],nrow=1)
ggsave("metabolic_GSEA_1.pdf",width =16,height =4.2)

save(metabolic_GSEA_BRCA1,file="metabolic_GSEA_1.rda")

###########################
b<-order(GSEA_metabolic_type$NES,decreasing =F)[1:8]

gseaplot<-list()
for (i in 1:length(b)) {
  gseaplot[[i]]<-gseaplot2(GSEA_metabolic_type, geneSetID = b[i], title =GSEA_metabolic_type$Description[b[i]])
}


metabolic_GSEA_BRCA2<-plot_grid(gseaplot[[1]],gseaplot[[2]],
                              gseaplot[[3]],gseaplot[[4]],nrow=1)
ggsave("metabolic_GSEA_2.pdf",width =16,height =4.2)

save(metabolic_GSEA_BRCA2,file="metabolic_GSEA_2.rda")


metabolic_GSEA<-plot_grid(metabolic_GSEA_BRCA1,
                          metabolic_GSEA_BRCA2,nrow=2)
ggsave("metabolic_GSEA.pdf",width =16,height=8.5)

##############################################
B<-GSEA_metabolic_type@result
kk<-GSEA_metabolic_type

sortkk<-kk[order(kk$NES,decreasing=T)]


p1<-gseaplot2(kk,row.names(sortkk)[c(1,2,3,4)],color="green",pvalue_table =F,rel_heights=c(1.5,0.5,0),subplots = 1:2)
p1
metabolic_GSEA_3<-p1
ggsave("metabolic_GSEA_3.pdf",width=8.0,height=6.0,units="in")
save(metabolic_GSEA_3,file="metabolic_GSEA_3.rda")

##################################################

xx<-GSEA_metabolic_type
res<-xx@result
res<-res[order(res[,2],decreasing=F),]
res$Significant<-ifelse(res$pvalue<0.05,"S","NS")
res$Enriched<-ifelse(res$NES>0,"E","D")
res<-tidyr::unite(res, "Sig_Enriched",Significant,Enriched)

p<-ggplot(res,aes(x=NES, y=-log10(pvalue)))+scale_x_continuous(limits=c(-5,5))+
  geom_point(aes(color=Sig_Enriched),size=6.0)+scale_color_manual(values=c("green","darkmagenta","orange","red"))+
  theme_bw(base_size=12)+theme(legend.position="None")+
  geom_text_repel(data=subset(res,abs(NES)>3.0&pvalue<0.01),family="serif",
                  aes(label=Description),size=4.5, box.padding=unit(0.45, "lines"),
                  point.padding=unit(0.50,"lines"),
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

GSEA_metabolic_vol<-p4A

save(GSEA_metabolic_vol,file="GSEA_metabolic_vol.rda")
ggsave("GSEA_metabolic_vol.pdf",width=7,height=7,units="in")


#########################################