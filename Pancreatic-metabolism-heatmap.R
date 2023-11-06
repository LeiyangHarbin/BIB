
################################################
rm(list = ls())
library(GEOquery)
library(ComplexHeatmap)
library(ggplot2)
library(survminer)
load("pancreatic cancer.rda")
load("metabolism_NES-1.rda")
load('cluster-2-train.rda')
label<-read.table('pred-2-train.txt')
colnames(label)<-'cluster'
rownames(label)<-rownames(label_random)
label$cluster[which(label$cluster==1)]<-2
label$cluster[which(label$cluster==0)]<-1

NES_train<-metabolism_NES[,rownames(label)]

C1<-as.matrix(table(NES_train[1,]))
C2<-names(C1[which(C1[,1]>1),])
C3<-as.character(NES_train[1,])
C4<-C3%in%C2
dat<-NES_train[,-which(C4==TRUE)]
NES_train<-dat

Cluster11<-rownames(label)[which(label[,1]==1)]
Cluster22<-rownames(label)[which(label[,1]==2)]

Cluster1<-colnames(NES_train)[colnames(NES_train)%in%Cluster11]
Cluster2<-colnames(NES_train)[colnames(NES_train)%in%Cluster22]

NES_Cluster1<-NES_train[,Cluster1]
NES_Cluster2<-NES_train[,Cluster2]

Sigmark<-function(plot1){
  plot1[which(plot1[,1]>0.05|plot1[,1]==1),2]="ns"
  plot1[which(0.01<plot1[,1]&plot1[,1]<0.05),2]="*"
  plot1[which(0.001<plot1[,1]&plot1[,1]<0.01),2]="**"
  plot1[which(plot1[,1]<0.001),2]="***"
  return(plot1[,2])
}

random_wilcox<-data.frame()
for(i in 1:nrow(NES_train)){
  random_wilcox[i,1]<-wilcox.test(NES_Cluster1[i,],NES_Cluster2[i,],
                                  exact=FALSE,correct=FALSE)$p.value
  random_wilcox[i,2]<-mean(NES_Cluster1[i,])
  random_wilcox[i,3]<-mean(NES_Cluster2[i,])
  if (random_wilcox[i,2]>random_wilcox[i,3]){
    random_wilcox[i,4]<-"High"
  }else if (random_wilcox[i,2]<=random_wilcox[i,3]){
    random_wilcox[i,4]<-"Low"
  }
}
random_wilcox[,5]<-Sigmark(random_wilcox)
row.names(random_wilcox)<-row.names(NES_train)
colnames(random_wilcox)<-c("P-value","Cluster 1","Cluster 2","label","Signifi")

random_wilcox<-random_wilcox[order(random_wilcox$label),]

random_wilcox<-data.frame(random_wilcox)
random_wilcox1<-subset(random_wilcox,P.value<=1)
pathway1high_name<-subset(rownames(random_wilcox1),random_wilcox1[,4] == "High")
pathway2high_name<-subset(rownames(random_wilcox1),random_wilcox1[,4] == "Low")

dat_cluster<-cbind(NES_train[,Cluster1],NES_train[,Cluster2])
dat_cluster_pathway<-rbind(dat_cluster[pathway1high_name,],dat_cluster[pathway2high_name,])

phe<-Pancreatic_cancer
phe_dat<-phe[colnames(dat_cluster),c("oncotree_code","sample_type","curated_subtype_display","os_status","death_or_last_contact_age")]
phe_dat<-na.omit(phe_dat)

phe_dat[which(phe_dat[,4]=="alive"),6]="Alive"
phe_dat[which(phe_dat[,4]=="dead"),6]="Dead"
phe_dat[,4]<-phe_dat[,6]

A1<-median(phe_dat[,5])
B1<-which(phe_dat[,5]<A1)
B2<-which(phe_dat[,5]>=A1)
phe_dat[B1,5]<-"Young"
phe_dat[B2,5]<-"Old"

phe_dat$subtype<-NA
phe_dat[intersect(rownames(phe_dat),Cluster1),"subtype"]<-"Cluster 1"
phe_dat[intersect(rownames(phe_dat),Cluster2),"subtype"]<-"Cluster 2"

dat_cluster_pathway<-dat_cluster_pathway[,rownames(phe_dat)]


library(RColorBrewer)
library(circlize)
library(Hmisc)
col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("#00AFFF", "white", "#FC0E00")
)

my_col1<-c("#8DD3C7","#FFFFB3")
names(my_col1)<-names(table(phe_dat$curated_subtype_display))



col_an<-columnAnnotation(Subtype = factor(phe_dat$subtype),
                         Age = factor(phe_dat$death_or_last_contact_age),
                         Status = factor(phe_dat$os_status),
                         Curated.subtype = factor(phe_dat$curated_subtype_display),
                         Sample.type = factor(phe_dat$sample_type),
                         Oncotree = factor(phe_dat$oncotree_code),
                         col = list(Subtype = c("Cluster 1" = "#da191c", "Cluster 2" = "#2E9FDF"),
                                    Age = c("Old" = "lightcyan","Young" = "sienna1"),
                                    Status = c("Dead" = "#da191c", "Alive" = "#2E9FDF"),
                                    Curated.subtype = my_col1,
                                    Sample.type =  c("Metastasis" = "lightcyan","Primary" = "sienna1"),
                                    Oncotree = c("PAAD" = "#da191c", "PANET" = "#2E9FDF")),
                         annotation_legend_param = list(
                           title_gp = gpar(fontsize = 13,fontfamily = "serif",fontface = "bold"),
                           labels_gp = gpar(fontsize = 10,fontfamily = "serif"),
                           grid_height = unit(0.6, "cm"), grid_width = unit(0.6, "cm")),
                         border = F,
                         annotation_name_side = "right",
                         annotation_name_gp = gpar(fontsize = 13,fontfamily = "serif",fontface = "bold")
                         
)



ha11=rowAnnotation(Pvalue=row_anno_text(as.matrix(random_wilcox1[rownames(dat_cluster_pathway),5]), rot =0,location=unit(0, "mm")))



A1=length(which(random_wilcox$label=="High"))
A2=length(which(random_wilcox$label=="Low"))
ha_cn1=Heatmap(dat_cluster_pathway,row_names_side="left",cluster_rows = F,cluster_columns =F,
               row_split = c(rep("High",A1),rep("Low",A2)),
               top_annotation=col_an,
               use_raster = F,
               #width = unit(20,"cm"),
               #height = unit(17,"cm"),
               show_column_names=F,show_row_names =T,col=col_fun,
               #row_split = rep(c(1,2),c(length(pathway1high_name),length(pathway2high_name))),
               row_title="Metabolism",#column_title="BRCA",
               row_title_gp=gpar(fontsize=15,fontface = "bold",fontfamily="serif"),
               row_names_gp = gpar(fontsize=8,fontface = "bold",fontfamily="serif"),
               column_title_gp=gpar(fontsize=15,fontface = "bold",fontfamily="serif"),
               heatmap_legend_param=list(title="NES",
                                         title_gp = gpar(fontsize = 11,fontface = "bold",fontfamily="serif"),
                                         labels_gp = gpar(fontsize = 8,fontface = "bold",fontfamily="serif"),
                                         legend_height = unit(6, "cm"),
                                         grid_width = unit(0.6, "cm"))
)


ha_cn1A<-ha_cn1+ha11
draw(ha_cn1A)

pdf("Pancreatic_train_metabolism_heatmap.pdf",width=17,height=15)
draw(ha_cn1A)
dev.off()

######################################################