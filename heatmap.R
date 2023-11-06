
#################################################
rm(list=ls())
library(ggplot2)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(limma)
library(RColorBrewer)
load("machine_learning_result_all.rda")

result$Cindex<-round(result$Cindex,4)
result$Model<-str_replace(result$Model,"Enet\\[¦Á=0\\]","Ridge")
result$Model<-str_replace(result$Model,"Enet\\[¦Á=1\\]","Lasso")

a<-strsplit2(result$Model,split = " \\+ ") %>% data.frame()
index<-which(a[,1]==a[,2])
index<-c(39,40,99,100,139,140,179,180)
result<-result[-index,]


heatmap<-spread(result,ID,Cindex) %>% tibble::column_to_rownames(var = "Model")
heatmap<-heatmap[,c(2,1)]
colnames(heatmap)<-c("Train","Test")

col_fun <- colorRamp2(
  c(0.5, 0.75, 1), 
  c("#58aaa1", "#FFFFFF", "#e7bc6a")
)
col<-c("#FB8072","#B3DE69")

col_an<-columnAnnotation(Cohort = factor(colnames(heatmap)),
                         col = list(
                           Cohort = c("Train" = col[2],"Test"= col[1]
                                      )),
                         annotation_legend_param = list(
                           title_gp = gpar(fontsize = 13,fontfamily = "serif",fontface = "bold"),
                           labels_gp = gpar(fontsize = 10,fontfamily = "serif"),
                           grid_height = unit(0.6, "cm"), grid_width = unit(0.6, "cm")),
                         border = F,
                         annotation_name_side = "right",
                         annotation_name_gp = gpar(fontsize = 13,fontfamily = "serif",fontface = "bold")
)

## å¹³å‡C-index
mean <- apply(heatmap, 1, mean) 
mean_sort <- sort(mean, decreasing = T)
heatmap <- heatmap[names(mean_sort), ] 
## plot
### barplot
row_bar = rowAnnotation(bar = anno_barplot(mean_sort, bar_width = 0.8, border = FALSE,
                                           gp = gpar(fill = "#aa9b81", col = NA),
                                           add_numbers = F, # ä¸æ˜¾ç¤ºæ•°å­?
                                           width = unit(2, "cm")),
                        show_annotation_name = F)

pdf("machine_complexheatmap-1.pdf",width =6,height = 20)
Heatmap(as.matrix(heatmap),
        col = col_fun,
        show_column_names = F,
        top_annotation = col_an,
        right_annotation = row_bar,
        row_title = NULL,
        column_title = NULL,
        column_split = factor(colnames(heatmap), levels = colnames(heatmap)),
        row_split = factor(rownames(heatmap), levels = rownames(heatmap)),
        column_gap = unit(2, "mm"),
        column_title_side = "top",
        column_title_gp = gpar(fontsize = 17,fontfamily = "serif",fontface = "bold"),
        row_title_gp = gpar(fontsize = 17,fontfamily = "serif",fontface = "bold"),
        cluster_columns = F,
        cluster_rows = F,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 8,fontfamily = "serif",fontface = "bold"),
        heatmap_legend_param = list(title = "Cindex",
                                    title_gp = gpar(fontsize = 13,fontfamily = "serif",fontface = "bold"),
                                    labels_gp = gpar(fontsize = 10,fontfamily = "serif",fontface = "bold"),
                                    legend_height = unit(8, "cm"),
                                    grid_width = unit(0.8, "cm")),
        border = F,
        cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
          grid.text(label = format(heatmap[i, j], digits = 3, nsmall = 4),
                    x, y, gp = gpar(fontsize = 10))
        }
)
dev.off()



annotation_col <- data.frame(
  Cohort = colnames(heatmap)
)
ann_colors = list(
  Cohort = c("train" = col[1],"test"= col[2]))

pdf("machine_heatmap.pdf",width = 8,height = 20)
pheatmap(heatmap,
         col = col_fun,
         cluster_cols = F,
         show_colnames = F,
         display_numbers = T,
         number_format = "%.4f",
         fontsize_row = 11,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         gaps_col = c(1,2),
         heatmap_legend_param = list(title = "Cindex",
                                     title_gp = gpar(fontsize = 13),
                                     labels_gp = gpar(fontsize = 10),
                                     legend_height = unit(6, "cm"),
                                     grid_width = unit(0.5, "cm")))
dev.off()

##################################################