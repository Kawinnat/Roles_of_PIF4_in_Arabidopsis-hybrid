## Overall transcriptomic profiles #1B

```R
# load required packages
library(psych)
pairs.panels(input_cor, scale=TRUE,cex=2)
```



## Creating the bar graph #1C

```R
# load required packages
library(ggpubr)
library(ggplot2)
library(tidyverse)
library(hashmap)

# import data
c24col_c24pi44_ZT0 <- read.csv("TPM-table.csv")
gene_anno <- read.csv("gene_name_33622.txt", sep = " ") # annotation of gene names (gene_id/gene_name)

# create long format
c24col_c24pi44_ZT0_lf <- c24col_c24pi44_ZT0 %>% rownames_to_column() %>%
  gather(.,key=line, value="TPM", -rowname) %>% 
  separate( col = line, into = c("line","temp", "time","das"), sep = "_")

c24col_c24pi44_ZT0_lf$line2 <- factor(c24col_c24pi44_ZT0_lf$line, levels = c("c24col","c24pif4"))
sp <- levels(c24col_c24pi44_ZT0_lf$line2) 
B <- hashmap(sp, c("C24xCol","C24xpif4-101"))
c24col_c24pi44_ZT0_lf$Hybrid <- B[[c24col_c24pi44_ZT0_lf$line2,]]

c24col_c24pi44_ZT0_lf$temp2 <- factor(c24col_c24pi44_ZT0_lf$temp, levels = c("22c","27c"))
temp_cond <- levels(c24col_c24pi44_ZT0_lf$temp2) 
C <- hashmap(temp_cond, c("22°C","27°C"))
c24col_c24pi44_ZT0_lf$Temp <- C[[c24col_c24pi44_ZT0_lf$temp2,]]

# create a function for the bar graph 
makeBARplot2 <- function(gene){
  gene_name <- gene_anno %>% filter(rowname==gene) %>% select(Gene.Name) 
  ggplot(c24col_c24pi44_ZT0_lf %>% filter(rowname==gene), aes(x=Hybrid,y=TPM,fill=Hybrid)) + geom_bar(stat = "identity") + facet_grid(das~Temp) + theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust = 1),axis.title = element_blank(),) +ggtitle(paste(gene," ","(" ,gene_name$Gene.Name,")", sep = "" ))+scale_fill_manual(values=c("#38761dff","#FF4500"))}

# plot the graph
p1 <- makeBARplot2("AT4G16780")
p2 <- makeBARplot2("AT4G28720")
p3  <- makeBARplot2("AT4G32280")  

ggarrange(plotlist = list(p1,p2,p3),ncol = 3,nrow = 2,common.legend = T)

# to save as a picture
#png(filename = "PIF4_targets_bar_graph.png", res = 300, width = 2000, height = 1500)
#ggarrange(plotlist = list(p1,p2,p3),ncol = 3,common.legend = T)
#dev.off()
```



## Differential expressed gene analyses

```R
# load required packages
library(DESeq2)
result_all <- list()
result_sig <- list() #padj < 0.05, lgfc >1
result_names_up <-list()
result_names_down <-list()

# An example of DEGs determinaltion of c24col between 27c/22c 
input <- c24col_c24pi44_ZT0 %>% select(grep("c24col", colnames(.)))
coldata_name <- coldata %>% filter(col_name%in% colnames(input))
coldata_name$temp <- factor(coldata_name$temp, levels = c("22c","27c"))
ddsFullCountTable <- DESeqDataSetFromMatrix(countData = input,colData = coldata_name,design = ~ rep  +temp)
dds <- DESeq(ddsFullCountTable)
result_all <-as.data.frame(results(dds)) %>% rownames_to_column() %>% left_join(gene , by = "rowname") %>% column_to_rownames()
result_sig <-result_all %>% rownames_to_column() %>%filter(padj < 0.05 & (log2FoldChange>1 | log2FoldChange<(-1))) 
result_names_up[paste("up",sep = "_")] <- result_all %>% rownames_to_column() %>% filter(padj < 0.05 & (log2FoldChange>1 )) 
result_names_down[paste("down",sep = "_")] <- result_all %>% rownames_to_column() %>% filter(padj < 0.05 & (log2FoldChange<(-1))) 
```



## Creating the heatmap #2A

```R
# load required packages
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(fastcluster)

# create column annotation
col_ann <- HeatmapAnnotation(df = data.frame( Age = c("4 DAS","4 DAS","7 DAS","7 DAS"), Hybrid = c("C24xCol","C24xpif4-101","C24xCol","C24xpif4-101") ), col = list(Age =  c("7 DAS" = "lightblue","4 DAS"="pink"), Hybrid=c(  "C24xCol"="#38761dff" , "C24xpif4-101"="#FF4500")), show_annotation_name = F,annotation_legend_param =list(direction = "vertical",legend_width=unit(4,"cm"),title_position = "topleft") )

# plot heatmap
heatmap <- 	Heatmap(all_clusters , col = colorRamp2(c(-1.5,0,1.5), c("#4B0082", "white", "orange")), name = "log2(27°C/22°C)", cluster_columns = FALSE, show_column_dend = FALSE, cluster_rows = FALSE, show_row_names = FALSE, show_column_names = FALSE, show_heatmap_legend = T,heatmap_legend_param = list(direction = "horizontal",legend_width=unit(4,"cm"),title_position = "topcenter"),border = T,top_annotation = col_ann)+ 
  			Heatmap(all_clusters[,5], col = colorRamp2(c(1,2), c("red","white")), cluster_columns = FALSE, show_column_dend = FALSE, cluster_rows = FALSE, show_row_names = FALSE, show_column_names = FALSE,na_col = "white",width = unit(6, "mm"), show_heatmap_legend = F,border = T)

#png(filename = "heatmap_DEGs_c24pif4_c24col.png",res = 300, width = 2000,height = 3500)
draw(heatmap , row_gap = unit(2.5, "mm"), row_split = all_clusters$cluster,heatmap_legend_side = "bottom",annotation_legend_side = "bottom")
#dev.off()
```



## GO enrichment analyses

```R
library(gprofiler2)
# extract gene name from each cluster
input_go  <- list()
clust_n <- c(1:6) %>% as.character()
for (i in 1:6) {
  input_go[paste("cluster",i,sep = "_")] <- all_clusters %>% filter(cluster==clust_n[i]) %>% select(rowname) 
}

# go term analyses
go_term <- gost(input_go, organism = "athaliana",correction_method = "fdr", custom_bg =input_table%>% rownames_to_column() %>% select(rowname) %>%unlist%>%as.character(),evcodes = T,sources = "GO:BP")

#2B
ggplot(GO_filtered_top5 ,aes(x=temp,y=reorder(term_name,-log10(p_value)),color=-log10(p_value),size=gene_ratio))+ geom_point() + scale_color_continuous(low = "blue", high = "red")+
  theme(plot.title = element_text( size=14),
       axis.title.x = element_blank(),
       axis.title.y =  element_blank(),
       axis.text.x = element_blank(),
       axis.text  = element_text( size=12))  + xlab("Gene ratio")+ ggtitle("Biological processes")+ facet_grid(rows = vars(query), scales = "free", space = "free",switch = "y")+ scale_y_discrete(position = "right")
```



