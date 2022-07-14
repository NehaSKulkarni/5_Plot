options(repos = list(CRAN="http://cran.rstudio.com/"))
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")}
BiocManager::install(c("biomaRt","GEOquery","ComplexHeatmap","gridExtra"),update = F)
list.of.packages <- c(c("Hmisc","Tmisc","tidyverse"))
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){install.packages(new.packages)}
library(gridExtra)
library(biomaRt)
library(Hmisc)
library(Tmisc)
library(GEOquery)
library(tidyverse)
library(ComplexHeatmap)


#markers <- c("POU5F1","NANOG","BRCA2")
args <- commandArgs(trailingOnly = T)
markers <- unlist(strsplit(args,","))
print(markers)

compute_var_genes=function(dataset,my_mean=1,cv=0.5,png=F,name=NULL){
  dataset=dataset[apply(dataset, MARGIN = 1, function(x) any(x > 0)), ]

  means <- rowMeans(dataset)
  vars <- apply(dataset,1,var)
  cv2 <- vars/means^2
  Data=data.frame(means, log(cv2))
  colnames(Data)=c("means","cv")

  Data_ordered=Data[with(Data, order(means)), ]
  nls_fit <- nls(cv ~ a+b*(means^c), Data, start = list(a = 1,b=1,c=1))
  i <- order(Data$means)

  Data_ordered_filtered=Data_ordered[which(Data_ordered$means>my_mean & Data_ordered$cv>predict(nls_fit)[i]+cv),]
  list2=as.vector(rownames(Data_ordered_filtered))
  if(png==T) {
    png(paste0(name,"_cv",cv,"_mean_",my_mean,"_Var_genes.png",sep=""),width=4000,height=4000,res=250,pointsize=12)
    par(mar = c(6, 6, 4, 4) + 0.1,mgp=c(2,0.65,0),cex=0.8)
    plot(Data_ordered,cex=0.8,pch=21,bg="grey",col=ifelse((Data_ordered$means>my_mean)&(predict(nls_fit)[i]+cv<Data_ordered$cv),"red","blue"))
    lines(Data$means[i], predict(nls_fit)[i], col = "red")
    dev.off()
  }
  return(dataset[list2,])

}
make_pca<-function(counts,var_genes,type="PCA"){
  pca_df <-prcomp(t(counts[var_genes,]))
  pca_df <- as.data.frame(pca_df$x[,1:2])
  print(identical(rownames(pca_df),rownames(t(counts))))
  pca_df <- cbind(pca_df,as.data.frame(t(counts)))

  return(pca_df)

}
plot_markers <- function(df,markers,return_grid=T,type=c("tsne","pca","umap"),ncol = 2){
  if(type == "pca"){
    colnames(df)[1:2] <- c("PC1","PC2")
  }
  else if(type == "tsne"){
    colnames(df)[1:2] <- c("tSNE1","tSNE2")
  }

  markers <- markers[markers %in% colnames(df)]

  if(type == "pca"){
    plot_data_column = function (data, column){
      ggplot(df)+geom_point(aes_string("PC1","PC2",colour= column))
    }
  }
  else if(type == "tsne"){
    plot_data_column = function (data, column){
      ggplot(df)+geom_point(aes_string("tSNE1","tSNE2",colour= column))
    }
  }

  else if(type == "umap"){
    plot_data_column = function (data, column){
      ggplot(df)+geom_point(aes_string("UMAP1","UMAP2",colour= column))
    }
  }

  myplots <- lapply(markers, plot_data_column, data = pca_df)

  # if T return a grid instead of list
  if(return_grid){
    n <- length(myplots)
    nCol <- floor(sqrt(n))
    grid<-do.call("grid.arrange", c(myplots, ncol=ncol))
    return(grid)
  }

  #return list
  return(myplots)
}

#load in counts
capacitation_counts <- read.table("GSE123055_counts.txt",header = T)
colnames(capacitation_counts) <- gsub(".txt","",colnames(capacitation_counts))

#load metadata and clean up
gds <- getGEO("GSE123055")
capacitation_annotations<- gds[[1]]@phenoData@data
capacitation_annotations <- capacitation_annotations[,c(1,19,42,43,44)]
colnames(capacitation_annotations) <- c("title","name","line","cell_type","transition_day")
capacitation_annotations$name <- gsub(capacitation_annotations$name,pattern = "H9-ctrl",replacement = "H9.ctrl")
capacitation_annotations$title <- gsub(gsub(capacitation_annotations$title,pattern = " ",replacement = "_"),pattern = ":_",replacement = "_")
capacitation_annotations$transition_state<- sapply(capacitation_annotations$cell_type,function(x){strsplit(x,split = " hPSC ")[[1]][2]})
rownames(capacitation_annotations) <- capacitation_annotations$title
capacitation_annotations$line <- sapply(capacitation_annotations$line,function(x){if(grepl("HNES1",x)){"HNES1"} else if(grepl("cR",x) & !grepl("parental",x)){"cR-H9-EOS"} else{"H9-EOS"}})
capacitation_annotations$transition_day <- factor(sapply(capacitation_annotations$transition_day,
                                                         function(x){if(x %in% c("22","23","24","25","26","27","28","22+")){"22+"} else if (x == "N/A"){"control"} else{x}}),levels = c("0","1","2","3","7","10","22+","control"))

# get gene names and mean transcript lengths
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl",host = "useast.ensembl.org"))
genes <- capacitation_counts$ID
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol","gene_biotype","transcript_length"),values=genes,mart= mart)
lengths <- aggregate(G_list[,4], list(G_list$ensembl_gene_id), mean)
colnames(lengths) <- c("ensembl_gene_id","length")
G_list <- G_list[!duplicated(G_list$ensembl_gene_id),1:3]
G_list <- merge(G_list,lengths,by = "ensembl_gene_id")
capacitation_counts <- merge(G_list,capacitation_counts,by.x = "ensembl_gene_id",by.y = "ID")
capacitation_counts$hgnc_symbol <- apply(capacitation_counts, 1,function(row){ if(row[2] == ""){ row[1]} else{row[2]}})
rownames(capacitation_counts) <- make.unique(capacitation_counts$hgnc_symbol)
capacitation_counts <- capacitation_counts[capacitation_counts$gene_biotype == "protein_coding",]
capacitation_lengths <- capacitation_counts[,4,drop = F]
capacitation_counts <- capacitation_counts[,-c(1:4)]
for(cell in colnames(capacitation_counts)){ 
  colnames(capacitation_counts)[which(colnames(capacitation_counts) == cell)] <- capacitation_annotations[capacitation_annotations$name == cell,"title"]
}

#log2FPKM
capacitation_fpkm <- counts2fpkm(as.matrix(capacitation_counts),capacitation_lengths$length)
capacitation_fpkm <- log2(capacitation_fpkm+1)


#PCA
var_genes <- rownames(compute_var_genes(capacitation_fpkm))
capacitation_pca<- make_pca(capacitation_fpkm,var_genes)
capacitation_pca <- merge(capacitation_pca,capacitation_annotations,by="row.names")
capacitation_pca <- column_to_rownames(capacitation_pca,"Row.names")

pdf("Capacitation_PCA.pdf",width = 10,height = 2*(length(markers)+2))
plot_markers(capacitation_pca,markers = c("line","transition_state","transition_day",markers),type = "pca")
dev.off()

#heatmap
sample_order <- rownames(arrange(capacitation_annotations,transition_day,cell_type))
capacitation_fpkm <- capacitation_fpkm[,sample_order]
col_split <- capacitation_annotations$transition_day
row_split <- capacitation_annotations$line

#heatmap averaged
#average fpkm values between groups(cell type and timepoint)
merged <- t(capacitation_fpkm) %>% as.data.frame() %>% rownames_to_column("sample_name") %>% left_join(rownames_to_column(capacitation_annotations,"sample_name"),by = "sample_name")
meaned <- merged %>% pivot_longer(!c(sample_name,title,name,line,cell_type,transition_day,transition_state),names_to = "gene") %>%
  group_by(transition_day,line,gene) %>% summarise(mean = mean(value)) %>% ungroup() %>% mutate(group = paste(line,transition_day,sep = "_")) %>%
  select(c(gene,mean,group)) %>% pivot_wider(names_from = group,values_from = mean) %>% column_to_rownames("gene")
timepoint <- sapply(colnames(meaned), function(x){strsplit(x,"_")[[1]][2]})
timepoint <- sapply(timepoint,function(x){if(x != "control"){paste0("d",x)} else{"H9"}})
timepoint <- factor(timepoint,levels = unique(timepoint))
cell_type <- sapply(colnames(meaned), function(x){strsplit(x,"_")[[1]][1]})


# Heatmap D0 - D10 capacitation
pdf("Capacitation_Heatmap_D0-D10.pdf",height = 1.5*(length(markers)/2),width = 4.5)
Heatmap(as.matrix(meaned[markers,-c(13:14)]),cluster_columns  = F,cluster_rows = F,cluster_column_slices = F,row_names_side = "left",
        column_split = timepoint[-c(13:14)],column_labels = rep("",ncol(meaned[markers,-c(13:14)])),heatmap_legend_param = list(title = "log2FPKM",title_gp = gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)))
dev.off()


pdf("Capacitation_Bar_Plot.pdf",height = 4*(length(markers)/2),width = 10)
bar_markers <- meaned[markers,] 
bar_markers <- bar_markers[,!grepl("22+",colnames(bar_markers))]
bar_markers <- bar_markers %>% rownames_to_column("Gene") %>% pivot_longer(-Gene,names_to = "line",values_to = "expression") %>% separate(line,into = c("line","day"), sep ="_") %>% mutate(day = replace(day,day == "control","H9")) %>% mutate(day = factor(day,levels = unique(day)))
bar_markers$day <- sapply(bar_markers$day,function(x){if(grepl("H9",x)){"H9"} else {paste0("d",x)}})
bar_markers$day <- factor(bar_markers$day,levels = unique(bar_markers$day))
bar_markers$line <- factor(bar_markers$line,levels = c("cR-H9-EOS","HNES1","H9-EOS"))
ggplot(bar_markers) + geom_col(aes(x = day, y = expression, fill = line),position = "dodge") + xlab(element_blank()) + ylab("log2(FPKM)") + theme(legend.title = element_blank()) + facet_wrap(~Gene,ncol = 2,scales = "free",)
dev.off()

#load in petropolous and filter for epi and ETS genes

#load in data
lanner_counts <- read.table("Counts_Lanner_Yan_Blackley.txt",header = T)
lanner_annotations <- read.csv("Sample_Info.csv",row.names = 1) 
colnames(lanner_annotations) <- capitalize(tolower(colnames(lanner_annotations)))
lanner_annotations <- lanner_annotations[lanner_annotations$Dataset == "LANNER",]
lanner_annotations$Lineage<- sapply(lanner_annotations$Lineage,function(x){if(x == "epiblast"){"EPI"} else if(x == "primitive_endoderm"){"PrE"} else if(x == "trophectoderm"){"TE"} else {x}})
rownames(lanner_annotations) <- lanner_annotations$Sample_name
lanner_annotations <- lanner_annotations[,-1]
lanner_annotations$Timepoint<- sapply(rownames(lanner_annotations),function(x){strsplit(x,"_")[[1]][1]})

# get gene names and mean transcript lengths
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl",host = "useast.ensembl.org"))
genes <- lanner_counts$EnsemblID
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol","transcript_length"),values=genes,mart= mart)
lengths <- aggregate(G_list[,3], list(G_list$ensembl_gene_id), mean)
colnames(lengths) <- c("ensembl_gene_id","length")
G_list <- G_list[!duplicated(G_list$ensembl_gene_id),1:2]
G_list <- merge(G_list,lengths,by = "ensembl_gene_id")
lanner_counts <- merge(G_list,lanner_counts,by.x = "ensembl_gene_id",by.y = "EnsemblID")
lanner_counts$hgnc_symbol <- apply(lanner_counts, 1,function(row){ if(row[2] == ""){ row[1]} else{row[2]}})
rownames(lanner_counts) <- make.unique(lanner_counts$hgnc_symbol)
lanner_lengths <- lanner_counts[,3,drop = F]
lanner_counts <- lanner_counts[,colnames(lanner_counts) %in% rownames(lanner_annotations)]

#log2FPKM
lanner_fpkm <- counts2fpkm(as.matrix(lanner_counts),lanner_lengths$length)
lanner_fpkm <- log2(lanner_fpkm+1)

EPI_lanner <- rownames(lanner_annotations[lanner_annotations$Lineage == "EPI",])
lanner_EPI_ETV<- as.data.frame(lanner_fpkm[markers,EPI_lanner])
lanner_ets_means <- lanner_EPI_ETV %>% rownames_to_column("gene") %>% 
  pivot_longer(-gene,names_to = "cells",values_to = "FPKM") %>%
  group_by(gene) %>% summarise(mean = mean(FPKM)) %>% column_to_rownames("gene")
lanner_ets_means <- lanner_ets_means[markers,,drop=F]


ha_bar = rowAnnotation(log2FPKM = anno_barplot(lanner_ets_means, bar_width = 1,ylim = range(8)))
ha_row = rowAnnotation(log2FPKM =anno_boxplot(as.matrix(lanner_EPI_ETV)))

#D1-D3 and Lanner EPI
marker_counts <- meaned[markers,1:8]
bulk_and_lanner <- cbind(lanner_ets_means,marker_counts)


pdf("Petropolous_and_capacitation_D1-D3.pdf",height = 1*(length(markers)/2),width = 4.5)
Heatmap(as.matrix(bulk_and_lanner),cluster_columns  = F,cluster_rows = F,cluster_column_slices = F,row_names_side = "left",
        column_split = unlist(list(factor("EPI"),timepoint[1:8])),column_gap = unit(c(5,1,1,1),"mm"),column_labels = rep("",ncol(bulk_and_lanner)),heatmap_legend_param = list(title = "log2FPKM",title_gp = gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)))
dev.off()

print("DONE")


                                            