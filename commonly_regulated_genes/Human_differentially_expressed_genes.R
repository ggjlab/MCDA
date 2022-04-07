library(Seurat)
library(reshape2)
library(dplyr)
library(Matrix)

#################################################### 
### use the python script to calculate the pca vote: Human_cell_type_pairs.ipynb


#################################################### 
### get the cell type pairs
final<-read.csv("final_n5_pc40_sub150.csv",stringsAsFactors = F,row.names = 1)
final$predictedvalue<-rep(0)
for (i in 1:dim(final)[1]) {
  z=final[i,]
  final$predictedvalue[i]<-apply(z[1:length(z)], 1, max)
}
final$predicted<-rep(0)
for (i in 1:dim(final)[1]) {
  z=final[i,]
  final$predicted[i]<-colnames(z)[apply(z, 1, which.max)]
}

anno<-read.csv("adata_down_anno.csv",stringsAsFactors = F,row.names = 1)
tmp <- anno[rownames(final),]
final$orignal<-tmp$celltype_tissue

### remove immune cells
celltype<-unique(union(final$predicted,final$orignal))
remain<-grep(celltype,pattern ="*Neutr" ,value = T ,invert=T)
remain<-grep(remain,pattern ="*T cell" ,value = T ,invert=T)
remain<-grep(remain,pattern ="*Macrophage" ,value = T ,invert=T)
remain<-grep(remain,pattern ="*Conventional" ,value = T ,invert=T)
remain<-grep(remain,pattern ="*B cell" ,value = T ,invert=T)
remain<-grep(remain,pattern ="*B.cell" ,value = T ,invert=T)
remain<-grep(remain,pattern ="*Megakaryocyte" ,value = T ,invert=T)
remain<-grep(remain,pattern ="*CD8" ,value = T ,invert=T)
remain<-grep(remain,pattern ="*Mast" ,value = T ,invert=T)
remain<-grep(remain,pattern ="*Erythroid" ,value = T ,invert=T)
remain<-grep(remain,pattern ="*NK cell" ,value = T ,invert=T)
remain<-grep(remain,pattern ="*Erythroid" ,value = T ,invert=T)
remain<-grep(remain,pattern ="*Plasmocyte" ,value = T ,invert=T)
remain<-grep(remain,pattern ="*macrophage" ,value = T ,invert=T)
remain<-grep(remain,pattern ="*Monocyte" ,value = T ,invert=T)
remain<-grep(remain,pattern ="*Leukocyte" ,value = T ,invert=T)
remain<-grep(remain,pattern ="*Myeloid" ,value = T ,invert=T)
final_noimmune<-final[final$orignal%in%remain,]
final_noimmune<-final_noimmune[final_noimmune$predicted%in%remain,]
final_noimmune<-final_noimmune[,colnames(final_noimmune)%in%union(remain,c("orignal","predictedvalue","predicted"))]
final_noimmune$predicted<-as.factor(final_noimmune$predicted)
final_noimmune$orignal<-as.factor(final_noimmune$orignal)

aa<-data.frame(table(final_noimmune$predicted,final_noimmune$orignal))
num1<-table(final_noimmune$orignal)
aa$sumall<-rep(0)
for (i in 1:dim(aa)[1]) {
  type1<-aa$Var2[i]
  aa$sumall[i]<-num1[type1]
}
aa$ratio_sumall<-aa$Freq/aa$sumall
aa %>% group_by(Var2) %>% top_n(1, ratio_sumall) ->aa_sort
head(aa_sort)


#################################################### 
### prepare the seurat object
setwd("muscle/integrated")
load("Adult-Muscle_UMAP.RData") ## seurat objects of all the human muscle cells in HCL
adult<-pbmc
load("Fetal-Muscle_UMAP.RData")
fetal<-pbmc

pbmc<-merge(adult,fetal)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", 
                      scale.factor = 10000)

### check the names
cl<-aa_sort[,1:2]
cl<-data.frame(cl)
colnames(cl)<-c("level1","level2")
cl$level1<-as.character(cl$level1)
cl$level2<-as.character(cl$level2)
cl$level1<-gsub(cl$level1,pattern = "\\.",replacement = " ")
cl$level2<-gsub(cl$level2,pattern = "\\.",replacement = " ")
cl$level2<-gsub(cl$level2,pattern = "\\(",replacement = " ")
cl$level2<-gsub(cl$level2,pattern = "\\)",replacement = " ")
cl$level2<-gsub(cl$level2,pattern = "\\-",replacement = " ")
pbmc@meta.data$celltype_tissue<-gsub(pbmc@meta.data$celltype_tissue,pattern = "\\(",replacement = " ")
pbmc@meta.data$celltype_tissue<-gsub(pbmc@meta.data$celltype_tissue,pattern = "\\)",replacement = " ")
pbmc$Celltype<-gsub(pbmc$celltype_tissue,pattern = "\\-",replacement = " ")
setdiff(union(cl$level1,cl$level2),unique(pbmc$Celltype))
setdiff(unique(pbmc$Celltype),union(cl$level1,cl$level2))


### find the markers in the single cell level
library(reshape2)
Idents(pbmc)<-pbmc$Celltype
zz<-data.frame( table(pbmc$Celltype))
zz<-data.frame(zz[,-1],row.names = zz$Var1)
colnames(zz)<-"Freq"
anno<-mergeanno
Idents(pbmc)<-pbmc$Celltype

j=1
for (i in 1:dim(cl)[1]) {
  marker1<-FindMarkers(pbmc,cl[i,j],ident.2 =cl[i,j+1],logfc.threshold = 0,min.pct = 0)
  marker1$gene<-rownames(marker1)
  marker1$cluster<-paste(cl[i,j],cl[i,j+1],sep = "__")
  marker1$p_val_adjbh <- p.adjust(marker1$p_val,method = "BH")
  #marker1 <- marker1[marker1$p_val_adjbh<0.1,]
  marker1 <- marker1[abs(marker1$avg_log2FC)>=0.25,]
  xx <- marker1[3:4]
  marker1$minpct <- apply(xx,1,function(t)max(t))
  marker1 <- marker1[marker1$minpct>=0.1,]
  if(i==1){
    marker_all<-marker1
  }
  else{
    marker_all<-rbind(marker_all,marker1)
  }
  print(i)
  
}

dim(marker_all)
tissue <- "Muscle"
marker_all$tissue <- rep(tissue)
outpath <- './marker/'
save(marker_all,file =paste0(outpath,tissue,"_markerall.RData") )



######################################################  get the upregulated and downregulated gene lists.
setwd(outpath)
source("functions.R")

## merge the marker gene list for all tissues.
files <- list.files(pattern = "markerall.RData",path = "./marker")
tissues <- colsplit(files,pattern = "_",names = c("a1","a2"))
tissues <- tissues$a1
for (i in 1:10) {
  load(files[i])
  assign(tissues[i],marker_all)
  if(i==1){
    marker_all_all <- marker_all
  }else{
    marker_all_all <- rbind(marker_all_all,marker_all)
  }
}

head(marker_all_all)

marker_all_all$species <- rep(species)
marker_all_all$gene_cluster <- paste(marker_all_all$gene,marker_all_all$cluster,sep = "--")
marker_all_all <- marker_all_all[!duplicated(marker_all_all$gene_cluster),]
marker_all_all <- marker_all_all[marker_all_all$p_val_adjbh<0.1,]
marker_all_all <- check_pairname(marker_all_all)
df_n_all <- marker_to_df(marker_all_all,"negative")
df_p_all <- marker_to_df(marker_all_all,"positive")


## get cell type pair information
species<-"Mouse"
cl_pair <- data.frame(celltype_pair=colnames(df_n_all))
cl_pair$species <- rep(species)
cl_pair$Subphyla <- rep('Invertebrate')

## upregulated regulated genes
gene_sum<-data.frame(Gene=rownames(df_n_all),num=rowSums(df_n_all<0))
gene_sum<-gene_sum[order(gene_sum$num,decreasing = T),]
summary(gene_sum$num) 
gene_sum$percentage <- gene_sum$num/dim(cl_pair)[1]
gene_sum_n <- gene_sum

## downregulated regulated genes
gene_sum<-data.frame(Gene=rownames(df_p_all),num=rowSums(df_p_all>0))
gene_sum<-gene_sum[order(gene_sum$num,decreasing = T),]
gene_sum$percentage <- gene_sum$num/dim(cl_pair)[1]
gene_sum_p <- gene_sum

## transform genes into human genes
orthfile <- "./homologousgenes/homo/Human_Hydra_one-one.orth"
target_gene <- data.frame(gene_sum_n)
target_gene$Gene <- rownames(target_gene)
gene_sum_n_orth <- Orth_Convert(target_gene)
gene_sum_n_orth <- gene_sum_n_orth[order(gene_sum_n_orth$num,decreasing = T),]

target_gene <- data.frame(gene_sum_p)
target_gene$Gene <- rownames(target_gene)
gene_sum_p_orth <- Orth_Convert(target_gene)
gene_sum_p_orth <- gene_sum_p_orth[order(gene_sum_p_orth$num,decreasing = T),]

table(gene_sum_p_orth$percentage>=0.2)
table(gene_sum_n_orth$percentage>=0.2)

outpath2 <- "./diff_BHadjust/species_seurat/"
write.csv(gene_sum_p_orth,file = paste0(outpath2,species,"_gene_sum_p_orth.csv"),quote = F)
write.csv(gene_sum_n_orth,file = paste0(outpath2,species,"_gene_sum_n_orth.csv"),quote = F)
write.csv(gene_sum_p,file = paste0(outpath2,species,"_gene_sum_p.csv"),quote = F)
write.csv(gene_sum_n,file = paste0(outpath2,species,"_gene_sum_n.csv"),quote = F)

