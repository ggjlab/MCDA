library(Seurat)
library(reshape2)
options(stringsAsFactors=F)
source("functions.R")

############################################ Find Differentially expressed genes between cell type pairs

## load and normalize the data
load("Hydra_dge_anno.RData")
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", 
                      scale.factor = 10000)

## load cell type information
cl<-read.csv("hydra_celltype_use.csv",stringsAsFactors = F)  # information about cell type pairs with the development branch in each row， generated using the Hydra_cell_type_pairs.ipynb
ct<-union(cl$level1,union(cl$level2,union(cl$level3,union(cl$level4,cl$level5))))

## check the cell type names
setdiff(ct,pbmc$cluster.long.portal)
setdiff(pbmc$cluster.long.portal,ct)

## find the markers in the single cell level
library(reshape2)
Idents(pbmc)<-pbmc$cluster.long.portal
for (i in dim(cl)[1]:1) {
  for (j in 1:(length(setdiff(cl[i,],""))-1)) {
    marker1<-FindMarkers(pbmc,cl[i,j],ident.2 =cl[i,j+1],logfc.threshold = 0,min.pct = 0)
    marker1$gene<-rownames(marker1)
    marker1$cluster<-paste(cl[i,j],cl[i,j+1],sep = "__")
    marker1$p_val_adjbh <- p.adjust(marker1$p_val,method = "BH")
    marker1 <- marker1[abs(marker1$avg_log2FC)>=0.25,]
    xx <- marker1[3:4]
    marker1$minpct <- apply(xx,1,function(t)max(t))
    marker1 <- marker1[marker1$minpct>=0.1,]
    if(j==1){
      marker11<-marker1
    }
    else{
      marker11<-rbind(marker11,marker1)
    }
    
  }
  if(i==dim(cl)[1]){
    marker_all<-marker11
  }
  else{
    marker_all<-rbind(marker_all,marker11)
  }
  print(i)
  
}

species <- "Hydra"
marker_all$species <- rep(species)
marker_all$gene_cluster <- paste(marker_all$gene,marker_all$cluster,sep = "--")
marker_all <- marker_all[!duplicated(marker_all$gene_cluster),]
marker_all <- marker_all[marker_all$p_val_adjbh<0.1,]
marker_all_all <- check_pairname(marker_all)
df_n_all <- marker_to_df(marker_all_all,"negative")
df_p_all <- marker_to_df(marker_all_all,"positive")

## save the marker gene list
outpath1 <- "./diff_BHadjust/species_seurat/"
write.csv(marker_all_all,file = paste0(outpath,species,"_marker_all.csv"), quote = F)
write.csv(df_p_all,file = paste0(outpath,species,"_stem_matrix_logFC.csv"),quote = F)
write.csv(df_n_all,file = paste0(outpath,species,"_diff_matrix_logFC.csv"),quote = F)

## get cell type pair information
cl_pair <- data.frame(celltype_pair=colnames(df_n_all))
cl_pair$species <- rep(species)
cl_pair$Subphyla <- rep('Invertebrate')

## upregulated regulated genes
a <- colsplit(rownames(df_n_all),names = c("a1","a2"),pattern = "-")
rownames(df_n_all) <- a$a1
gene_sum<-data.frame(Gene=rownames(df_n_all),num=rowSums(df_n_all<0))
gene_sum<-gene_sum[order(gene_sum$num,decreasing = T),]
summary(gene_sum$num) 
gene_sum$percentage <- gene_sum$num/dim(cl_pair)[1]
gene_sum_n <- gene_sum

## downregulated regulated genes
a <- colsplit(rownames(df_n_all),na
a <- colsplit(rownames(df_p_all),names = c("a1","a2"),pattern = "-")
rownames(df_p_all) <- a$a1
gene_sum<-data.frame(Gene=rownames(df_p_all),num=rowSums(df_p_all>0))
gene_sum<-gene_sum[order(gene_sum$num,decreasing = T),]
gene_sum$percentage <- gene_sum$num/dim(cl_pair)[1]
gene_sum_p <- gene_sum

## transform genes into human genes
orthfile <- "homologous _genes/one_one/Human_Hydra_one-one.orth"
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


#### Do the same process with other species, 

