library(Seurat)
library(reshape2)
options(stringsAsFactors=F)

#### calculate the common upregulated and downregulated genes
### upregulated genes
files <- list.files(path = "./",pattern = "gene_sum_n_orth.csv",recursive = T)
files
## species except human
for (i in c(1,3:7)) {
  species <- unlist(strsplit(files[i],'/'))[1]
  sp <- read.csv(files[i],row.names = 1)
  sp <- sp[c(2,4)]
  sp <- sp[sp$percentage>=0.2,]
  sp <- sp[!duplicated(sp$Human_id),]
  rownames(sp) <- sp$Human_id
  colnames(sp) <- paste(species,colnames(sp),sep = "_")
  if(i==1){
    diff_gene <- sp
  }else{
    diff_gene <- merge(sp,diff_gene,by='row.names',all=T)
    diff_gene <- data.frame(diff_gene[,-1],row.names = diff_gene$Row.names)
  }
}

## human
i=2
species <- unlist(strsplit(files[i],'/'))[1]
sp <- read.csv(files[i],row.names = 1)
sp <- sp[c(1,3)]
sp <- sp[sp$percentage>=0.2,]
colnames(sp) <- paste(species,colnames(sp),sep = "_")
diff_gene <- merge(sp,diff_gene,by='row.names',all=T)
diff_gene <- data.frame(diff_gene[,-1],row.names = diff_gene$Row.names)
diff_gene[is.na(diff_gene)] <- 0
diff_gene_use <- diff_gene[,grep(colnames(diff_gene),pattern = "*percentage")]
diff_gene_use$num <- rowSums(diff_gene_use>0)
diff_gene_use$gene <- rownames(diff_gene_use)
diff_gene_use3 <- diff_gene_use[diff_gene_use$num>=3,]


###  upregulated genes
files <- list.files(path = "./",pattern = "gene_sum_p_orth.csv",recursive = T)
## species except human
for (i in c(1,3:7)) {
  species <- unlist(strsplit(files[i],'/'))[1]
  sp <- read.csv(files[i],row.names = 1)
  sp <- sp[c(2,4)]
  sp <- sp[sp$percentage>=0.2,]
  sp <- sp[!duplicated(sp$Human_id),]
  rownames(sp) <- sp$Human_id
  colnames(sp) <- paste(species,colnames(sp),sep = "_")
  if(i==1){
    stem_gene <- sp
  }else{
    stem_gene <- merge(sp,stem_gene,by='row.names',all=T)
    stem_gene <- data.frame(stem_gene[,-1],row.names = stem_gene$Row.names)
  }
}

## human
i=2
species <- unlist(strsplit(files[i],'/'))[1]
sp <- read.csv(files[i],row.names = 1)
sp <- sp[c(1,3)]
sp <- sp[sp$percentage>=0.2,]
colnames(sp) <- paste(species,colnames(sp),sep = "_")
stem_gene <- merge(sp,stem_gene,by='row.names',all=T)
stem_gene <- data.frame(stem_gene[,-1],row.names = stem_gene$Row.names)
stem_gene[is.na(stem_gene)] <- 0
stem_gene_use <- stem_gene[,grep(colnames(stem_gene),pattern = "*percentage")]
stem_gene_use$num <- rowSums(stem_gene_use>0)
stem_gene_use$gene <- rownames(stem_gene_use)
stem_gene_use3 <- stem_gene_use[stem_gene_use$num>=3,]



### Exclude genes that are both up- and down-regulated.
genes_exclude <- intersect(rownames(diff_gene_use3),rownames(stem_gene_use3))
diff_gene_use3_use <- diff_gene_use3[!diff_gene_use3$gene%in%genes_exclude,]
stem_gene_use3_use <- stem_gene_use3[!stem_gene_use3$gene%in%genes_exclude,]

### upregulated tfs
human_TF <- read.delim("./TF/Homo_sapiens_TF.txt")
diff_gene_use3_use_tf <- diff_gene_use3[diff_gene_use3$gene%in%human_TF$Symbol,]


### save the files.
write.csv(diff_gene_use3_use_tf,file = "diff_species_nointersected_morethan2_tf.csv",quote = F)
write.csv(diff_gene_use3_use,file = "diff_species_nointersected_morethan2.csv",quote = F)
write.csv(stem_gene_use3_use,file = "stem_species_nointersected_morethan2.csv",quote = F)
write.csv(diff_gene,file = "diff_species_all.csv",quote = F)
write.csv(stem_gene,file = "stem_species_all.csv",quote = F)









