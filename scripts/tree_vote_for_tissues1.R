
setwd("/media/ggj/NEW/DifferentiationForce/Human/tissue/muscle/integrated")
library(Seurat)
load("/media/ggj/NEW/DifferentiationForce/Human/tissue/muscle/Adult-Muscle_UMAP.RData")
adult<-pbmc
adult.marker<-pbmc.markers
load("/media/ggj/NEW/DifferentiationForce/Human/tissue/muscle/Fetal-Muscle_UMAP.RData")
fetal<-pbmc
fetal.marker<-pbmc.markers
head(fetal@meta.data)
rm(pbmc)
tissue1<-"Fetal-Muscle"
tissue2<-"Adult-Muscle"

DimPlot(adult,label = T,pt.size = 0.8)
DimPlot(adult,label = T,pt.size = 0.8,group.by = "tissue_names")
DimPlot(adult,label = T,pt.size = 0.8,group.by = "anno")
DimPlot(fetal,label = T,pt.size = 0.8,group.by = "tissue_names")
DimPlot(fetal,label = T,pt.size = 0.8,group.by = "anno")

fetal_anno<-FetchData(fetal,c("orig.ident","ident"))
fetal_anno$tissue<-paste(tissue1)
adult_anno<-FetchData(adult,c("orig.ident","ident"))
adult_anno$tissue<-paste(tissue2)
anno_tissue<-rbind(fetal_anno,adult_anno)
anno_tissue$cluster<-paste(anno_tissue$tissue,anno_tissue$ident,sep = "_")
head(anno_tissue)


### prepare the cell types
annoall<-read.csv("/media/ggj/NEW/DifferentiationForce/Human/tissue/mergeClusterAnno.csv",stringsAsFactors = F,row.names = 1)
head(annoall)
library(reshape2)
annoall$cluster<-rownames(annoall)
annoall <- annoall[grep(rownames(annoall),pattern = "*Muscle",value = T),]
annoall[22,"celltype_tissue"] <- "B cell_Plasmocyte(Adult-Muscle1)"
zz<-colsplit(annoall$celltype_tissue,pattern = "\\(",names = c("a1","a2"))
head(zz)
annoall$celltype<-zz$a1
head(annoall)

annoall$tissue<-gsub(zz$a2,pattern = "\\)",replacement = "")
annoall$tissue <- gsub(annoall$tissue,pattern = "1",replacement = "")
table(annoall$tissue)

head(annoall)
write.csv(annoall,file = "./mergeClusterAnno_detail.csv",quote = F)

### one to one
annoall_tissue<-annoall[annoall$tissue%in%union(tissue1,tissue2),]
head(annoall_tissue)
head(anno_tissue)
setdiff(anno_tissue$cluster,annoall_tissue$cluster)

anno_tissue$celltype_tissue<-rep(0)
anno_tissue$celltype<-rep(0)

for (i in 1:dim(annoall_tissue)[1]) {
  a1<-annoall_tissue$cluster[i]
  a2<-annoall_tissue$celltype_tissue[rownames(annoall_tissue)==a1]
  a3<-annoall_tissue$celltype[rownames(annoall_tissue)==a1]
  anno_tissue$celltype[anno_tissue$cluster==a1]<-rep(a3)
  anno_tissue$celltype_tissue[anno_tissue$cluster==a1]<-rep(a2)
  print(i)
}
head(anno_tissue)
write.csv(anno_tissue,file = "muscle.integrated_anno.csv",quote = F)

## hcl cells
cell<-read.csv("Muscle.integrated_cell.csv",stringsAsFactors = F,header = T)
head(cell)
tmp<-colsplit(cell$index,pattern = "\\-",names = c("a1","a2"))
tmp1<-colsplit(cell$index,pattern = "\\.",names = c("a1","a2"))
cell$cell1<-tmp$a1
cell$batch<-tmp1$a1
table(cell$batch)

anno_tissue$cell <- rownames(anno_tissue)
tmp<-colsplit(anno_tissue$cell,pattern = "\\.",names = c("a1","a2"))
anno_tissue$batch<-tmp$a1
table(anno_tissue$batch)



#cell$cell1<-gsub(cell$cell1,pattern = "Adultmuscle_4",replacement = "Adultmuscle_3")
tmp1<-colsplit(cell$cell1,pattern = "\\.",names = c("a1","a2"))
cell$batch1<-tmp1$a1
table(cell$batch1)

cell_lack<-setdiff(cell$cell1,rownames(anno_tissue))
table(cell$lack)
cell_lack_yuan<-cell$index[cell$cell1%in%cell_lack]
write.csv(cell_lack_yuan,file = "cell_lack.csv",quote = F,row.names = F)



celluse<-read.csv("Muscle.integrated_celluse.csv",stringsAsFactors = F,row.names = 1)
tmp<-colsplit(celluse$index,pattern = "\\-",names = c("a1","a2"))
celluse$cell1<-tmp$a1

head(celluse)


anno_tissue_sorted<-anno_tissue[celluse$cell1,]
head(anno_tissue_sorted)
table(is.na(anno_tissue_sorted))
write.csv(anno_tissue_sorted,file = "anno_tissue_sorted.csv",quote = F,row.names = F)


#################################################### turn to the python script to calculate the pca vote: tree_vote_for_tissues2.ipynb


####################################################3
###################################################### check the cell types

final<-read.csv("final_n5_pc40_sub150.csv",stringsAsFactors = F,row.names = 1)
dim(final)
final[1:4,1:4]

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


final[1:4,1:4]

anno<-read.csv("adata_down_anno.csv",stringsAsFactors = F,row.names = 1)
tmp <- anno[rownames(final),]
final$orignal<-tmp$celltype_tissue
table(final$predicted,final$orignal)
aa<-data.frame(table(final$predicted,final$orignal))

t<-table(final$predicted,final$orignal)
write.csv(t,file = "n5_pc40_sub150.csv",quote = F)
pp<-read.csv("n5_pc40_sub150.csv",row.names = 1,stringsAsFactors = F)
head(pp)

pp_n5_pc40_sub150_all<-pp

pp<-pp_n5_pc40_sub150_all

head(final)

## remove immune cells
dim(final)
celltype<-unique(union(final$predicted,final$orignal))
length(celltype)
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


length(remain)
remain
final_noimmune<-final[final$orignal%in%remain,]
unique(final_noimmune$orignal)
final_noimmune<-final_noimmune[final_noimmune$predicted%in%remain,]
unique(final_noimmune$predicted)
final_noimmune<-final_noimmune[,colnames(final_noimmune)%in%union(remain,c("orignal","predictedvalue","predicted"))]
dim(final_noimmune)
dim(final)
final_noimmune$predicted<-as.factor(final_noimmune$predicted)
final_noimmune$orignal<-as.factor(final_noimmune$orignal)

t1<-table(final_noimmune$predicted,final_noimmune$orignal)
final_noimmune04<-final_noimmune[final_noimmune$predictedvalue>0.4,]
t2<-table(final_noimmune04$predicted,final_noimmune04$orignal)

zz<-t2/t1

aa1<-data.frame(table(final_noimmune$predicted,final_noimmune$orignal))
aa2<-data.frame(table(final_noimmune04$predicted,final_noimmune04$orignal))

identical(aa1$Var1,aa2$Var1)
identical(aa2$Var1,aa2$Var1)
aa<-aa1
aa$Freq04<-aa2$Freq
head(aa)
aa$ratio<-aa$Freq04/aa$Freq

num1<-table(final_noimmune04$orignal)
num2<-table(final_noimmune$orignal)

head(aa)
aa$sumall<-rep(0)
aa$sum04<-rep(0)

for (i in 1:dim(aa)[1]) {
  type1<-aa$Var2[i]
  aa$sumall[i]<-num2[type1]
  aa$sum04[i]<-num1[type1]
  
}

tail(aa)
aa$ratio_sumall<-aa$Freq/aa$sumall
aa$ratio_sum04<-aa$Freq04/aa$sum04
aa$ratio_04all<-aa$Freq04/aa$sumall

write.csv(aa,file = "ratio.csv",quote = F)


head(aa)
aa
library(dplyr)
aa %>% group_by(Var2) %>% top_n(1, ratio_sumall) ->aa_sort
head(aa_sort)



fetaldata<-as.data.frame(as.matrix(fetal@assays$RNA@counts))
fetaldata[1:4,1:4]
#colnames(fetaldata)<-gsub(colnames(fetaldata),pattern = "Adultmuscle",replacement = "Fetalmuscle")
adultdata<-as.data.frame(as.matrix(adult@assays$RNA@counts))
mergedata<-merge(fetaldata,adultdata,by="row.names",all = T)
mergedata<-data.frame(mergedata[,-1],row.names = mergedata$Row.names)
mergedata[is.na(mergedata)]<-0
mergedata[1:4,1:4]

library(Seurat)
library(Matrix)
pbmc<-CreateSeuratObject(counts=Matrix(as.matrix(mergedata),sparse=T), min.cells = 3, min.features = 0,project = "muscle")
#pbmc<-merge(fetal,adult)
dim(pbmc) #17998 10783

VlnPlot(pbmc, features= c("nFeature_RNA","nCount_RNA"))
# 23873   307
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", 
                      scale.factor = 10000)
#pbmc@assays$RNA@data<-pbmc@assays$RNA@counts
#pbmc@assays$RNA@data<-pbmc@assays$RNA@counts
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
#mergeanno<-anno_tissue_bei
mergeanno<-anno_tissue
pbmc<-AddMetaData(pbmc,metadata = mergeanno)
head(pbmc@meta.data)
table(is.na(pbmc@meta.data))
### check the names
cl<-aa_sort[,1:2]
head(cl)
cl<-data.frame(cl)
colnames(cl)<-c("level1","level2")

cl$level1<-as.character(cl$level1)
cl$level2<-as.character(cl$level2)
cl$level1<-gsub(cl$level1,pattern = "\\.",replacement = " ")
cl$level2<-gsub(cl$level2,pattern = "\\.",replacement = " ")
cl$level2<-gsub(cl$level2,pattern = "\\(",replacement = " ")
cl$level2<-gsub(cl$level2,pattern = "\\)",replacement = " ")
cl$level2<-gsub(cl$level2,pattern = "\\-",replacement = " ")
head(cl)
cl$level2
#cl[5,1]<-"SOX4+ stem cell  Fetal muscle " 

pbmc@meta.data$celltype_tissue<-gsub(pbmc@meta.data$celltype_tissue,pattern = "\\(",replacement = " ")
pbmc@meta.data$celltype_tissue<-gsub(pbmc@meta.data$celltype_tissue,pattern = "\\)",replacement = " ")
head(pbmc@meta.data)
pbmc$Celltype<-gsub(pbmc$celltype_tissue,pattern = "\\-",replacement = " ")

setdiff(union(cl$level1,cl$level2),unique(pbmc$Celltype))
setdiff(unique(pbmc$Celltype),union(cl$level1,cl$level2))

### find the markers in the single cell level
library(reshape2)
Idents(pbmc)<-pbmc$Celltype
zz<-data.frame( table(pbmc$Celltype))
zz<-data.frame(zz[,-1],row.names = zz$Var1)
colnames(zz)<-"Freq"
head(zz)

anno<-mergeanno
#Idents(pbmc)<-pbmc$Celltype

#rm(list = setdiff(ls(),c("pbmc","sample_data","pseudo_anno","anno","dge","cl")))

i=1
j=1
n<-length(cl[i,])
marker1<-FindMarkers(pbmc,cl[i,j],ident.2 =cl[i,j+1])
marker1$gene<-rownames(marker1)
marker1$cluster<-paste(cl[i,j],cl[i,j+1],sep = "_")
marker1<-marker1[marker1$p_val<0.05,]
if(j==1){
  marker<-marker1
}else{
  marker<-rbind(marker,marker1)
}
print(j)

marker11<-marker[,c("gene","cluster","avg_logFC")]


## make the positive chart
marker11_p<-marker11[marker11$avg_logFC>0,]
df_p<-dcast(marker11_p,gene~cluster,value.var="avg_logFC")
a1<-colnames(df_p)[-1]
df_p<-data.frame(df_p[,-1],row.names = df_p$gene)
colnames(df_p)<-a1
df_p[is.na(df_p)]<-0
ds_p<-data.frame(matrix(ncol =4 ,nrow =dim(df_p)[1]))
colnames(ds_p)<-c("gene","num","sum","celltype")
ds_p$gene<-rownames(df_p)
ds_p$num<-rowSums(df_p!=0)
ds_p$sum<-rep(dim(df_p)[2])
ds_p$celltype<-rep(cl[i,j+1])

## make the negative chart
marker11_n<-marker11[marker11$avg_logFC<0,]
df_n<-dcast(marker11_n,gene~cluster,value.var="avg_logFC")
a1<-colnames(df_n)[-1]
df_n<-data.frame(df_n[,-1],row.names = df_n$gene)
colnames(df_n)<-a1
df_n[is.na(df_n)]<-0
ds_n<-data.frame(matrix(ncol =4 ,nrow =dim(df_n)[1]))
colnames(ds_n)<-c("gene","num","sum","celltype")
ds_n$gene<-rownames(df_n)
ds_n$num<-rowSums(df_n!=0)
ds_n$sum<-rep(dim(df_p)[2])
ds_n$celltype<-rep(cl[i,j+1])
print(i)

marker_all<-marker
df_p_all<-df_p
ds_p_all<-ds_p
df_n_all<-df_n
ds_n_all<-ds_n






for (i in 2:dim(cl)[1]) {
  marker1<-FindMarkers(pbmc,cl[i,j],ident.2 =cl[i,j+1])
  marker1$gene<-rownames(marker1)
  marker1$cluster<-paste(cl[i,j],cl[i,j+1],sep = "_")
  marker1<-marker1[marker1$p_val<0.1,]
  if(j==1){
    marker<-marker1
  }
  else{
    marker<-rbind(marker,marker1)
  }
  marker<-marker[!marker$cluster%in%unique(marker_all$cluster),]
  marker11<-marker[,c("gene","cluster","avg_logFC")]
  
  
  ## make the positive chart
  marker11_p<-marker11[marker11$avg_logFC>0,]
  df_p<-dcast(marker11_p,gene~cluster,value.var="avg_logFC")
  a1<-colnames(df_p)[-1]
  df_p<-data.frame(df_p[,-1],row.names = df_p$gene)
  colnames(df_p)<-a1
  df_p[is.na(df_p)]<-0
  ds_p<-data.frame(matrix(ncol =4 ,nrow =dim(df_p)[1]))
  colnames(ds_p)<-c("gene","num","sum","celltype")
  ds_p$gene<-rownames(df_p)
  ds_p$num<-rowSums(df_p!=0)
  ds_p$sum<-rep(dim(df_p)[2])
  ds_p$celltype<-rep(cl[i,j+1])
  
  ## make the negative chart
  marker11_n<-marker11[marker11$avg_logFC<0,]
  df_n<-dcast(marker11_n,gene~cluster,value.var="avg_logFC")
  a1<-colnames(df_n)[-1]
  df_n<-data.frame(df_n[,-1],row.names = df_n$gene)
  colnames(df_n)<-a1
  df_n[is.na(df_n)]<-0
  ds_n<-data.frame(matrix(ncol =4 ,nrow =dim(df_n)[1]))
  colnames(ds_n)<-c("gene","num","sum","celltype")
  ds_n$gene<-rownames(df_n)
  ds_n$num<-rowSums(df_n!=0)
  ds_n$sum<-rep(dim(df_n)[2])
  ds_n$celltype<-rep(cl[i,j+1])
  
  
  marker_all<-rbind(marker_all,marker)
  df_p_all<-merge(df_p_all,df_p,by="row.names",all = T)
  df_p_all<-data.frame(df_p_all[,-1],row.names = df_p_all$Row.names)
  ds_p_all<-rbind(ds_p_all,ds_p)
  
  df_n_all<-merge(df_n_all,df_n,by="row.names",all = T)
  df_n_all<-data.frame(df_n_all[,-1],row.names = df_n_all$Row.names)
  ds_n_all<-rbind(ds_n_all,ds_n)
  print(i)
  
  
  
}

df_p_all[is.na(df_p_all)]<-0
df_n_all[is.na(df_n_all)]<-0
dim(df_p_all)
dim(df_n_all)
dim(marker_all)
head(ds_n_all)





##### sum the cell types
geneuse<-unique(ds_n_all$gene)


gene_sum<-data.frame(Gene=rownames(df_n_all),num=rowSums(df_n_all<0))
head(gene_sum)
gene_sum<-gene_sum[order(gene_sum$num,decreasing = T),]
summary(gene_sum$num) 

write.csv(gene_sum,file ="gene_sum.csv",quote = F)
target_gene<-gene_sum[gene_sum$num>ceiling(dim(df_n_all)[2]*0.2),]
#target_gene<-gene_sum[gene_sum$num>10,]
dim(target_gene)

rm(pbmc)
rm(list = ls(pattern = "Adultmuscle"))
rm(list = ls(pattern = "Fetalmuscle"))
rm(fetaldata,mergedata,adultdata)
rm(adult,fetal)

save.image("./muscle.RData")


dt<-df_n_all[target_gene$Gene,]
library(RColorBrewer)
brewer.pal.info
display.brewer.all()
brewer.pal(9,"Blues")
brewer.pal(9,"BuGn")
pdf("muscle_allbatches_heatmap.pdf",height = 40,width = 40)
a<-pheatmap::pheatmap(dt,color =rev(brewer.pal(9,"Blues")),cluster_rows = T, fontsize_row = 6)
dev.off()
str(a)


muscle_df_n_all<-df_n_all
muscle_ds_n_all<-ds_n_all
muscle_n_target_gene<-target_gene
muscle_n_gene_sum<-gene_sum


gene_sum<-data.frame(Gene=rownames(df_p_all),num=rowSums(df_p_all>0))
head(gene_sum)
gene_sum<-gene_sum[order(gene_sum$num,decreasing = T),]
write.csv(gene_sum,file = "fetal_gene_sum.csv",quote = F)
target_gene<-gene_sum[gene_sum$num>ceiling(dim(df_p_all)[2]*0.2),]

muscle_df_p_all<-df_p_all
muscle_ds_p_all<-ds_p_all
muscle_p_target_gene<-target_gene
muscle_p_gene_sum<-gene_sum


save(muscle_df_n_all,muscle_ds_n_all,muscle_n_target_gene,muscle_n_gene_sum,
     muscle_df_p_all,muscle_ds_p_all,muscle_p_target_gene,muscle_p_gene_sum,
     file = "muscle_gene.RData")


