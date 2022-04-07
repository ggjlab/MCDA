
############################################################################## calculate entropy using three methods for KO and WT embryos
##################### SLICE

#### input
load("embryo125_seuratobject.RData") ## load the seurat object (named the seurat object as "pbmc" here)
dataset <- pbmc@assays$RNA@data
dataset <- as_matrix(dataset)
dataset <- as.data.frame(dataset)
library(SLICE)
sc <- construct(exprmatrix=dataset, 
                cellidentity=annouse$orig.ident,
                projname="SLICE")

#### load the kappa similarity matrix
load("./SLICE-master/SLICE.a12242016/data/mm_km.Rda")

#### bootstrap calculation of scEntropy
sc <- getEntropy(sc, km=km,                             # use the pre-computed kappa similarity matrix of mouse genes
                 calculation="bootstrap",               # choose the bootstrap calculation
                 B.num=25,                             # 100 iterations
                 exp.cutoff=1,                          # the threshold for expressed genes
                 B.size=1000,                           # the size of bootstrap sample
                 clustering.k=floor(sqrt(1000/2)),      # the number of functional clusters  
                 random.seed=123)                    # set the random seed to reproduce the results in the paper

slice_entropy <- sc@entropies

#### visualization
annouse <- cbind(anno[rownames(entropy),],entropy)
head(annouse)


library(ggpubr)
library(ggplot2)
p1 <- ggplot(annouse,aes(embryo,scEntropy.bootstrap,fill=orig.ident))+geom_boxplot(alpha=0.7)+
  scale_fill_manual(values = c("orange","dark green"))+
  stat_compare_means(aes(label = ..p.signif..),method="wilcox.test")+
  theme(axis.text.x = element_text( size = 20,angle = 90),
        axis.text.y = element_text( size = 20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20) )
ggsave("SLICE_entropy.pdf")


##################### StemID (modified function for large datasets)

#### load class definition and functions
source("./StemID-master/RaceID2_StemID_class.R")
as_matrix <- function(mat){
  
  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  
  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x
  
  for (i in seq_along(val)){
    tmp[row_pos[i],col_pos[i]] <- val[i]
  }
  
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}

#### input data
load("embryo125_seuratobject.RData") ## load the seurat object (named the seurat object as "pbmc" here)
#zz_all <- read.csv("./colsum.csv",stringsAsFactors=F) #for large datasets, calculate the colsums (sum of umi per cell) in python. df.sum(axis=1) first
zz_all <- Matrix::colSums(pbmc@assays$RNA@counts)

mm<-median(zz_all$x)
tissue <- unique(pbmc@meta.data$batch)
Idents(pbmc)<-"batch"
i=1
pbmc1<-subset(pbmc,idents = tissue[i])
x<-pbmc1@assays$RNA@counts
x <- as_matrix(x)
zz<-apply(x,2,sum,na.rm=TRUE)
normalized_mm <- as.data.frame(t(t(x)/zz)*mm + .1 )
normalized_mm[1:5,1:5]
normalized_mm_sum<-apply(normalized_mm,2,sum)
probs   <- t(t(normalized_mm)/normalized_mm_sum)
nrowuse <- nrow(normalized_mm)
entropy <- -apply(probs*log(probs)/log(nrowuse),2,sum)
entropy<-data.frame(entropy,row.names = colnames(pbmc1))
entropy_all <- entropy

for(i in 2:length(tissue)){
  pbmc1<-subset(pbmc,idents = tissue[i])
  ## input data
  x<-pbmc1@assays$RNA@counts
  x <- as_matrix(x)
  zz<-apply(x,2,sum,na.rm=TRUE)
  normalized_mm <- as.data.frame(t(t(x)/zz)*mm + .1 )
  normalized_mm[1:5,1:5]
  normalized_mm_sum<-apply(normalized_mm,2,sum)
  probs   <- t(t(normalized_mm)/normalized_mm_sum)
  nrowuse <- nrow(normalized_mm)
  entropy <- -apply(probs*log(probs)/log(nrowuse),2,sum)
  entropy<-data.frame(entropy,row.names = colnames(pbmc1))
  entropy_all <- rbind(entropy_all,entropy)
  rm(pbmc1,probs,x,normalized_mm)
  gc()
  print(i)
}

colnames(entropy_all)<-"StemID"

library(Seurat)
anno <- FetchData(pbmc,vars=c("cluster_final","type"))
anno<-cbind(anno[rownames(entropy_all),],entropy_all)

library(ggplot2)
library(ggpubr)
anno$cluster_final <- factor(anno$cluster_final)
ggplot(anno,aes(cluster_final,StemID,fill=type,col=type))+geom_boxplot(alpha=0.9)+
  stat_compare_means(aes(label = ..p.signif..),method="wilcox.test")+
  scale_fill_manual(values = c("#fda45e","#6596ef"))+scale_color_manual(values = c("#fda45e","#6596ef"))+
  theme(axis.text.x = element_text( size = 20))+
  theme(axis.text.x = element_text( size = 20),
        axis.text.y = element_text( size = 20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20) )+
  theme_bw() +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))
ggsave("StemID_entropy.pdf")


##################### CCAT
library(Seurat)
library(SCENT)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(qlcMatrix)
library(ggplot2)
library(ggpubr)

#### input
load("embryo125_seuratobject.RData") ## load the seurat object (named the seurat object as "pbmc" here)
dt<-pbmc@assays$RNA@data
dt <- dt[,rownames(anno)]


#### using homologous gene mapping
dduse<-read.csv("homologous _genes/all_all/Mouse_human_all_all.csv") # gene list download from the ensemble website.
dduse <- dduse[dduse$Human_id!="",]
dduse <- dduse[dduse$Species_id!="",]

### make the adjacency matrix
genelist <- dduse
genelist$num <- rep(1)
dt <- dcast(genelist,Human_id~Species_id,value.var = "num")
dt <- na.omit(dt)
dt <- data.frame(dt[,-1],row.names = dt$Human_id)
dt[is.na(dt)] <- 0
dt <- t(dt)


### nomalized to 1 for each human gene
hydra_human <- apply(dt,2, function(x) (x/sum(x))) 
hydra_human[1:7,1:7]
tmp <- hydra_human[,colSums(hydra_human)<1]
dim(tmp)


### prepare the hydra*cell matrix
cpm <- pbmc@assays$RNA@data
cpm1 <- cpm[intersect(rownames(cpm),rownames(hydra_human)),]
hydra_human1 <- hydra_human[intersect(rownames(cpm),rownames(hydra_human)),]
library(Matrix)
cpm1 <- Matrix(as.matrix(cpm1),sparse=T)
cpm1 <- t(cpm1)

### matrix multiplication to weight the hydra genes to human genes
library(Matrix)
hydra_human1 <- Matrix(as.matrix(hydra_human1),sparse=T)
cell_human <- cpm1 %*% hydra_human1
save(cell_human,file = "cell_human.RData")


### calculate entropy
cell_human1 <- t(cell_human) 
gene<-mapIds(x=org.Hs.eg.db,keys =rownames(cell_human1),
             column='ENTREZID', keytype='SYMBOL', multiVals = "first" )
z<-data.frame(unlist(gene))
z<-na.omit(z)
colnames(z)<-"id"
cell_human1<-cell_human1[rownames(z),]
rownames(cell_human1)<-z$id
ccat.v <- CompCCAT(exp = cell_human1, ppiA = net13Jun12.m);
head(ccat.v)

### visualization
anno <- FetchData(pbmc,vars = c("ident","type"))
anno$caat<-ccat.v
anno$ident <- factor(anno$ident)
ggplot(anno1,aes(ident,caat,fill=type,col=type))+geom_boxplot(alpha=0.9)+ylab("CCAT")+
  stat_compare_means(aes(label = ..p.signif..),method="wilcox.test")+
  scale_fill_manual(values = c("#fda45e","#6596ef"))+scale_color_manual(values = c("#fda45e","#6596ef"))+
  #theme(axis.text.x = element_text( size = 20))+
  theme(axis.text.x = element_text( size = 20),
        axis.text.y = element_text( size = 40),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20) )+
  theme_bw() +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))

ggsave("CCAT_entropy.pdf")

