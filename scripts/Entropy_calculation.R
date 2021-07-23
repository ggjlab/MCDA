
############################################################################## calculate entropy using four method
##################### SLICE

#### input
cpm <- apply(datause,2, function(x) (x/sum(x))*10000) # datause: raw count datasets with row as gene, column as cell.
cpm<-as.data.frame(cpm)
library(SLICE)
sc <- construct(exprmatrix=cpm, 
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


##################### StemID

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
load("embryo125_pbmc.RData") ## load the seurat object
zz_all <- read.csv("./colsum.csv",stringsAsFactors=F) #for large datasets, calculate the colsums (sum of umi per cell) in python. df.sum(axis=1) first
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
  #gene <- read.csv("/media/ggj/BACKUP21/DifferentiationForce/DMCA_new/merge/scanpy/gene34212.csv",stringsAsFactors = F,row.names = 1)
  x<-pbmc1@assays$RNA@counts
  #x <- x[rownames(gene),]
  #x1 <-as.data.frame(as_matrix(x))
  x <- as_matrix(x)
  #zz_all <- read.csv("/media/ggj/BACKUP21/DifferentiationForce/DMCA_new/merge/entropy/stemid/colsum.csv",stringsAsFactors=F)
  #mm<-median(zz_all$X0)
  #mm
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

dim(entropy_all)
head(entropy_all)
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


#### input
load("embryo125_pbmc.RData") ## load the seurat object
dt<-pbmc@assays$RNA@data
dt <- dt[,rownames(anno)]

#### convert gene id
orth<-read.delim("./human_mouse_geneid.txt.txt",stringsAsFactors = F) # download from ensemble
humangene <- orth[orth$Mouse.gene.name%in%rownames(dt),]
humangene<-humangene[!duplicated(humangene$Gene.name),]
dt <- dt[humangene$Mouse.gene.name,]
rownames(dt) <- humangene$Gene.name
gene<-mapIds(x=org.Hs.eg.db,keys =rownames(dt),
             column='ENTREZID', keytype='SYMBOL', multiVals = "first" )
z<-data.frame(unlist(gene))
z<-na.omit(z)
colnames(z)<-"id"
cpm1<-dt[rownames(z),]
rownames(cpm1)<-z$id
cpm1[1:3,1:3]
max(cpm1)
ccat.v <- CompCCAT(exp = cpm1, ppiA = net13Jun12.m)

anno <- FetchData(pbmc,vars = c("ident","type"))
anno$caat<-ccat.v
write.csv(anno,file = "embryo_caat.csv",quote = F)

library(ggplot2)
library(ggpubr)
anno$ident <- factor(anno$ident)
ggplot(anno,aes(ident,caat,fill=type,col=type))+geom_boxplot(alpha=0.9)+
  stat_compare_means(aes(label = ..p.signif..),method="wilcox.test")+
  scale_fill_manual(values = c("#fda45e","#6596ef"))+scale_color_manual(values = c("#fda45e","#6596ef"))+
  theme(axis.text.x = element_text( size = 20))+
  theme(axis.text.x = element_text( size = 20),
        axis.text.y = element_text( size = 20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20) )+
  theme_bw() +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))
ggsave("CCAT_entropy.pdf")


##################### cytotrace
library(Seurat)
library(CytoTRACE)
load("embryo125_pbmc.RData") ## load the seurat object
Idents(pbmc) <- "orig.ident"
table(Idents(pbmc))

wt <- subset(pbmc,idents = "WT")
wt <- as.data.frame(as.matrix(wt@assays$RNA@counts))
wt1 <- rowSums(wt>0)>=3
wt <- wt[wt1,]

ko <- subset(pbmc,idents = "KO")
ko <- as.data.frame(as.matrix(ko@assays$RNA@counts))
ko1 <- rowSums(ko>0)>=3
ko <- ko[ko1,]

datasets1 <- list(wt,ko)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
result1 <- iCytoTRACE(datasets1,ncores = 18)
anno <- FetchData(pbmc,vars = c("ident","type"))
anno$CytoTRACE<-result1$CytoTRACE
write.csv(anno,file = "embryo_CytoTRACE.csv",quote = F)

library(ggplot2)
library(ggpubr)
anno$ident <- factor(anno$ident)
ggplot(anno,aes(ident,caat,fill=type,col=type))+geom_boxplot(alpha=0.9)+
  stat_compare_means(aes(label = ..p.signif..),method="wilcox.test")+
  scale_fill_manual(values = c("#fda45e","#6596ef"))+scale_color_manual(values = c("#fda45e","#6596ef"))+
  theme(axis.text.x = element_text( size = 20))+
  theme(axis.text.x = element_text( size = 20),
        axis.text.y = element_text( size = 20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20) )+
  theme_bw() +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))
ggsave("CytoTRACE_entropy.pdf")


