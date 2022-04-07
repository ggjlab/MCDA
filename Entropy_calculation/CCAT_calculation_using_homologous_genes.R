

############### use hydra as an example
setwd("./CCAT_homo/hy")
library(reshape2)
library(Seurat)
library(Matrix)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(SCENT)
library(qlcMatrix)
library(ggplot2)
library(Seurat)
library(ggsignif)
library(ggpubr)
library(viridis)
options(stringsAsFactors = F)



################ 1. use homologous genes
### load the orthoFinder resuly
dduse <- read.csv("./homologousgenes/Hydra_human_all.csv",header = T)
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
tmp <- hydra_human[,colSums(hydra_human)<1]




################  2. prepare the hydra*cell matrix
### load the datasets
load('data/Hydra_dge_anno.RData')
cpm <- pbmc@assays$RNA@data
tmp <- rownames(cpm)
tmp <- colsplit(tmp,pattern = "-",names = c("a1","a2",'a3'))
length(unique(tmp$a1))
cpm <- cpm[1:33376,] ## delete the MT genes
tmp <- rownames(cpm)
tmp <- colsplit(tmp,pattern = "-",names = c("a1","a2",'a3'))
length(unique(tmp$a1))
rownames(cpm) <- tmp$a1

cpm1 <- cpm[intersect(rownames(cpm),rownames(hydra_human)),]
hydra_human1 <- hydra_human[intersect(rownames(cpm),rownames(hydra_human)),]
cpm1 <- Matrix(as.matrix(cpm1),sparse=T)
cpm1 <- t(cpm1)




################ 3. matrix multiplication to weight the hydra genes to human genes
library(Matrix)
hydra_human1 <- Matrix(as.matrix(hydra_human1),sparse=T)
cell_human <- cpm1 %*% hydra_human1
save(cell_human,file = "cell_human.RData")





################ 4. calculate CCAT
cell_human1 <- t(cell_human) 
gene<-mapIds(x=org.Hs.eg.db,keys =rownames(cell_human1),
             column='ENTREZID', keytype='SYMBOL', multiVals = "first" )
z<-data.frame(unlist(gene))
z<-na.omit(z)
colnames(z)<-"id"
cell_human1<-cell_human1[rownames(z),]
rownames(cell_human1)<-z$id

## modified CCAT function
CompCCAT1 <- function (exp.m, ppiA.m) 
{
  if (max(exp.m) > 100) {
    exp.m <- log2(exp.m + 1)
  }
  classMATRIX <- class(exp.m)
  commonEID.v <- intersect(rownames(ppiA.m), rownames(exp.m))
  print(length(commonEID.v))
  k.v <- rowSums(ppiA.m[match(commonEID.v, rownames(ppiA.m)), 
                        ])
  if (classMATRIX == "matrix") {
    ccat.v <- as.vector(cor(exp.m[match(commonEID.v, rownames(exp.m)), 
                                  ], k.v))
  }
  else if (classMATRIX == "dgCMatrix") {
    ccat.v <- as.vector(corSparse(exp.m[match(commonEID.v, 
                                              rownames(exp.m)), ], Matrix(matrix(k.v, ncol = 1))))
  }
  return(ccat.v)
}

ccat.v <- CompCCAT1(exp = cell_human1, ppiA = net13Jun12.m);
ccat <- data.frame(ccat_orthoFinder=ccat.v,row.names =colnames(cell_human1))
anno <- cbind(ccat,anno[rownames(ccat),])




################ 5. visulization
anno$stage <- anno$level
table(anno$stage)
stage<-unique(anno$stage)
stage <- paste0("level",1:5)
anno$stage <- factor(anno$stage,levels = stage)
compare_list <- function(stage_list){
  l1s <- list()
  m=1
  while (m<(length(stage_list))) {
    l1s[[m]] <- c(stage_list[m],stage_list[m+1])
    m=m+1
  }
  l1s[[m]] <- c(stage_list[1], stage_list[m])
  return(l1s)
}
comparsion <- compare_list(stage)
a1 <- max(anno$ccat_orthoFinder)
ggplot(anno,aes(stage,ccat_orthoFinder,fill=stage))+geom_violin(trim=F,color="grey",alpha=0.8)+
  geom_boxplot(width=0.2,color="grey",outlier.colour = "grey",outlier.size = 0.05,position = position_dodge(0.9))+scale_fill_viridis_d()+
  theme_classic()+
  stat_compare_means(method="wilcox.test",comparisons = comparsion,label = "p.format",label.y =c(rep(a1+0.01,length(comparsion)-1),a1+0.04) )+
  theme(axis.text.x = element_text( size = 15,angle = 90),
        axis.text.y = element_text( size = 15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),)+
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 0.5),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))

filename <- "Hydra_ccat_orthoFinder.pdf"
ggsave(filename,width = 8,height = 8)
write.csv(anno,file = "Hydra_ccat_orthoFinder_use.csv",quote = F)

