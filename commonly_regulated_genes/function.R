

#### tongyi the name of cell type
check_pairname <- function(marker_all_all){
  
  marker_all_all$cluster <- gsub(marker_all_all$cluster,pattern = " ",replacement = ".")
  marker_all_all$cluster <- gsub(marker_all_all$cluster,pattern = "\\(",replacement = ".")
  marker_all_all$cluster <- gsub(marker_all_all$cluster,pattern = "\\)",replacement = ".")
  marker_all_all$cluster <- gsub(marker_all_all$cluster,pattern = "\\-",replacement = ".")
  marker_all_all$cluster <- gsub(marker_all_all$cluster,pattern = "\\&",replacement = ".")
  marker_all_all$cluster <- gsub(marker_all_all$cluster,pattern = "/",replacement = ".")
  marker_all_all$cluster <- gsub(marker_all_all$cluster,pattern = "\\+",replacement = ".")
  marker_all_all$cluster <- gsub(marker_all_all$cluster,pattern = "/",replacement = ".")
  return(marker_all_all)
  
}


marker_to_df <- function(marker_all,cate){
  if(cate=="negative"){
    marker11_n<-marker_all[marker_all$avg_log2FC<0,]
  }else if( cate=="positive"){
    marker11_n<-marker_all[marker_all$avg_log2FC>0,]
  }
  marker11_n <- marker11_n[marker11_n$p_val_adjbh<=0.1,]
  df_n_all<-dcast(marker11_n,gene~cluster,value.var="avg_log2FC")
  a1<-colnames(df_n_all)[-1]
  df_n_all<-data.frame(df_n_all[,-1],row.names = df_n_all$gene)
  colnames(df_n_all)<-a1
  df_n_all[is.na(df_n_all)]<-0
  
  return(df_n_all)
}


marker_to_genesum <- function(df_n,cate){

  if(cate=="negative"){
    gene_sum<-data.frame(Gene=rownames(df_n),num=rowSums(df_n<0))
    gene_sum<-gene_sum[order(gene_sum$num,decreasing = T),]
    
  }else 
    if( cate=="positive"){
      gene_sum<-data.frame(Gene=rownames(df_n),num=rowSums(df_n>0))
      gene_sum<-gene_sum[order(gene_sum$num,decreasing = T),]
      
    }
  return(gene_sum)
}


Orth_Convert <- function(x){
  orth<-read.csv(orthfile,header = F,stringsAsFactors = F)
  #orth <- data.frame(V1=orth$V5,V2=orth$V4)
  orth$V2<-gsub(orth$V2,pattern ="_" ,replacement = "-")
  orth_target<-orth[orth$V2%in%target_gene$Gene,]
  orth_target1<-target_gene[target_gene$Gene%in%orth_target$V2,]
  orth_target<-orth_target[,1:2]
  colnames(orth_target)<-c("Human_id","Gene")
  orth_target<-merge(orth_target,orth_target1,by="Gene",all = T)
  #orth_target<-orth_target[order(orth_target$num,decreasing = T),]
  return(orth_target) 
}
