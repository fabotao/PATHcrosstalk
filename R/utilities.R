## Extract and reorganize score into dataframe for analysis from PDS Rdata
# Input: 
# pds: Object from quantify_pathways_deregulation(); 
# data: the n x m mRNA expression data matrix for PDS calculation, where n is the number of genes and m the number of samples

extract.score <- function(PDS, samples){
  
  score <- PDS$scores
  df <- data.frame(sample=samples, stringsAsFactors = F)
  for(pathway in names(score)){
    df[,pathway] <- c(score[[pathway]])
  }
  return(df)
}



## Organize crosstalk matrix for heatmap illustration
crosstalk.mat.surv <- function(PDS, crosstalk.PDS, samples, normals, pathway.names, surv.dat, by.surv=NULL){
  # Extract PDS score
  score <- extract.score(PDS, samples)
  score.tumor <- score[!normals,]
  cross.score <- extract.score(cross.PDS, samples)
  cross.score.tumor <- cross.score[!normals,]
  
  score.all <- merge(score.tumor, cross.score.tumor, by='sample')
  
  data <- merge(surv.dat, score.all, by.x=by.surv, by.y='sample')
  data <- data[,2:dim(data)[2]]
  
  time=data[,1]
  fail=data[,2]
  
  p.func <- function(dat){
    p.vals <- c()
    for(i in 3:length(dat)){
      uni.mod=coxph(Surv(time,fail)~dat[,i])
      p.val <- (summary(uni.mod)$coef)[,'Pr(>|z|)']
      p.vals <- c(p.vals, p.val)
    }
    
    names(p.vals) <- names(dat)[3:length(dat)]
    return(p.vals)
  }
  p.all <- p.func(data)
  
  # 排序整合P值矩阵
  p.mat <- matrix(NA, ncol=length(pathway.names), nrow=length(pathway.names))
  
  p.orign <- p.all[pathway.names]
  p.orign <- p.orign[order(p.orign, decreasing=F)]
  
  for(i in 1:length(p.orign)){
    for(j in 1:length(p.orign)){
      i_inter_j <- paste0(names(p.orign)[i], '_intersect_', names(p.orign)[j])
      j_inter_i <- paste0(names(p.orign)[j], '_intersect_', names(p.orign)[i])
      i_diff_j <- paste0(names(p.orign)[i], '_diff_', names(p.orign)[j])
      if(i==j){
        p.mat[i,j] <- p.orign[i]
      }else if(i_diff_j %in% names(p.all)){
        p.mat[i,j] <- p.all[i_diff_j]
      }else{
        p.mat[i,j] <- NA
      }
    }
  }
  p.matrix <- p.mat
  
  dimnames(p.matrix) <- list(names(p.orign), names(p.orign))
  return(p.matrix)
}




