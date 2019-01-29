
## Generate subpathways after removing overlapping genes
# between pathway pairs. 
crosstalk.pathway <- function(PDS, pathway.names=NULL){
  genes <- PDS$genesinpathway
  path.all <- names(genes)
  
  # If the pathways input is null, use all pathways 
  # in PDS to construct pathways with crosstalk accommodated
  if(is.null(pathway.names)){
    pathway.names <- path.all
    print('Use all pathways in PDS data.')
  }else if(length(pathway.names)<2){
    stop('Not enough pathways, please provide at least two pathways!')
    
  }else if(length(setdiff(pathway.names, path.all))>=1){
    print('Pathways input is out of the range of the PDS data provided!')
    print('Use the common pathways...')
    pathway.names <- intersect(pathway.names, path.all)
  }

  len <- length(pathway.names)
  pair.all <- combn(1:len, 2)

  gs <- c()
  pathwaynames <- c()
  for(i in 1:dim(pair.all)[2]){

    path1 <- pathway.names[pair.all[1,i]]
    path2 <- pathway.names[pair.all[2,i]]
    id1 <- which(path.all %in% path1)
    id2 <- which(path.all %in% path2)
    gene1 <- genes[[id1]]
    gene2 <- genes[[id2]]
    gene.both <- intersect(gene1, gene2)
    gene.all <- union(gene1, gene2)
    if(length(gene.both)>=3){
      p.name1 <- paste0(path1,'_diff_',path2)
      p.name2 <- paste0(path2,'_diff_',path1)
      p.name.inter <- paste0(path1,'_intersect_',path2)
      
      gs <- c(gs, list(t(t(setdiff(gene1, gene2)))), list(t(t(setdiff(gene2, gene1)))), list(t(t(gene.both))))
      pathwaynames <- c(pathwaynames, list(p.name1), list(p.name2), 
                        list(p.name.inter))
    }
    
    
  }
  pathway.dat <- list()
  pathway.dat$gs <- gs
  pathway.dat$pathwaynames <- pathwaynames
  print(paste0('Pathway No. is: ', length(pathway.dat$gs)))
  return(pathway.dat)
  
}



