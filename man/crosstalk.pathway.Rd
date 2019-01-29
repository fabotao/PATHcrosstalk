\name{crosstalk.pathway}
\alias{crosstalk.pathway}
\title{
Quantify deregulation of pathways in cancer
}
\description{
Pathifier is an algorithm that infers pathway deregulation scores for each tumor sample on the basis of expression data. This score is determined, in a context-specific manner, for every particular dataset and type of cancer that is being investigated. The algorithm transforms gene-level information into pathway-level information, generating a compact and biologically relevant representation of each sample.
}
\usage{
crosstalk.pathway(PDS, pathway.names=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PDS}{
 The pathway deregulation scores (PDS) for each tumor sample inferred by Pathifier algorithm on the basis of expression data
}
  \item{pathway.names}{
A vector of pathways for corsstalk correction (subpathways after removing overlapping genes with another pathway)
}
}

\value{
%%  ~Describe the value returned
\item{gs}{The list of genes (by EntrezID) in each pathway}
\item{pathwaynames}{The KEGG map ID of the pathways}
}

\author{
Botao Fa <faber@sjtu.edu.cn>; Zhangsheng Yu
Maintainer: Botao Fa <faber@sjtu.edu.cn>
}

\examples{
library(PATHcrosstalk)
data(KEGG) # The 322 pathways of the KEGG database 
data(Sheffer, package='pathifier') # The colorectal data of Sheffer et al.
data(Symbol_EntrezID) # The Gene Symbol and EntrezID mapping data for convertion
genes <- sheffer$allgenes
ind <- which(genes %in% symbol.entrezid$GeneSymbol)
sheffer.map <- merge(data.frame(GeneSymbol=genes), symbol.entrezid, by='GeneSymbol', sort=F)

PDS<-quantify_pathways_deregulation(sheffer$data[ind,], sheffer.map$EntrezID,
  kegg$gs, kegg$pathwaynames, sheffer$normals, attempts = 100,
  logfile="sheffer.kegg.log", min_exp=sheffer$minexp, min_std=sheffer$minstd)
  
path.crosstalk <- crosstalk.pathway(PDS, pathway.names=c("00564", "00260", "00340", "00360", "00511", "04145", "04614", "04960", "05014", "05030", "05130"))

crosstalk.PDS<-quantify_pathways_deregulation(sheffer$data[ind,], sheffer.map$EntrezID,
  path.crosstalk$gs, path.crosstalk$pathwaynames, sheffer$normals, attempts = 100,
  logfile="sheffer.kegg.log", min_exp=sheffer$minexp, min_std=sheffer$minstd)
  
corss.mat <- crosstalk.mat.surv(PDS, crosstalk.PDS, sheffer$samples, sheffer$normals, pathway.names=c("00564", "00260", "00340", "00360", "00511", "04145", "04614", "04960", "05014", "05030", "05130"), surv.dat=sheffer.surv, by.surv='samples')

}