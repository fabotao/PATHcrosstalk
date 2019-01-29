\name{crosstalk.mat.surv}
\alias{crosstalk.mat.surv}
\title{
Quantify deregulation of pathways in cancer
}
\description{
Pathifier is an algorithm that infers pathway deregulation scores for each tumor sample on the basis of expression data. This score is determined, in a context-specific manner, for every particular dataset and type of cancer that is being investigated. The algorithm transforms gene-level information into pathway-level information, generating a compact and biologically relevant representation of each sample.
}
\usage{
crosstalk.mat.surv(PDS, crosstalk.PDS, samples, normals, pathway.names, surv.dat, by.surv=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PDS}{
 The pathway deregulation scores (PDS) for each tumor sample inferred by Pathifier algorithm on the basis of expression data
}
  \item{crosstalk.PDS}{
 The pathway deregulation scores (PDS) of pathways with overlapping genes after crosstalk corection for each tumor sample inferred by Pathifier algorithm on the basis of expression data
}
  \item{samples}{
 The sample identifier for each tumor sample in expression data
}
  \item{normals}{
 A list of m logicals, true if a normal sample, false if tumor.	
}
  \item{pathway.names}{
 A vector of KEGG pathways (Map ID) which will be used to identify subpathways after crosstalk correction. 
}
  \item{surv.dat}{
 A data.frame consists of tumor sample identifier and corresponding survival data (time (length of follow up), event(death status)) of the patients where tumor samples were obtained.  
}
  \item{by.surv}{
 Specifications of the columns in surv.dat data.frame used for merging with pathway deregulation score extracted.
}
}

\value{
%%  ~Describe the value returned
\item{p.matrix}{The crosstalk matrix consisting of p.values representing the significance of subpathways after crosstalk correction. This matrix can be plotted with heatmap functions where each cell [i,j] represent the significance (eg. survival significance) of pathways i after removing gene-overlap with pathway j}
}

\author{
Botao Fa <faber@sjtu.edu.cn>; Zhangsheng Yu
Maintainer: Botao Fa <faber@sjtu.edu.cn>
}

\examples{
library(PATHcrosstalk)
library(pathifier)
data(KEGG, package='PATHcrosstalk') # The 322 pathways of the KEGG database 
data(Sheffer, package='pathifier') # The colorectal data of Sheffer et al.
data(symbol.EntrezID) # The Gene Symbol and EntrezID mapping data for convertion
genes <- sheffer$allgenes
gene.df=data.frame(GeneSymbol=genes, stringsAsFactors=FALSE)
sheffer.map <- merge(gene.df, symbol.entrezid, by='GeneSymbol', sort=FALSE)
ind <- which(genes %in% symbol.entrezid$GeneSymbol)

PDS<-quantify_pathways_deregulation(sheffer$data[ind,], sheffer.map$EntrezID,
  kegg$gs, kegg$pathwaynames, sheffer$normals, attempts = 100,
  logfile="sheffer.kegg.log", min_exp=sheffer$minexp, min_std=sheffer$minstd)
  
path.crosstalk <- crosstalk.pathway(PDS, pathway.names=c("00564", "00260", "00340", "00360", "00511", "04145", "04614", "04960", "05014", "05030", "05130"))

crosstalk.PDS<-quantify_pathways_deregulation(sheffer$data[ind,], sheffer.map$EntrezID,
  path.crosstalk$gs, path.crosstalk$pathwaynames, sheffer$normals, attempts = 100,
  logfile="sheffer.kegg.log", min_exp=sheffer$minexp, min_std=sheffer$minstd)
  
corss.mat <- cross.mat.surv(PDS, crosstalk.PDS, sheffer$samples, sheffer$normals, pathway.names=c("00564", "00260", "00340", "00360", "00511", "04145", "04614", "04960", "05014", "05030", "05130"), surv.dat=sheffer.surv, by.surv='samples')

}