\name{PATHcrosstalk-package}
\alias{PATHcrosstalk-package}
\alias{PATHcrosstalk}
\docType{package}
\title{
Identify deregulation of pathways in cancer after crosstalk effect acommodated 
}
\description{
Pathway information is always treated as independent, however, pathways can crosstalk with each other by sharing overlapping genes. This effect can complicate the identification of pathway-level biomarkers in statistical models like machine learning. Here we disscuss this effect based on the algorithm of Pathifier, which infers pathway deregulation scores for each tumor sample on the basis of expression data. This  algorithm transforms gene-level information into pathway-level information, generating a compact and biologically relevant representation of each sample. With PATHcrosstalk, we can evaludate the deregulation of subpathways after crosstalk correction (removing overlapping genes with other pathways).
}
\details{
\tabular{ll}{
Package: \tab PATHcrosstalk\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2019-01-22\cr
License: \tab Artistic-1.0\cr
}
}
\author{
Botao Fa <faber@sjtu.edu.cn>; Zhangsheng Yu
Maintainer: Botao Fa <faber@sjtu.edu.cn>
}
\references{

}

\keyword{ package }

\examples{
library(PATHcrosstalk)
library(pathifier)
data(KEGG, package='PATHcrosstalk') # The 322 pathways of the KEGG database 
data(Sheffer, package='pathifier') # The colorectal data of Sheffer et al.
data(symbol.EntrezID) # The Gene Symbol and EntrezID mapping data for convertion
data(sheffer.surv)  # Survival information of the colorectal data of Sheffer et al.

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
  
cross.mat <- cross.mat.surv(PDS, crosstalk.PDS, sheffer$samples, sheffer$normals, pathway.names=c("00564", "00260", "00340", "00360", "00511", "04145", "04614", "04960", "05014", "05030", "05130"), surv.dat=sheffer.surv, by.surv='samples')
}
