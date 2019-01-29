\name{symbol.EntrezID}
\alias{symbol.EntrezID}
\docType{data}
\title{
Gene symbol and EntrezID mapping data
}
\description{
Gene symbol and EntrezID mapping data in KEGG data, aiming to convert gene identifiers between gene symbol and EntrezID when calculating PDS using Pathifier
}
\usage{data(symbol.EntrezID)}
\format{
  \describe{
    \item{\code{GeneSymbol}}{The gene symbol of genes in KEGG data}
    \item{\code{EntrezID}}{The EntrezID of the genes in KEGG data}
  }
}

\source{
Kanehisa M, Goto S, Sato Y, Furumichi M and Tanabe M. KEGG for integration and interpretation of large-scale molecular datasets. \emph{Nucleic Acids Res}, 2012, Vol 40(Database issue):D109-D114.
}
\examples{
data(symbol.EntrezID)
}
\keyword{datasets}
