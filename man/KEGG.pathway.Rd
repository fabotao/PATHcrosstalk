\name{KEGG}
\alias{KEGG}
\docType{data}
\title{
The 322 pathways of the KEGG database
}
\description{
The 322 pathways of the KEGG database
}
\usage{data(KEGG.pathway)}
\format{
  \describe{
    \item{\code{pathwaynames}}{The KEGG map ID of the pathways}
    \item{\code{gs}}{The list of genes (by EntrezID) in each pathway}
	\item{\code{pathwayfunction}}{The function (KEGG) of the pathways}
  }
}

\source{
Kanehisa M, Goto S, Sato Y, Furumichi M and Tanabe M. KEGG for integration and interpretation of large-scale molecular datasets. \emph{Nucleic Acids Res}, 2012, Vol 40(Database issue):D109-D114.
}
\examples{
data(KEGG)
}
\keyword{datasets}
