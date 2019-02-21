\name{AnnotateDECO}
\alias{AnnotateDECO}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
DECO annotating function
}
\description{
This function allows to annotate different ID types using \emph{Bioconductor} annotation packages (i.e. \emph{Homo.sapiens, Mus.musculus, org.Hs.eg.db, org.MeSH.Spo.972h.db, org.Mm.eg.db ...}).
}
\usage{
AnnotateDECO(ids, id.type, attributes = NA, pack.db = "Homo.sapiens",
                      rm.redundant = TRUE, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ids}{
a character vector including IDs to annotate.
}
  \item{id.type}{
a character indicating what type of ID we are managing. It should match with correct ID types including in annotation package ("Homo.sapiens" by default).
}
  \item{attributes}{
a character indicating which types of ID should annotate the original IDs. It should match with correct ID types including in annotation package (\emph{Homo.sapiens} by default).
}
  \item{pack.db}{
name of the Bioconductor annotation package to use.
}
  \item{rm.redundant}{
logical indicating if redundant annotation (i.e. two different HGNC symbols for the same EntrezID) should be removed or not.
}
  \item{verbose}{
logical indicating if messages should be displayed.
}
}
\value{
This function returns a \code{data.frame} containing original IDs and corresponding attributes to each one.
}
\author{
Francisco Jose Campos Laborie <fjcamlab@usal.es>
}
\examples{
##### Annotation example
# Select an appropiate 'id.type' depending on which annotation
# package will be used.
# library(Homo.sapiens)
# columns(Homo.sapiens)

# Selecting first 10 original IDs from ALCL dataset
data(ALCLdata)
ids <- rownames(assay(ALCL))[1:10]

## Not run as example
# annot <- AnnotateDECO(ids, id.type = "ENSEMBL",
#                       attributes = c("SYMBOL","CHR","GENENAME"),
#                       pack.db = "Homo.sapiens", rm.redundant = TRUE)

}

