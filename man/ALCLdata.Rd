\name{ALCLdata}
\alias{ALCLdata}
\docType{data}
\title{
Subset of microarray gene expression data from Anaplastic Large Cell Leukemia disease.
}
\description{
RMA normalized microarray gene expression data from clinical samples of Anaplastic Large Cell Leukemia (ALCL) patients. ALCL subtypes ALK+ (n=11) and ALK- (n=20). Corresponding to GSE65823 (GEO database ID) on Affymetrix HGU133Plus2.0 platform mapped to ENSEMBL genes with \emph{genemapperhgu133plus2cdf} CDF package from \bold{GATEexplorer} website. RMA (Robust Multi-Array Average) normalization of microarrays and filtering of genes were applied (from 20172 genes to 1000 genes). Used as example of group-VS-group (binary) comparison with outlier behavior within ALK- samples.
}
\usage{data(ALCLdata)}
\format{
Formal class 'ExpressionSet' [package "Biobase"] with 7 slots.
}
\value{
A matrix (1000 rows x 31 columns) containing the gene expression data (RMA normalized microarray data).
}
\source{
GEO database ID: GSE65823

\url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65823}

GATExplorer CDF genemappers:

\url{https://doi.org/10.1186/1471-2105-11-221}
}
\references{
Scarfo, Irene et al. \bold{Identification of a new subclass of ALK negative ALCL expressing aberrant levels of ERBB4 transcripts.} \emph{Blood} (2015). \url{http://dx.doi.org/10.1182/blood-2014-12-614503}.
}
\examples{
data(ALCLdata)

## list objects included
ls()
# [1] "ALCLdata"        "sub.ma.3r.1K"    "deco.results.ma"

# phenodata of each sample
head(colData(ALCL))

# classes vector
classes.ALCL <- colData(ALCL)[,"Alk.positivity"]
names(classes.ALCL) <- colnames(ALCL)

# gene expression matrix with 1000 ENSEMBL genes (mapped using GATExplorer genemappers)
# and 31 samples (11 ALK+ and 20 ALK-)
dim(assay(ALCL))
# [1] 1000    31

# density plot of RMA gene values from all dataset
plot(density(assay(ALCL)))

}
\keyword{dataset}
\keyword{transcriptomic}
