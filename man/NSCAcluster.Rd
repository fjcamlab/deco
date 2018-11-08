\name{NSCAcluster}
\alias{NSCAcluster}
\title{
NSCA-clustering function on an incidence matrix
}
\description{
This function provides a NSCA analysis and posterior biclustering (based on 'h statistic' or inner product from NSCA) of both features and samples from an incidence matrix (absolute frequency values).
}
\usage{
NSCAcluster(mx, data = NULL, id.names = NULL,
                k = NULL, label = NA,
                v = 80, method.dend = NULL,
                raw = NULL, UpDw = NULL, dir = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mx}{an incidence matrix containing absolute frequency values of any event.
}
  \item{data}{original omic data matrix. If any omic data is input, 'h statistic' would be computed and used as main statistic for biclustering analysis. Set as NULL by default, so inner product from NSCA is used to find out subclasses of samples.
}
  \item{id.names}{character vector whom names should correspond to rownames of incidenceMatrix to replace by original IDs (UP and DOWN labels were pasted if 'binary' analysis was carried out within DECO algorithm). If data is not related to DECO analysis, it should be NULL.
}
  \item{k}{numeric, it corresponds to number of subclasses to find out. Set as NULL by default, so this function will calculate better number of subclasses depending on Huber's Gamma coefficient. Further information could be found in our article.
}
  \item{label}{character indicating a label to paste in some informative slots for this group of samples.
}
  \item{v}{minimum variability threshold imposed to NSCA.
}
  \item{method.dend}{character indicating which agglomerative method should be used to generate sample dendrogram, from 'hclust' R function. If is NULL (by default), all possible methods would be tested and corresponding to the highest cophenetic correlation would be selected.
}
  \item{raw}{numeric vector indicating mean signal value from samples within other category compared if 'binary' analysis was carried out. This parameter is just used to 'h statistic' calculation.
}
  \item{UpDw}{character vector indicating direction of differential signal from samples within other category compared if 'binary' analysis was carried out. This parameter is just used to 'h statistic' calculation.
}
  \item{dir}{character indicating which features would be taken into account if 'binary' analysis was carried out. This parameter is just used to 'h statistic' calculation.
}
}
\value{
This function will return a list containing following slots:

\item{NSCA}{is a list with all information about NSCA analysis, where main results like Tau statistic and \emph{h statistic} matrix are contained in 'tau' and 'h' slots. If \emph{h statistic} was not calculated, a numeric matrix, called 'inner.prod', corresponding to original NSCA inner product would replace it.
}
\item{infoFeature}{is a data.frame containing all feature statistics derived from NSCA subclasses analysis.
}
\item{rankingFeature.h}{is a data frame containing mean \emph{h statistic} per subclass of samples, ranking for each subclass and \emph{h.Range} estimating discriminant power of each feature.
}
\item{infoSubclass}{is a summary about how many features and samples are included in each subclass.
}
\item{incidenceMatrix}{contains the original input matrix 'mx'.
}
\item{samplesSubclass}{data.frame which shows sample membership to a subclass.
}
\item{var}{numeric corresponding to the minimum variability explained by NSCA (defined by user as 'v').
}
\item{clus.method}{is the agglomeration method used to calculate the hierarchical clustering of samples.
}
\item{hclustSamp}{'hclust' R object from sample's dendrogram.
}
\item{hclustFeat}{'hclust' R object from feature's dendrogram.
}
\item{inner}{numeric matrix corresponding to original NSCA inner product result.
}
\item{dispersion}{numeric matrix corresponding to dispersion analysis for \emph{h statistic} calculation.
}
\item{Biplot.coordinates}{a matrix with two columns corresponding to first and second NSCA dimension, which disposes in rows all features (rows) and samples (columns) in order. First, we included all features and below all samples. It is intended to represent a 2D plot of NSCA relational space.
}
}
\references{
Lauro, N. and D'Ambra, L. (1984). \bold{L'analyse non symetrique des correspondances.} \emph{Data Analysis and Informatics}

Beh, E.J. and Lombardo, R. (2014). \bold{Correspondence Analysis. Theory, Practice and New Strategies}. \emph{John Wiley & Sons}
}
\author{
Francisco Jose Campos Laborie. <fjcamlab@usal.es>
}
\seealso{
\code{\link[stats]{hclust}}, \code{\link[deco]{decoRDA}}, \code{\link[deco]{decoNSCA}}
}
\examples{
#### Further information in DECO vignette(s).
## Run after decoRDA and decoNSCA steps on 'ALCLdata'
ALCLdata

# an example of absolute frequency matrix from ALCLdata analysis
# incidenceMatrix <- deco@NSCAcluster$Case$incidenceMatrix

# vector of names to replace.
# ids <- sapply(rownames(incidenceMatrix), function(x) unlist(strsplit(x,split = "_",
#              fixed = TRUE))[1])

# id.names <- deco@featureTable[match(ids, deco@featureTable[,c("ID")]),c("ID")]
# names(id.names) <- names(ids)

## NSCAcluster function
# res_nscaclus <- NSCAcluster(mx = incidenceMatrix, id.names = id.names,
#                             data = exprs(ALCLdata), k = NULL, label = NA,
#                             v = 80, method.dend = "ward.D")

}
