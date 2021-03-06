\name{plotHeatmapH}
\alias{plotHeatmapH}
\title{
Heatmap corresponding to h-statistic matrix provided by DECO analysis.
}
\description{
This function returns three plots showing the association among different sample subclasses found by DECO and any sample information of interest. This information can be the initial classes or any new included sample information. 
}
\usage{
plotHeatmapH(deco, info.sample = NA, info.feature = NA, 
      print.annot = FALSE, cex.legend = 1, 
      cex.names = 1)
}
\arguments{
  \item{deco}{
a "deco" R object generated by 'decoNSCA' R function.
}
  \item{info.sample}{
a data.frame with relevant information of samples to be plot on 'Heatmap'.
}
  \item{info.feature}{
a data.frame with relevant information of genes or features to be plot on 'Heatmap'.
}
  \item{print.annot}{
logical indicating if annotation (gene SYMBOL) should be printed instead original IDs. This option needs previous annotation within 'decoRDA' function or late annotation binding information to @featureTable slot of 'deco' R object.
}
  \item{cex.legend}{
numerical value giving the amount of magnification of legend.
}
  \item{cex.names}{
numerical value giving the amount of magnification to gene names or sample names.
}
}
\details{
Heatmap plot showing h-statistic matrix (H) and biclustering of features-samples.
}
\value{
A customized heatmap plot showing h-statistic matrix and biclustering of features-samples.
}
\author{
Francisco Jose Campos Laborie. <fjcamlab@gmail.com>
}
\seealso{
\code{\link[deco]{decoReport}}, \code{\link[deco]{decoNSCA}}
}
\examples{
#### Further information in DECO vignette(s).
## Run after decoRDA and decoNSCA steps on 'ALCLdata'
data(ALCLdata)
ALCL

### Phenotype information
info.sample <- as.data.frame(colData(ALCL))[,8:9]
rownames(info.sample) <- rownames(colData(ALCL))

#### Heatmap with h-statistic matrix and biclustering of features-samples.
## Not run as example
# plotHeatmapH(deco = deco.results.ma, info.sample = info.sample, 
#              cex.names = 0.3, print.annot = FALSE)
}

