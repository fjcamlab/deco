\name{featureTable}
\alias{featureTable}
\title{
Accessor to 'featureTable' slot of 'deco-class' objects.}
\usage{
featureTable(object)
}
\arguments{
  \item{object}{a deco-class object.
}
}
\description{
The 'featureTable' slot is the main output table, including all the feature statistics and rankings from both steps of DECO: RDA and NSCA.
}
\value{
This function returns the featureTable within a deco-class object.
}
\author{
Francisco Jose Campos Laborie. <fjcamlab@gmail.com>
}
\seealso{
\code{\link[deco]{decoRDA}}
}
\examples{
## Loading pre-calculated objects...
data(ALCLdata)

featTable <- featureTable(deco.results.ma)

head(featTable)

}
