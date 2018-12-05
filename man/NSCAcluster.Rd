\name{NSCAcluster}
\alias{NSCAcluster}
\title{
Accessor to 'NSCAcluster' slot of 'deco-class' objects.}
\usage{
NSCAcluster(object)
}
\arguments{
  \item{object}{a deco-class object.
}
}
\description{
The 'NSCAcluster' slot is the main output from NSCA step. I would be divided in 'Control' and 'Case' if only two categories of samples were input to first 'decoRDA()' function, or would include only one 'All' slot for multiclass or unsupervised analysis. 
}
\value{
This function returns the NSCAcluster list within a deco-class object.
}
\author{
Francisco Jose Campos Laborie. <fjcamlab@gmail.com>
}
\seealso{
\code{\link[deco]{decoRDA}}
}
\examples{
##
## Loading pre-calculated objects...
data(ALCLdata)

resNSCA <- NSCAcluster(deco.results.ma)

head(resNSCA)

# Information within the 'case' category
str(resNSCA$Case)

}
