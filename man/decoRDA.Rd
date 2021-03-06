\name{decoRDA}
\alias{decoRDA}
\title{
Subsampling function to find out differential events along samples from an omic dataset
}
\description{
RDA R function within 'deco' R package. This function allows user to subsampling a omic data matrix to find out significant differential features characterizing any hidden subclass or group of samples included.}
\usage{
decoRDA(data, classes = NA, control = NA, r = NULL, q.val = 0.01,
              iterations = 10000, bpparam = SerialParam(), annot = FALSE, 
              id.type = NA, attributes = NA, rm.xy = FALSE, 
              pack.db = "Homo.sapiens")
}
\arguments{
  \item{data}{
matrix of normalized omic data with 'f' features (rows) by 'n' samples (columns).
}
  \item{classes}{
factor or character vector of 's' size with two ('Supervised' or 'Binary' analysis) or more classes ('Multiclass' analysis) to contrast. Names of vector must match with colnames of 'data'. If any vector of classes is given, 'Unsupervised' analysis will be run with all samples.
}
  \item{control}{
character label defining 'control' or '0' class for LIMMA contrast design.
}
  \item{r}{
number of samples of each class to take for subsampling. Same number of samples will be taken from each class.
}
  \item{q.val}{
p.value threshold for each iteration of subsampling process. It corresponds to \emph{adj.P.value} of LIMMA and is set as default as 0.01.
}
  \item{iterations}{
number of iterations to subsample 'data', set as default as 10000. If this value exceeds number of all possible combinations, it would be replaced.
}
  \item{bpparam}{registered parallelization using BiocParallel R functions: SerialParam(), MulticoreParam(), SnowParam(), etc. Further information available in vignette.
}
  \item{annot}{
logical value indicating if annotation of IDs provided by 'data' as rownames should be done.
}
  \item{id.type}{
character indicating what ID is provided by user in 'data'. This value must match any of distinct possibles values of R package used to annotate.
}
  \item{attributes}{
character vector with distinct fields to return by annotation package.
}
  \item{rm.xy}{
logical indicating if features located on chromosomes 'X' or 'Y' have to be removed. Features depending on these location could be an artefact of unbalanced contrasts by gender of patients or samples, so removing it should clear final results due to an appropiated subsampling of samples.
}
  \item{pack.db}{
character naming annotation package to be used. \emph{Homo.sapiens} annotation package set as default, try first 'AnnotateDECO' R function with any other annotation package from Bioconductor to find out any problem.
}}
\details{
The \bold{RDA step} is primarly conditioned by contrast design. Our capacity to highlight \emph{majority} or \emph{minority} features will vary depending on how much we want focusing the analyisis on classes or individual samples, or what is the same, the \bold{granularity of RDA}.

In this way, it is necessary to highlight that two \bold{RDA strategies} could be: (a) a \bold{majority} one trying to improve our generic differential expression with higher subsampling size; and (b) a \bold{minority} strategy to find out all possible hidden subclasses. By default, \code{decoRDA()} defines an optimal subsample size \emph{r}. Let \emph{n} be the number of samples included in the analysis or \deqn{n_1,n_2,...,n_k} in case of two or more \emph{k} classes, then based on previous analysis done by Babu et al.: \deqn{r = \sqrt{n}}
or in case of 'Supervised' or 'Multiclass' design: \deqn{r = \sqrt{min(n_1,n_2,...,n_k)}}

The analysis can be set to find differences among two populations or classes and to find differences within the whole cohort of samples. Besides, \code{decoRDA()} can work under both scenarios: \bold{binary} analysis with two classes or \bold{unsupervised} analysis contrasting all samples. The \emph{classes} input is used to define it. If the user introduces any vector with labels of two classes for each sample, a \bold{binary} analysis will be run. Otherwise, without any \emph{classes} vector all samples will be compared with each other.
}
\value{
Returns a list containing:
  \item{data}{
input matrix of normalized data with 'f' features (rows) by 's' samples (columns).
}
  \item{subStatFeature}{
table of differential features with statistical and annotation information, if it was required.
}
  \item{incidenceMatrix}{
absolute frequency matrix of 'd' features by 's' samples size that summarizes differential events for each feature per sample. It will be essential for Non-Symmetrical Correspondence Analysis.
}
  \item{classes}{
factor or character vector of 's' size provided by user.
}
  \item{resampleSize}{
number of samples of each class to take for subsampling provided by user in 'r'.
}
  \item{control}{
character label defining 'control' or '0' class for LIMMA contrast design.
}
  \item{pos.iter}{
number of subsampling iterations showing at least 1 differential event.
}
  \item{q.val}{
\emph{adj.p.value} threshold used in each LIMMA iteration within subsampling.
}
}
\author{
Francisco Jose Campos Laborie. <fjcamlab@usal.es>
}
\seealso{
\code{\link[deco]{decoNSCA}}, \code{\link[limma]{voom}}
}
\examples{
#### ALCL EXAMPLE (Scarfo et al., 2015. Blood) ####

########################
# Loading example data #
########################
## Data from two subtypes (ALK+ and ALK-) of Anaplastic Large Cell Leukemia (ALCL).
data(ALCLdata)

## Classes vector to run a binary analysis to compare both classes.
classes.ALCL <- colData(ALCL)[,"Alk.positivity"]
names(classes.ALCL) <- colnames(ALCL)

#######################################################################
# RUNNING SUBSAMPLING OF DATA: BINARY design (two classes of samples) #
#######################################################################
# library(Homo.sapiens) # for gene annotation

## Not run as example
# sub.ma.3r.1K <- decoRDA(data = assay(ALCL), classes = classes.ALCL, q.val = 0.01,
#               rm.xy = TRUE, r = NULL, control = "pos", annot = TRUE, 
#               id.type = "ENSEMBL", iterations = 1000,
#               pack.db = "Homo.sapiens")


## Slots included in returned R object.

# names(sub.ma.3r.1K)
# [1] "data"            "results"         "subStatFeature"  "incidenceMatrix" "classes"
# [6] "resampleSize"    "control"         "pos.iter"        "q.val"           "call"


## Top-10 RDA features.
# head(sub.ma.3r.1K$subStatFeature, 10)


#################################################################
# RUNNING SUBSAMPLING OF DATA: UNSUPERVISED design (no classes) #
#################################################################

## Not run as example
# sub.ma.3r.1K.uns <- decoRDA(data = assay(ALCL), q.val = 0.01,
#                      rm.xy = TRUE, r = 3, annot = TRUE, 
#                      id.type = "ENSEMBL", iterations = 1000,
#                      pack.db = "Homo.sapiens")



#################################################################
# RUNNING SUBSAMPLING OF DATA: MULTICLASS design (no classes) #
#################################################################

# 3 classes: ALK+, PTCL and ALK- without PTCLs
multiclasses.ALCL <- factor(apply(
    as.data.frame(colData(ALCL)[, c("Alk.positivity", "Type")]), 1,
    function(x) paste(x, collapse = ".")
))
head(multiclasses.ALCL)

## Not run as example
# sub.ma.3r.1K.multi <- decoRDA(data = assay(ALCL), classes = multiclasses.ALCL,
#                       q.val = 0.01, rm.xy = TRUE, r = 3, annot = TRUE, 
#                       id.type = "ENSEMBL", iterations = 1000,
#                       pack.db = "Homo.sapiens")
}
