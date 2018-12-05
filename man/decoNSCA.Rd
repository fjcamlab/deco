\name{decoNSCA}
\alias{decoNSCA}
\title{
Non-Symmetrical Correspondence Analysis (NSCA) applied to a subsampled dataset
}
\description{
This function needs R object from 'decoRDA' function to operate. It will applied NSCA to the frequency matrix of differential events, so it could calculate 'h' statistic per feature per samples to find out new subclasses of samples.
}
\usage{
decoNSCA(sub, v = 80, k.control = NULL,
         k.case = NULL, rep.thr = 3,
         bpparam = SerialParam(),
         samp.perc = 0.05,
         method = "ward.D")
}
\arguments{
  \item{sub}{an output R object from 'decoRDA' R function.
}
  \item{v}{minimum percentage of variability that NSCA have to explain.
}
  \item{k.control}{number of subclasses to find out within 'Control' samples (binary contrast) or all samples (unsupervised contrast).
}
  \item{k.case}{number of subclasses to find out within 'Case' samples (binary contrast).
}
  \item{rep.thr}{numeric, it should correspond to the minimum number of differential events that should have each sample within 'samp.perc' percentage at least. Both 'rep.thr' and 'samp.perc' compounds a threshold to remove noisy features not explicity associated to any subgroup of samples.
}
  \item{bpparam}{registered parallelization using BiocParallel R functions: SerialParam(), MulticoreParam(), SnowParam(), etc. Further information available in vignette.
}
  \item{samp.perc}{numeric, it should correspond to the minimum percentage of samples having to be affected with a minimum 'rep.thr' number of differential events per feature. Both 'rep.thr' and 'samp.perc' compounds a threshold to remove noisy features not explicity associated to any subgroup of samples.
}
  \item{method}{character indicating which agglomerative method should be used to generate sample dendrogram. If is NULL (by default), all possible methods would be tested and corresponding to the highest cophenetic correlation would be selected.
}
}
\details{
Once we obtain frequency matrix of DE events or \emph{incidenceMatrix}, we applied a \bold{NSCA} procedure. NSCA let us to analyse all dependendent structures between differential features and samples derived from the same relational space (dataset).

In this context, NSCA stablishes an assymetric and directional association between features (response) and samples (predictor), so each sample contributes in a particular way to the differential expression of a singular feature. Thus, an outlier subsets of samples for a feature could explain most of its differential expression signal, allowing us to assign this feature with these samples (see figure of \emph{incidenceMatrix} above). Further information about how NSCA functions and other properties are detailed by Lombardo et al.


}
\value{
Returns a "deco" object with following slots:
  \item{featureTable}{
table of results with gene or feature statistical information about NSCA merged with subsampling information contained in previous 'diffgenes' value from 'decoRDA' function results.
}
  \item{NSCAcluster}{
list from 'NSCAcluster' R function containing information about NSCA. If 'Binary' analysis was previously set, two lists corresponding to each class are contained.
}
  \item{incidenceMatrix}{
absolute frequency matrix of 'd' genes by 's' samples size that summarizes DE events for each gene per sample. It will be essential for Non-symmetrical correspondence analysis.
}
  \item{classes}{
factor or character vector of 's' size provided by user.
}
  \item{control}{
label if control class was defined.
}
  \item{pos.iter}{
lowest number of repeats after applying 'decoNSCA' threshold.
}
  \item{q.val}{
adjusted.p.value threshold defined to each subsampling iteration with LIMMA within 'decoRDA' function.
}
  \item{subsampling.call}{
call to decoRDA previously run.
}
  \item{deco.call}{
call to decoNSCA previously run.
}
}
\references{
Lauro, N. and D'Ambra, L. (1984). \bold{L'analyse non symetrique des correspondances.} \emph{Data Analysis and Informatics}

Beh, E.J. and Lombardo, R. (2014). \bold{Correspondence Analysis. Theory, Practice and New Strategies}. \emph{John Wiley & Sons}
}
\author{
Francisco Jose Campos Laborie. <fjcamlab@gmail.com>
}
\seealso{
\code{\link[deco]{decoRDA}}, \code{\link[stats]{hclust}}
}
\examples{
## User has to provide a RDA R object, running decoRDA() previously
data(ALCLdata)

## Slots included in the returned RDA object are:
names(sub.ma.3r.1K)
# [1] "data"            "results"         "subStatFeature"  "incidenceMatrix" "classes"
# [6] "resampleSize"    "control"         "pos.iter"        "q.val"           "call"

# Computing in shared memory
# all cores by deault
bpparam <- MulticoreParam()

#########################################################################################
# RUNNING NSCA STEP: Looking for subclasses within a category/class of samples compared #
#########################################################################################

deco.results.ma <- decoNSCA(sub = sub.ma.3r.1K, v = 80, 
        method = "ward.D", bpparam = bpparam,
        k.control = 3, k.case = 3, samp.perc = 0.05, rep.thr = 5)


## Class 'deco'
class(deco.results.ma)
# [1] "deco"


## Slots included within 'deco' R object

# slotNames(deco.results.ma)
# [1] "featureTable"     "NSCAcluster"      "incidenceMatrix"  "classes"          "pos.iter"
# [6] "control"          "q.val"            "rep.thr"          "samp.perc"        "subsampling.call"
#[11] "nsca.call"


## Top-10 features from DECO analysis based on "Standard.Chi.Square"
head(featureTable(deco.results.ma)[order(featureTable(deco.results.ma)$Standard.Chi.Square,
     decreasing = TRUE),], 10)


## Matrix of 'h statistic' calculated by DECO per category of samples
# if binary analysis was carried out.
dim(NSCAcluster(deco.results.ma)$Control$NSCA$h)
dim(NSCAcluster(deco.results.ma)$Case$NSCA$h)


## Top-10 discriminant features for CASE samples (ALK-) based on
# 'h statistic per subclass found.
head(NSCAcluster(deco.results.ma)$Case$rankingFeature.h, 10)


### Sample and subclass information could be found in two slots within 'NSCAcluster':
## General information about subclasses
# NSCAcluster(deco.results.ma)$Case$infoSubclass

## Sample membership to a subclass.
# NSCAcluster(deco.results.ma)$Case$samplesSubclass


## Both 'hclust' dendrogram information of CASE samples and features.
## They include Huber's Gamma coefficient value and Cophenetic
## correlation between dendrogram distances
# and distance matrix.

#names(NSCAcluster(deco.results.ma)$Case$hclustSamp)
#[1] "dend"    "coph"    "cluster" "huber"

#names(NSCAcluster(deco.results.ma)$Case$hclustFeat)
#[1] "dend"    "coph"    "cluster" "huber"

###### Get summary info about "DECO" analysis:
###
summary(deco.results.ma)
# Decomposing Heterogeneous Cohorts from Omic profiling: DECO
# Summary:

# Analysis design: Binary
# Classes compared:
# neg pos
#  20  11

#            RDA.q.value Minimum.repeats Percentage.of.affected.samples NSCA.variability
# Thresholds        0.01           10.00                           5.00             84.1

# Number of features out of thresholds: 255
# Feature profile table:
# Complete Majority Minority
#       12       79      164
# Number of samples affected: 31
# Number of positive RDA comparisons: 1955
# Number of total RDA comparisons: 10000

###### Get info about "DECO" analysis
show(deco.results.ma)

#########################################################################
# RUNNING NSCA STEP: Looking for subclasses within all samples compared #
#########################################################################

## Not run as example
# deco.results.ma.uns <- decoNSCA(sub = sub.ma.3r.1K.uns, v = 80, method = "ward.D",
#                         k.control = 3, k.case = 3, samp.perc = 0.05, rep.thr = 3)

########################################
# RUNNING NSCA STEP: Multiclass design #
########################################

## Not run as example
# deco.results.ma.multi <- decoNSCA(sub = sub.ma.3r.1K.multi, v = 80, method = "ward.D",
#                         k.control = 3, k.case = 3, samp.perc = 0.05, rep.thr = 3)

}
