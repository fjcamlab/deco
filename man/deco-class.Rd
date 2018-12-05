\name{deco-class}
\Rdversion{1.1}
\docType{class}
\alias{deco-class}
\alias{summary,deco-method}
\alias{show,deco-method}
\alias{featureTable,deco-method}
\alias{NSCAcluster,deco-method}

\title{Class \code{"deco"}}
\description{
A list R object returned by \code{decoNSCA} function from DECO R package.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("deco", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{featureTable}:}{Object of class \code{"data.frame"}: table of results with feature statistical information about NSCA merged with RDA information contained in 'subStatFeature' table from \emph{decoRDA} function results.
 }
    \item{\code{NSCAcluster}:}{Object of class \code{"list"}: object from \emph{NSCAcluster} intern R function containing information about NSCA. If 'Binary' analysis was previously set, two lists corresponding to each class are contained (Control and Case). }
    \item{\code{incidenceMatrix}:}{Object of class \code{"data.frame"}: absolute frequency matrix of 'm' features by 'n' samples size that summarizes differential events per feature per sample. It is essential for Non-symmetrical correspondence analysis.}
    \item{\code{classes}:}{Object of class \code{"factor"}: vector with class labels.}
    \item{\code{pos.iter}:}{Object of class \code{"numeric"}: number of iterations with DE signal.}
    \item{\code{control}:}{Object of class \code{"character"}: control label.}
    \item{\code{q.val}:}{Object of class \code{"numeric"}: adjusted.p.value threshold to LIMMA.}
    \item{\code{rep.thr}:}{Object of class \code{"numeric"}: number of minimum repeats per sample to establish a differential event threshold.}
    \item{\code{samp.perc}:}{Object of class \code{"numeric"}: minimum amount of samples showing \code{rep.thr} number of differential events to consider one feature as relevant (non-noisy) after RDA or subsampling step.}
    \item{\code{subsampling.call}:}{Object of class \code{"numeric"}: number of minimum repeats per sample to establish a differential event threshold.}
    \item{\code{nsca.call}:}{Object of class \code{"numeric"}: number of minimum repeats per sample to establish a differential event threshold.}
  }
}
\section{Methods}{ 
  \describe{ 
    \item{summary}{\code{summary(object = "deco")}: Prints a summary of the 
      information inside a 'deco' R object.} 
    }
  \describe{ 
    \item{show}{\code{show(object = "deco")}: Prints the 
      information inside a 'deco' R object.} 
    }
  \describe{ 
    \item{featureTable}{\code{featureTable(object = "deco")}: Extract the 
      featureTable containing all feature statistics.} 
    }
  \describe{ 
    \item{NSCAcluster}{\code{NSCAcluster(object = "deco")}: Extract the 
      NSCAcluster list containing all information and tables derived from decoNSCA step.} 
    }
} 
\author{
Francisco Jose Campos Laborie. <fjcamlab@usal.es>
}
\examples{
showClass("deco")

showMethods("deco")
}
\keyword{classes}
