############################################ 
#' Class deco
#'
#' Class \code{deco}.
#'
#' @name deco-class
#' @rdname deco-class
#' @exportClass deco
setClass(Class = "deco", slots = c(featureTable = "data.frame", NSCAcluster = "list", 
    incidenceMatrix = "data.frame", classes = "factor", pos.iter = "numeric", 
    control = "character", q.val = "numeric", rep.thr = "numeric", samp.perc = "numeric", 
    subsampling.call = "call", nsca.call = "call"), sealed = FALSE)

############################################ 
#' Method featureTable
#' @name featureTable
#' @rdname deco-class
#' @exportMethod featureTable
setGeneric("featureTable", function(object) standardGeneric("featureTable"))

#'
#' @rdname deco-class
#' @aliases featureTable,deco-method
setMethod("featureTable", signature = "deco", definition = function(object) object@featureTable)

############################################ 
#' Method NSCAcluster
#' @name NSCAcluster
#' @rdname deco-class
#' @exportMethod NSCAcluster
setGeneric("NSCAcluster", function(object) standardGeneric("NSCAcluster"))

#'
#' @rdname deco-class
#' @aliases NSCAcluster,deco-method
setMethod("NSCAcluster", signature = "deco", definition = function(object) object@NSCAcluster)

############################################ 
#'
#' @rdname deco-class
#' @aliases show,deco-method
setMethod("show", "deco", function(object) {
    cat("\nDecomposing Heterogeneous Cohorts from Omic profiling: DECO\nSummary:\n")
    cat("\nAnalysis design: ")
    if (length(object@NSCAcluster) == 2) {
        cat("Binary\nClasses compared:")
        print(table(object@classes))
    } else if (all(is.na(object@classes))) {
        cat("Unsupervised\n")
    } else {
        cat("Multiclass\nClasses compared:")
        print(table(object@classes))
    }
    
    cat("\n")
    thr <- data.frame(`RDA q.value` = object@q.val, `Minimum repeats` = object@rep.thr, 
        `Percentage of affected samples` = object@samp.perc * 100, `NSCA variability` = object@NSCAcluster[[1]]$var)
    rownames(thr) <- "Thresholds"
    printCoefmat(thr, digits = 3)
    cat("\nNumber of features out of thresholds:", dim(object@featureTable)[1], 
        "\n")
    cat("Number of samples affected:", dim(object@incidenceMatrix)[2], "\n")
    cat("Number of positive RDA comparisons:", object@pos.iter, "\n")
    cat("Number of total RDA comparisons:", round(object@featureTable[1, c("Repeats")]/object@featureTable[1, 
        c("FR.Repeats")], digits = 0), "\n\n")
    cat("RDA call:\n")
    print(object@subsampling.call)
    cat("NSCA call:\n")
    print(object@nsca.call)
})

############################################ 
#' Method summary

#'
#' @rdname deco-class
#' @aliases summary,deco-method
setMethod("summary", "deco", function(object, ...) {
    cat("Decomposing Heterogeneous Cohorts from Omic profiling: DECO\nSummary:\n")
    cat("\nAnalysis design: ")
    if (length(object@NSCAcluster) == 2) {
        cat("Binary\nClasses compared:")
        print(table(object@classes))
    } else if (all(is.na(object@classes))) {
        cat("Unsupervised\n")
    } else {
        cat("Multiclass\nClasses compared:")
        print(table(object@classes))
    }
    
    cat("\n")
    thr <- data.frame(`RDA q.value` = object@q.val, `Minimum repeats` = object@rep.thr, 
        `Percentage of affected samples` = object@samp.perc * 100, `NSCA variability` = object@NSCAcluster[[1]]$var)
    rownames(thr) <- "Thresholds"
    printCoefmat(thr, digits = 3)
    cat("\nNumber of features out of thresholds:", dim(object@featureTable)[1], 
        "\n")
    if ("Profile" %in% colnames(object@featureTable)) {
        cat("Feature profile table:")
        print(table(object@featureTable$Profile))
    }
    cat("Number of samples affected:", dim(object@incidenceMatrix)[2], "\n")
    cat("Number of positive RDA comparisons:", object@pos.iter, "\n")
    cat("Number of total RDA comparisons:", round(object@featureTable[1, c("Repeats")]/object@featureTable[1, 
        c("FR.Repeats")], digits = 0), "\n\n")
})
