############################################
setClass(
    Class = "deco",
    slots = c(
        featureTable = "data.frame", NSCAcluster = "list",
        incidenceMatrix = "data.frame", classes = "factor",
        pos.iter = "numeric", control = "character",
        q.val = "numeric", rep.thr = "numeric", samp.perc = "numeric",
        subsampling.call = "call", nsca.call = "call"
    ),
    sealed = FALSE
)

############################################
setGeneric("featureTable", function(object) standardGeneric("featureTable"))

setMethod("featureTable",
    signature = "deco", definition =
        function(object) object@featureTable
)

###
setGeneric("NSCAcluster", function(object) standardGeneric("NSCAcluster"))

setMethod("NSCAcluster",
    signature = "deco", definition =
        function(object) object@NSCAcluster
)

################################################
################################################
# R show function for 'deco' object

setMethod("show", "deco", function(object) {
    cat("\nDecomposing Heterogeneous Cohorts from Omic profiling: DECO\nSummary:\n")
    cat("\nAnalysis design: ")
    if (length(object@NSCAcluster) == 2) {
        cat("Binary\nClasses compared:")
        print(table(object@classes))
    }
    else if (all(is.na(object@classes))) {
        cat("Unsupervised\n")
    }
    else {
        cat("Multiclass\nClasses compared:")
        print(table(object@classes))
    }

    cat("\n")
    thr <- data.frame(
        "RDA q.value" = object@q.val, "Minimum repeats" = object@rep.thr,
        "Percentage of affected samples" = object@samp.perc * 100,
        "NSCA variability" = object@NSCAcluster[[1]]$var
    )
    rownames(thr) <- "Thresholds"
    printCoefmat(thr, digits = 3)
    cat("\nNumber of features out of thresholds:", dim(object@featureTable)[1], "\n")
    cat("Number of samples affected:", dim(object@incidenceMatrix)[2], "\n")
    cat("Number of positive RDA comparisons:", object@pos.iter, "\n")
    cat(
        "Number of total RDA comparisons:", round(
            object@featureTable[1, c("Repeats")] / object@featureTable[1, c("FR.Repeats")],
            digits = 0
        ),
        "\n\n"
    )
    cat("RDA call:\n")
    print(object@subsampling.call)
    cat("NSCA call:\n")
    print(object@nsca.call)
})

################################################
################################################
# R summary function for 'deco' object

setMethod("summary", "deco", function(object, ...) {
    cat("Decomposing Heterogeneous Cohorts from Omic profiling: DECO\nSummary:\n")
    cat("\nAnalysis design: ")
    if (length(object@NSCAcluster) == 2) {
        cat("Binary\nClasses compared:")
        print(table(object@classes))
    }
    else if (all(is.na(object@classes))) {
        cat("Unsupervised\n")
    }
    else {
        cat("Multiclass\nClasses compared:")
        print(table(object@classes))
    }

    cat("\n")
    thr <- data.frame(
        "RDA q.value" = object@q.val, "Minimum repeats" = object@rep.thr,
        "Percentage of affected samples" = object@samp.perc * 100,
        "NSCA variability" = object@NSCAcluster[[1]]$var
    )
    rownames(thr) <- "Thresholds"
    printCoefmat(thr, digits = 3)
    cat("\nNumber of features out of thresholds:", dim(object@featureTable)[1], "\n")
    if ("Profile" %in% colnames(object@featureTable)) {
        cat("Feature profile table:")
        print(table(object@featureTable$Profile))
    }
    cat("Number of samples affected:", dim(object@incidenceMatrix)[2], "\n")
    cat("Number of positive RDA comparisons:", object@pos.iter, "\n")
    cat(
        "Number of total RDA comparisons:", round(
            object@featureTable[1, c("Repeats")] / object@featureTable[1, c("FR.Repeats")],
            digits = 0
        ),
        "\n\n"
    )
})
