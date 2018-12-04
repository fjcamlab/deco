########################################################################################
### DECO ###
# algorithm to find differential significant markers in heterogeneous cohorts
# including a subsampling approach and posterior integration with a
# Non-Symmetrical Correspondence Analysis to create the new heterogeneity statistic
# (h-statistic).
#
# Authors: Francisco J. Campos-Laborie, Jose Manuel Sanchez-Santos & Javier De Las Rivas
# Bioinformatics & Functional Genomics Group
# Cancer Research Center (CiC-IBMCC, CSIC/USAL)
# Salamanca, Spain
# http://www.cicancer.org/
# http://bioinfow.dep.usal.es
#
########################################################################################

###################################
###         decoRDA           ###
###################################
## Iterative function to discover differentially expressed features
## using LIMMA R package (eBayes method) along the data.

decoRDA <- function(data, classes = NA, control = NA, r = NULL,
                    q.val = 0.01, iterations = 10000, cpus = 2, parallel = FALSE,
                    annot = FALSE, id.type = NA, attributes = NA,
                    rm.xy = FALSE, pack.db = "org.Hs.eg.db") {
    call <- match.call()

    # Assessing 'data' input.
    msg <- NA
    if (is(data, "SummarizedExperiment") | is(data, "MultiAssayExperiment")) {
        if (!exists("SummarizedExperiment")) {
            msg <- c("ERROR: 'SummarizedExperiment' is not in library.")
        } else {
            requireNamespace("SummarizedExperiment",
                             quietly = TRUE)
        }
        if (is(assay(data), "matrix")) {
            data <- assay(data)
        } else {
            msg <- c("ERROR: 'data' object must include only one 'assay'.")
        }
    }
    if (!is.matrix(data)) {
        msg <- c("Data must be a matrix or SummarizedExperiment object")
    }
    data <- as.matrix(data)


    if (!is.na(msg)) {
        stop(msg)
    }

    # Setting up temporary directory to write intermediate results.
    temp.path <- tempdir()
    message("'temp' folder will be created to store internal loop data.")

    # Annotating IDs to chromosome location (locus) in order to remove those placed on X or Y.
    # It is proposed to avoid significant results from unbalanced gender contrasts.
    if (rm.xy & annot) {
        if (!exists(pack.db)) {
            stop(
                "ERROR: Annotation package '",
                pack.db, "' has not been loaded."
            )
        } else {
            requireNamespace(pack.db, quietly = TRUE)
        }
        if (id.type %in% columns(x = get(pack.db))) {
            infogenes <- AnnotateDECO(
                ids = rownames(data),
                id.type = id.type, attributes = c("CHR", "CHR"),
                pack.db = pack.db, verbose = FALSE
            )
            chrTest <- infogenes[, "CHR"] %in% c("X", "Y")
            if (length(which(chrTest)) > 1) {
                data <- data[rownames(infogenes)[!chrTest], ]
                msg <- paste(
                    .timestamp(), "-- Features located in X or Y chromosome have been filtered:\n",
                    length(which(chrTest)), "features have been discarded."
                )
            } else {
                msg <- paste(.timestamp(), " All features are placed on X or Y chromosome,
                       filter will not be applied.")
            }
        } else {
            msg <- paste(.timestamp(), "It was not possible to assign chromosome location
                   to input IDs.")
        }
        message(msg)
    }

    ## Initial matrix for resampling design...
    results <- matrix(0, nrow = 1, ncol = 6)

    ## Defining classes and sample sizes
    unsup <- all(is.na(classes))
    multi <- FALSE
    if (!unsup) {
        combin <- .combCalcM(
            data, classes, control,
            iterations, multi, r, results
        )

        categories.control <- combin$categories.control
        categories.case <- combin$categories.case
        multi <- combin$multi
        cl1 <- combin$cl1
        cl2 <- combin$cl2
        n1 <- combin$n1
        n2 <- combin$n2
        results <- combin$results
        data <- combin$data
    } else {
        n1 <- dim(data)[2]
        n2 <- dim(data)[2]
        categories.control <- seq_len(dim(data)[2])
        categories.case <- seq_len(dim(data)[2])
        message("\n Classes vector not defined as input.
                UNSUPERVISED analysis will be carry out.\n")
    }

    ## Setting up subsampling size 'r'.
    if ((is.null(r) || r <= 0)) {
        r <- round(sqrt(min(c(n1, n2))), digits = 0)
    }
    if (r < 3) {
        r <- 3
        message(" Resampling size set to 3.\n")
    }
    if (r > n1 || r > n2) {
        stop("ERROR: Resampling size can not be
             higher than samples.")
    }

    ## Message to user
    if (unsup) {
        msg <- paste(
            .timestamp(), "-- Resampling design:\n Unsupervised analysis for",
            dim(data)[2], "samples.\n Resampling size:", r,
            "\n adj.p.value threshold:", q.val
        )
    } else if (!multi) {
        msg <- paste(
            .timestamp(), "-- Resampling design:\n Binary analysis for",
            dim(data)[2], "samples.\n Control or 0 --> '", cl1, "' with", n1,
            "samples.\n Case or 1    --> '", cl2, "' with",
            n2, "samples.\n", "Resampling size for both classes:", r, "\n adj.p.value threshold:", q.val
        )
    } else {
        msg <- paste(
            .timestamp(), "-- Resampling design:\n Multiple supervised analysis for",
            dim(data)[2], "samples and", length(levels(classes)), "classes:\n",
            paste(levels(classes), collapse = "/"),
            ".\n Resampling size for all classes:", r, "\n adj.p.value threshold:", q.val
        )
    }

    message(msg)

    ## Combinatorial calculations for binary and unsupervised designs
    results <- .combCalc(
        iterations, multi, classes, r, results,
        categories.control, categories.case,
        n1, n2
    )

    ## SUBSAMPLING STEP using LIMMA
    limmaRes <- .LIMMAcalc(
        data, results, classes, control, q.val,
        call, r, cpus, parallel, temp.path,
        multi, n1, n2
    )
    res <- limmaRes$res
    results <- limmaRes$results

    ## FREQUENCY MATRIX:
    # Reading intermediate files to produce a frequency matrix
    message(.timestamp(), " -- Summarizing results...")
    intermediate <- .freqMatrix(
      limmaRes, data,
      n1, ifelse(unsup, 0, n2),
      r, multi, unsup
      )

    # Matrix of combinations and number of events.
    ord <- order(results[, "nFeatures"], decreasing = TRUE)
    res$results <- results[ord, -(r * 2 + 2)]
    res$pos.iter <- length(limmaRes$limma1) / 2

    if (!unsup & !multi) {
        res$incidenceMatrix <- interleave(intermediate$UP, intermediate$DOWN,
            append.source = TRUE,
            sep = "deco", drop = FALSE
        )
    } else {
        res$incidenceMatrix <- intermediate$MULTI
    }

    if (length(intermediate$both.feat) > 0) {
        res$incidenceMatrix <- rbind(res$incidenceMatrix,
                                     intermediate$both.feat.mx)
    }

    res$incidenceMatrix <- res$incidenceMatrix[, colnames(
        res$incidenceMatrix
    ) %in% colnames(data)]

    ### CALCULATING FINAL STATISTICS
    ## Creating a table to contain all statistical information about features.
    j <- length(limmaRes$limma1) / 2
    top_eje <- intermediate$top_eje[2:length(rownames(intermediate$top_eje)), ]
    res$subStatFeature <- as.data.frame(matrix(ncol = 12,
                                               nrow = length(unique(top_eje[, c("ID")]))))
    colnames(res$subStatFeature) <- c(
        "ID", "UpDw", "Avrg.logFC", "Best.adj.P.Val",
        "Repeats", "FR.Repeats", "RF.Positive.Repeats",
        "Chi.Square", "P.Val.ChiSq",
        "ChiSq.adj.P.Val.FDR"
    )
    res$subStatFeature[, c("ID")] <- unique(top_eje[, c("ID")])
    res$subStatFeature[, 3:10] <- 0
    ncomb <- dim(results)[1]

    ## Summarizing statistics...
    if (dim(res$subStatFeature)[1] > 1) {
        cpusC <- cpus
    } else {
        cpusC <- 1
    }
    suppressWarnings(limma2 <- bplapply(seq_len(dim(res$subStatFeature)[1]),
        FUN = .statCalc,
        BPPARAM = MulticoreParam(cpusC),
        tab = res$subStatFeature,
        top_eje, ncomb, j
    ))

    res$subStatFeature <- as.data.frame(do.call(rbind, limma2))
    res$subStatFeature[, 3:10] <- suppressWarnings(
        apply(res$subStatFeature[, 3:10], 2, as.numeric)
    )

    message(
        .timestamp(), " -- Computed ",
        length(unique(top_eje[, c("ID")])), " features.\n"
    )

    ## Adjusting p.values from Chi.Square.
    res$subStatFeature[, c("ChiSq.adj.P.Val.FDR")] <- p.adjust(
        res$subStatFeature[, c("P.Val.ChiSq")],
        method = "fdr",
        n = length(res$subStatFeature[, c("P.Val.ChiSq")])
    )

    # Ordering columns...
    res$subStatFeature <- res$subStatFeature[, c(
        "ID", "UpDw", "Avrg.logFC", "Best.adj.P.Val", "Repeats", "FR.Repeats",
        "RF.Positive.Repeats", "Chi.Square", "P.Val.ChiSq", "ChiSq.adj.P.Val.FDR"
    )]

    ## Annotating features using 'AnnotateDECO' function.
    if (annot) {
        infogenes <- AnnotateDECO(
            ids = as.character(res$subStatFeature[, c("ID")]),
            id.type = id.type, attributes = attributes,
            pack.db = pack.db
        )
        res$subStatFeature <- cbind(
            res$subStatFeature,
            infogenes[as.vector(res$subStatFeature[, c("ID")]), ]
        )
    }
    # Making up all statistics.
    res$subStatFeature <- droplevels.data.frame(as.data.frame(res$subStatFeature))
    res$subStatFeature[, 3:10] <- suppressWarnings(
        apply(res$subStatFeature[, 3:10], 2, as.numeric)
    )
    if (!unsup & !multi) {
        res$subStatFeature[, c("UpDw")] <- intermediate$UpDw[match(
            as.character(res$subStatFeature[, c("ID")]), names(intermediate$UpDw)
        )]
    } else {
        res$subStatFeature <- res$subStatFeature[, colnames(res$subStatFeature) != "UpDw"]
    }

    # Standard.Chi.Square
    res$subStatFeature[, "Standard.Chi.Square"] <- ((as.numeric(res$subStatFeature[, c("Chi.Square")]) -
        (2 * as.numeric(res$subStatFeature[, c("Repeats")]))) /
        sqrt(4 * as.numeric(res$subStatFeature[, c("Repeats")])))

    res$subStatFeature <- res$subStatFeature[order(
        res$subStatFeature$Standard.Chi.Square,
        decreasing = TRUE
    ), ]
    rownames(res$subStatFeature) <-
      res$subStatFeature[, "ID"] <-
        as.character(res$subStatFeature[, "ID"])

    # Removing temporary dir and all intermediate files generated.
    unlink(temp.path, recursive = TRUE)

    message(.timestamp(), " -- Done.\n")

    return(res)
}


#################################
###        decoNSCA           ###
#################################
## Function to carry out factorial analysis with NSCA (Non-Symmetrical Correspondence Analysis) of
## the incidenceMatrix generated by 'DECO' subsampling function.

decoNSCA <- function(sub, v = 80, k.control = NULL,
                     k.case = NULL, rep.thr = 3,
                     samp.perc = 0.05,
                     method = "ward.D",
                     parallel = FALSE, cpus = 2) {
    call <- match.call()

    # Agglomeration method for hierarchical clustering of samples.
    if (all(!(is.null(method))) &&
        !(method %in% c(
            "ward.D", "ward.D2", "single", "complete",
            "average", "mcquitty", "median", "centroid"
        ))) {
        stop("ERROR: Input 'method.heatmap' have to be one of 'hclust'{stats} function methods:
         ward.D, ward.D2, single, complete, average, mcquitty, median or centroid")
    }

    data <- sub$data

    # Variable with original IDs used to substitute modified IDs of incidenceMatrix.
    message(.timestamp(), " -- Applying repeats threshold...\n")

    ## Threshold for Repeats
    Sub <- .repThr(
        sub, rep.thr,
        samp.perc, cpus
    )
    g.names <- Sub$g.names
    sub <- Sub$sub

    # Setting up both classes and control for 'binary' design.
    if (all(!(is.na(sub$classes))) & length(levels(sub$classes)) <= 2) {
        # Instructions just for 'binary' RDA design.
        if (is.na(sub$control)) {
            cl1 <- levels(sub$classes)[1]
            cl2 <- levels(sub$classes)[2]
            ind1 <- which(sub$classes == cl1)
            ind2 <- which(sub$classes == cl2)
        } else {
            cl1 <- sub$control
            cl2 <- levels(sub$classes)[which(levels(sub$classes) != sub$control)]
            ind1 <- which(sub$classes == cl1)
            ind2 <- which(sub$classes == cl2)
            sub$classes <- sub$classes[c(ind1, ind2)]
        }
        n1 <- length(ind1)
        n2 <- length(ind2)
        sub$incidenceMatrix <- sub$incidenceMatrix[, names(sub$classes)]

        # Calculating raw 'delta.signal' differences between classes for all differential features.
        d <- sub$data[rownames(sub$data) %in% rownames(
            sub$subStatFeature
        ), ]

        delta.signal <- rowMeans(d[, (n1 + 1):(n1 + n2)]) -
            rowMeans(d[, seq_len(n1)])

        sub$subStatFeature <- data.frame(sub$subStatFeature[order(sub$subStatFeature[, c("ID")]), ],
            delta.signal = delta.signal[order(names(delta.signal))]
        )
        # Standard deviation per class
        sd.Ctrl <- apply(sub$data[rownames(sub$subStatFeature), seq_len(n1)], 1, sd)
        sd.Case <- apply(sub$data[rownames(sub$subStatFeature), (n1 + 1):(n1 + n2)], 1, sd)

        # Classification of feature profiles: 'ideal', 'generic', 'specific, and 'both'.
        overlap <- overlapC(
          sub, cl1,
          ifelse(parallel, cpus, 1)
        )

        profile <- overlap$profile
        overlap <- overlap$overlap

        UpDw <- sub$subStatFeature$UpDw
        names(UpDw) <- as.character(sub$subStatFeature$ID)
    } else {
        n1 <- dim(data)[2]
        n2 <- dim(data)[2]
        SD <- apply(sub$data, 1,
                function(x) sd(x, na.rm = TRUE))
    }
    # Applying NSCA to all samples or both classes.
    z <- 1
    while (z %in% seq_len(2)) {
        # NSCA for 'unsupervised' and 'multiclass' design.
        if (all(is.na(sub$classes)) || length(levels(sub$classes)) > 2) {
            # Removing features or samples with any differential event.
            mx <- sub$incidenceMatrix[rowSums(sub$incidenceMatrix) > 0, ]
            mx <- mx[, colSums(sub$incidenceMatrix) > 0]
            message(.timestamp(), "-- Calculating all samples together...")
            # NSCA
            nsc.res <- NSCACluster(
                mx = mx, data = sub$data, k = k.control,
                method.dend = method, id.names = g.names, v = v
            )
            k <- k.control
            diffMX <- cbind(sub$subStatFeature, nsc.res$infoFeature[rownames(sub$subStatFeature), ])
            diffMX <- cbind(diffMX, sd = SD[rownames(sub$subStatFeature)])
            # Order for columns of final table with statistical information.
            coln <- c(
                "ID", "SYMBOL", "hgnc_symbol", "Repeats", "Repeats.index", "FR.Repeats", "Avrg.logFC",
                "Standard.Chi.Square", "P.Val.ChiSq", "ChiSq.adj.P.Val.FDR",
                "sd", "Tau.feature", "Dendrogram.group",
                "h.Best", "h.Range", "GENENAME"
            )
            ord <- na.omit(match(c(coln, colnames(diffMX)[!colnames(diffMX) %in% coln]), colnames(diffMX)))
            z <- 3
        }
        # NSCA for 'control' samples
        if (z == 1) {
            # Removing features or samples with any differential event.
            mx <- sub$incidenceMatrix[rowSums(sub$incidenceMatrix[, seq_len(n1)]) > 0, seq_len(n1)]
            mx <- mx[, colSums(sub$incidenceMatrix[, seq_len(n1)]) > 0]
            message(.timestamp(), "-- Calculating control group...")
            # NSCA
            nsc.res1 <- NSCACluster(
                mx = mx, data = sub$data, k = k.control,
                method.dend = method, UpDw = UpDw, dir = "UP",
                id.names = g.names, v = v, label = cl1,
                raw = rowMeans(sub$data[, colnames(mx)])
            )
            colnames(nsc.res1$infoFeature) <- paste(colnames(nsc.res1$infoFeature), "Ctrl", sep = ".")
            diffMX <- cbind(sub$subStatFeature,
                sd.Ctrl = sd.Ctrl[rownames(sub$subStatFeature)],
                nsc.res1$infoFeature[rownames(sub$subStatFeature), c(
                    "Tau.feature.Ctrl",
                    "Dendrogram.group.Ctrl", "h.Best.Ctrl", "h.Range.Ctrl"
                )]
            )
        }
        # NSCA for 'case' samples
        if (z == 2) {
            # Removing features or samples with any differential event.
            mx <- sub$incidenceMatrix[rowSums(sub$incidenceMatrix[, (n1 + 1):(n1 + n2)]) > 0, (n1 + 1):(n1 + n2)]
            mx <- mx[, colSums(sub$incidenceMatrix[, (n1 + 1):(n1 + n2)]) > 0]
            message(.timestamp(), " -- Calculating case group...")
            # NSCA
            nsc.res2 <- NSCACluster(
                mx = mx, data = sub$data, k = k.case, method.dend = method, UpDw = UpDw, dir = "DOWN",
                id.names = g.names, v = v, label = cl2, raw = rowMeans(sub$data[, colnames(mx)])
            )
            colnames(nsc.res2$infoFeature)[seq_len(length(colnames(nsc.res2$infoFeature)) - 1)] <-
                paste(colnames(nsc.res2$infoFeature)[seq_len(length(colnames(nsc.res2$infoFeature)) - 1)], "Case", sep = ".")
            diffMX <- cbind(diffMX,
                sd.Case = sd.Case[rownames(sub$subStatFeature)],
                overlap.Ctrl.Case = overlap[rownames(sub$subStatFeature)],
                nsc.res2$infoFeature[rownames(sub$subStatFeature), c(
                    "Tau.feature.Case", "Dendrogram.group.Case",
                    "h.Best.Case", "h.Range.Case"
                )]
            )
            # Standard.Chi.Square calculation.
            diffMX <- cbind(diffMX, Profile = profile[rownames(sub$subStatFeature)])
            # Order for columns of final table with statistical information.
            coln <- c(
                "ID", "SYMBOL", "hgnc_symbol", "UpDw", "Profile", "overlap.Ctrl.Case",
                "Repeats", "Repeats.index", "FR.Repeats",
                "delta.signal", "Avrg.logFC", "Standard.Chi.Square", "P.Val.ChiSq", "ChiSq.adj.P.Val.FDR",
                "sd.Ctrl", "Tau.feature.Ctrl", "Dendrogram.group.Ctrl", "h.Best.Ctrl",
                "h.Range.Ctrl", "sd.Case", "Tau.feature.Case", "Dendrogram.group.Case",
                "h.Best.Case", "h.Range.Case", "GENENAME"
            )
            ord <- na.omit(match(c(coln, colnames(diffMX)[!colnames(diffMX) %in% coln]), colnames(diffMX)))
        }
        z <- z + 1
    }
    # Making up of final statistical table.
    diffg <- data.frame(diffMX[, c(ord)], diffMX[, na.action(ord)])
    colnames(diffg) <- c(colnames(diffMX)[c(ord)], colnames(diffMX)[na.action(ord)])
    # Creating final result object.
    if (!(all(is.na(sub$classes))) & !length(levels(sub$classes)) > 2) {
        ord <- with(diffg, order(Profile, -Standard.Chi.Square))
        diffg <- diffg[ord, ]
        results <- new("deco",
            featureTable = as.data.frame(diffg[, !duplicated(colnames(diffg))]),
            NSCAcluster = list(Control = nsc.res1, Case = nsc.res2),
            incidenceMatrix = sub$incidenceMatrix,
            classes = as.factor(sub$classes), pos.iter = sub$pos.iter,
            control = as.character(sub$control), q.val = sub$q.val,
            rep.thr = rep.thr, samp.perc = samp.perc,
            subsampling.call = sub$call, nsca.call = call
        )
    } else {
        ord <- order(rowMeans(apply(cbind(
            diffg$Standard.Chi.Square, diffg$Standard.Chi.Square, diffg$Repeats,
            diffg$sd, diffg$h.Range
        ), 2, function(x) rank(-x))))
        diffg <- diffg[ord, ]
        results <- new("deco",
            featureTable = as.data.frame(diffg[, !duplicated(colnames(diffg))]),
            NSCAcluster = list(All = nsc.res), incidenceMatrix = sub$incidenceMatrix,
            classes = as.factor(sub$classes), pos.iter = sub$pos.iter,
            control = as.character(sub$control), q.val = sub$q.val,
            rep.thr = rep.thr, samp.perc = samp.perc,
            subsampling.call = sub$call, nsca.call = call
        )
    }

    return(results)
}
