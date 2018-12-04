###########################################################################################################
########## Intern functions ##########

# Printing messages

.timestamp <- function()
    format(Sys.time(), " %H:%M:%S")

# Temporary directory

.createTempFile <- function()
{
  uniqueName <- unlist(
    strsplit(tempfile(pattern = "DECO_temp"), 
                                split = "//")
    )[2]
  
  uniqueName <- paste(
    unlist(strsplit(uniqueName, split = "/")), 
    collapse = "_"
  )
  
  temp.path <- paste(getwd(), "/", uniqueName, sep = "")
  
  dir.create(temp.path)
  
  return(temp.path)
}

################################################
################################################
# localMinima and localMaxima

localMinima <- function(x) {
    y <- diff(c(.Machine$integer.max, x)) > 0L
    rle(y)$lengths
    y <- cumsum(rle(y)$lengths)
    y <- y[seq.int(1L, length(y), 2L)]
    if (x[[1]] == x[[2]]) {
        y <- y[-1]
    }

    return(y)
}

localMaxima <- function(x) {
    y <- diff(c(-.Machine$integer.max, x)) > 0L
    rle(y)$lengths
    y <- cumsum(rle(y)$lengths)
    y <- y[seq.int(1L, length(y), 2L)]
    if (x[[1]] == x[[2]]) {
        y <- y[-1]
    }

    return(y)
}

################################################
################################################
# Calculate quantile threshold per standard deviation

quantSD <- function(sdx, p = NULL) {
    if (is.null(p)) {
        p <- 0.75
    }
    if (quantile(sdx, probs = p) <= median(range(sdx))) {
        thr <- quantile(sdx, probs = p)
    } else {
        thr <- median(sdx)
    }
    return(thr)
}

################################################
################################################
# AnnotateDECO

AnnotateDECO <- function(ids, id.type, attributes = NA,
                         pack.db = "org.Hs.eg.db",
                         rm.redundant = TRUE, verbose = FALSE) {
    if (!exists(pack.db)) {
        stop("ERROR: Annotation package '", pack.db, "' is not in library.")
    } else {
        requireNamespace(pack.db, quietly = TRUE)
    }

    if (verbose) {
        message(.timestamp(), "-- Annotating IDs with Bioconductor...")
    }

    if (is.na(attributes)) {
        attributes <- c(
            "SYMBOL", "ENSEMBL", "ENTREZID",
            "CHR", "CHRLOC", "GENENAME"
        )
    }

    if (id.type == "PROBEID") {
        xx <- as.list(get(paste(unlist(strsplit(pack.db, split = ".db")),
            "ENSEMBL",
            sep = ""
        )))

        infogenes <- suppressMessages(AnnotationDbi::select(
            x = get(pack.db), obj = names(xx),
            columns = c(attributes, id.type),
            keys = keys(get(pack.db), keytype = "ENSEMBL"),
            keytype = "ENSEMBL"
        ))
    } else {
        infogenes <- suppressMessages(AnnotationDbi::select(
            x = get(pack.db), obj = ids,
            columns = attributes,
            keys = keys(get(pack.db), keytype = id.type),
            keytype = id.type
        ))
    }

    infogenes <- infogenes[apply(
        infogenes[, attributes], 1,
        function(x) all(!is.na(x))
    ), ]

    if (all(is.na(infogenes))) {
        stop("Annotation library does not found input IDs.")
    }

    infogenes <- infogenes[order(infogenes[, id.type]), ]
    infogenes <- infogenes[infogenes[, id.type] %in% ids, ]

    if (rm.redundant) {
        infogenes <- infogenes[!(duplicated(infogenes[, id.type])), ]
        rownames(infogenes) <- infogenes[, id.type]

        if (verbose) {
            message(.timestamp(), "-- Removed duplicated info from IDs.")
        }
    }
    return(infogenes)
}

##############################
##############################
## innerProductAssign

innerProductAssign <- function(inner, samples = NA,
                               control = NA, analysis) {
    samples <- samples[names(inner)]
    if (analysis == "Binary") {
        res <- 0
        inner0 <- inner
        for (j in seq_len(length(table(samples))))
        {
            if (j == 1) {
                inner <- inner0[names(samples[samples == control])]
            } else {
                inner <- inner0[names(samples[samples != control])]
            }

            error <- try(
                expr = {
                    d <- density.lf(x = as.matrix(inner))
                },
                silent = TRUE
            )
            if (inherits(error, "try-error")) {
                d <- density(x = inner)
            }
            x <- d$x[localMaxima(d$y)]

            error <- try(
                expr = {
                    p <- kmeans(inner, centers = x)
                },
                silent = TRUE
            )
            if (inherits(error, "try-error")) {
                error <- try(
                    expr = {
                        p <- kmeans(inner, centers = length(x) - 1, nstart = 10)
                    },
                    silent = TRUE
                )
            }

            res <- c(res, p$cluster + max(res))
        }
        res <- res[2:length(res)]
        res <- res[order(res, inner0)]
    } else {
        error <- try(
            expr = {
                d <- density.lf(x = as.matrix(inner))
            },
            silent = TRUE
        )
        if (inherits(error, "try-error")) {
            d <- density(x = inner)
        }
        x <- d$x[localMaxima(d$y)]

        error <- try(
            expr = {
                p <- kmeans(inner, centers = x)
            },
            silent = TRUE
        )
        if (inherits(error, "try-error")) {
            error <- try(
                expr = {
                    p <- kmeans(inner, centers = length(x), nstart = 10)
                },
                silent = TRUE
            )
        }
        res <- p$cluster[order(p$cluster, inner)]
    }
    return(list(cl = res, centroid = x))
}

##############################
##############################
## Distance function with correlation

distFunc <- function(x, use = "pairwise.complete.obs",
                     cor.method = "pearson") {
    co.x <- cor(x, use = use, method = cor.method)
    dist.co.x <- 1 - co.x
    return(as.dist(dist.co.x))
}


##############################
##############################
## Cophenetic correlation and clustering dendogram

cophDECO <- function(data, method.heatmap = "ward.D",
                     k = NULL, scale = FALSE,
                     verbose = TRUE, coph = TRUE) {
    if (is.null(method.heatmap)) {
        method.heatmap <- c(
            "ward.D", "ward.D2", "single",
            "complete", "average", "mcquitty",
            "median", "centroid"
        )
    }

    if (scale) {
        data <- scale(data)
    }

    d <- distFunc(x = as.matrix(data), cor.method = "pearson")
    d[is.na(d)] <- 0

    if (coph) {
        hmet <- vapply(method.heatmap, function(x) {
            sampleTree <- as.dendrogram(hclust(d, method = x))
            coph <- cor(c(d), c(cophenetic(sampleTree)),
                method = "pearson"
            )
        }, FUN.VALUE = numeric(1))
    } else {
        hmet <- rep(1, length(method.heatmap))
        names(hmet) <- method.heatmap
    }

    sampleTree <- hclust(d, method = names(hmet[order(hmet, decreasing = TRUE)])[1])
    sampleTree$height <- sampleTree$height[order(sampleTree$height)]

    if (is.null(k)) {
        p <- c()
        for (i in 2:(dim(data)[2] - 1))
        {
            cutt <- cutree(sampleTree, k = i)
            p[i - 1] <- .cluster.stat(d = d, clustering = cutt)$pearsongamma
        }
        p <- unlist(na.omit(p))
        cutt <- cutree(sampleTree, k = which(p == max(p)) + 1)
    }
    else {
        cutt <- cutree(sampleTree, k = k)
        p <- .cluster.stat(d = d, clustering = cutt)$pearsongamma
        if (verbose) {
            message(.timestamp(), "-- NOTE: Number of 'k' subclasses defined
              by user have been settle down.")
        }
    }
    return(list(
        dend = sampleTree, coph = hmet,
        cluster = cutt, huber = max(p)
    ))
}


##############################
##############################
## probability subsampling

subsamplingProb <- function(x, iter, n1, n2 = 0) {
    if (n2 > 0) {
        p1o <- x / n1
        p2o <- x / n2

        m <- p1o * p2o

        r1 <- factorial(n1) / (factorial(x) * factorial(n1 - x))
        r2 <- factorial(n2) / (factorial(x) * factorial(n2 - x))

        m1 <- m * p2o * iter
        m2 <- m * p1o * iter

        m12 <- m * (p2o * p1o)^(x - 1)

        all <- r1 * r2
    } else {
        p1o <- x / n1
        r1 <- factorial(n1) / (factorial(x) * factorial(n1 - x))
        m <- 0
        all <- r1
    }
    return(list(prob.s = m, prob.comb = iter / all, n.all = all))
}


##############################
##############################
## Overlapping omic signal

overlapFeature <- function(id, data, classes, control, analysis,
                           infoS = NA, plot = FALSE) {
    if (analysis == "Binary") {
        da <- density.lf(as.matrix(data[id, names(classes[classes == control])]),
            from = range(data)[1],
            to = range(data)[2], n = 1000, width = 1
        )
        db <- density.lf(as.matrix(data[id, names(classes[classes != control])]),
            from = range(data)[1],
            to = range(data)[2], n = 1000, width = 1
        )
        d <- data.frame(x = da$x, a = da$y, b = db$y)

        # calculate intersection densities
        d$w <- pmin(d$a, d$b)

        # integrate areas under curves
        total <- integrate.xy(d$x, d$a) + integrate.xy(d$x, d$b)
        intersection <- integrate.xy(d$x, d$w)
        col <- adjustcolor(c("navyblue", "darkred"), alpha.f = 0.3)
        pos <- 2:dim(d)[2]
        txt <- c("Control", "Case")
    }
    if (analysis == "Multiclass") {
        d <- vapply(levels(classes), function(x) {
            error <- try(expr = {
                d <- density.lf(as.matrix(data[id, names(classes[classes == x])]),
                    from = range(data)[1], to = range(data)[2],
                    n = 1000, width = 1
                )$y
            }, silent = TRUE)
            if (inherits(error, "try-error")) {
                d <- density(as.numeric(data[id, names(classes[classes == x])]), n = 1000)$y
            }
            return(d)
        }, FUN.VALUE = numeric(1000))
        error <- try(expr = {
            x <- density.lf(as.matrix(data[id, names(classes[classes == levels(classes)[1]])]),
                from = range(data)[1], to = range(data)[2],
                n = 1000, width = 1
            )$x
        }, silent = TRUE)
        if (inherits(error, "try-error")) {
            x <- density(as.numeric(data[id, names(classes[classes == levels(classes)[1]])]), n = 1000)$x
        }
        w <- apply(apply(
            expand.grid(levels(classes), levels(classes)), 1,
            function(y) pmin(d[, y[1]], d[, y[2]])
        ), 1, max)
        d <- data.frame(d, x = x, w = w)
        total <- sum(vapply(seq_len(length(levels(classes))),
            function(y) integrate.xy(d$x, d[, y]),
            FUN.VALUE = numeric(1)
        ))^2
        intersection <- integrate.xy(d$x, d$w)

        pos <- seq_len(length(levels(classes)))
        col <- adjustcolor(colorRampPalette(c(
            "navyblue",
            "darkorange", "darkred", "darkgreen"
        ))(
            length(levels(classes))), alpha.f = 0.3)
        txt <- levels(classes)
    }
    # compute overlap coefficient
    overlap <- length(levels(classes)) * intersection / total

    if (plot) {
        par(mar = c(10, 6, 10, 8))
        stackpoly.2(
            x = d[, pos], col = col, border = "grey", axis1 = TRUE,
            stack = FALSE, lwd = 1, lty = 3,
            xaxlab = round(d[seq(1, dim(d)[1], dim(d)[1] * 0.05), "x"], 1),
            xat = seq(1, dim(d)[1], dim(d)[1] * 0.05),
            xlab = "Raw Data Signal", ylab = "Density",
            main = paste(
                "Overlap:", round(overlap, 3),
                "\nLikelihood density estimation, bandwidth = 1"
            )
        )
        legend("topright",
            legend = txt, pch = 15, cex = 1.5, pt.cex = 1.5,
            col = col, bty = "n"
        )
    }

    return(overlap)
}

##############################
##############################
## Colors assignment

jColor <- function(info) {
    info <- as.matrix(info)
    nam <- rownames(info)
    double <- c("ivory2", "gold", "gray30", "cornflowerblue", "chocolate", "honeydew1")

    colnames(info) <- paste(rev(LETTERS[seq_len(dim(info)[2])]), ": ", colnames(info), sep = "")
    for (j in seq_len(dim(info)[2]))
        info[, j] <- paste(rev(LETTERS[seq_len(dim(info)[2])])[j], ": ", as.character(info[, j]), sep = "")
    rownames(info) <- nam
    info.sample.color <- info
    info.sample.color <- apply(info.sample.color, 2, as.character)
    if (dim(info)[2] > 1) {
        pos <- apply(info, 2, unique)
        if (is.matrix(pos)) {
            pos <- rep(dim(pos)[1] == 2, dim(pos)[2])
        } else {
            pos <- lapply(pos, length) == 2
        }
    }
    else {
        pos <- dim(apply(info, 2, unique))[1] <= 2
    }
    ty <- sort(unlist(apply(info[, !pos, drop = FALSE], 2, unique)))
    tyy <- sort(unlist(apply(info[, pos, drop = FALSE], 2, unique)))
    myPalette <- c(brewer.pal(name = "Set1", n = 9)[seq_len(6)], "white", brewer.pal(name = "Set1", n = 9)[8:9], "black")
    if (length(ty) > 10) {
        myPalette <- colorRampPalette(myPalette)(length(unlist(apply(info, 2, unique))))
    }
    if (length(which(pos)) < dim(info)[2]) {
        for (z in seq_len(length(ty)))
        {
            info.sample.color[info %in% ty[z]] <-
                as.character(rep(myPalette[z], length(info.sample.color[info == unlist(apply(info, 2, unique))[z]])))
            names(ty)[z] <- unique(as.character(rep(myPalette[z], length(info.sample.color[info == unlist(apply(info, 2, unique))[z]]))))
        }
    }
    rownames(info.sample.color) <- nam
    ty <- c(tyy, ty)
    if (any(pos)) {
        for (j in seq_len((length(which(pos)) * 2)))
            info.sample.color[info == ty[j]] <- rep(double, 10)[j]
        names(ty)[seq_len(j)] <- rep(double, 10)[seq_len(j)]
    }

    return(list(orig = info, col = info.sample.color, ty = ty))
}


#################################
#################################
### Combinatorial calculation considering
### multiclass experimental design

.combCalcM <- function(data, classes, control,
                       iterations, multi, r, results) {
    if (!(all(names(classes) %in% colnames(data))) | is.null(names(classes))) {
        stop("ERROR: Names of samples in classes
         vector not match with data.")
    }

    classes <- sort(as.factor(classes))

    # If multiclass vector has been proposed:
    if (length(levels(classes)) > 2) {
        control <- NA
        combMulti <- t(combn(levels(classes), 2))

        # Removing auto-combination
        combMulti <- combMulti[combMulti[, 1] != combMulti[, 2], ]

        # Calculating maximal resampling size per class
        if (is.null(r) || r <= 0) {
            r <- round(sqrt(min(table(classes))), digits = 0)
        }
        maxSub <- subsamplingProb(
            x = r, n1 = min(table(classes)),
            iter = 0
        )$n.all
        if (iterations == 0) {
            iterations <- maxSub
        }
        if (is.na(maxSub)) {
            stop("ERROR: Resampling size is higher than minimum class size.")
        }
        multIter <- round(iterations / dim(combMulti)[1], 0)
        if (multIter > maxSub) {
            multIter <- maxSub
        }
        message(
            " NOTE: More than two classes of samples.\n
            Multiclass analysis must run ",
            dim(combMulti)[1],
            " (or multiple) iterations at least:\n ",
            multIter * dim(combMulti)[1], " random iterations (",
            multIter, " rounds) will be calculated.\n"
        )

        c1 <- apply(combMulti, 1, function(y) names(classes[classes == y[1]]))
        c2 <- apply(combMulti, 1, function(y) names(classes[classes == y[2]]))
        combi <- matrix(data = 0, ncol = 2 * r, nrow = 1)

        # Calculating combinations of samples without mixing classes.
        for (i in seq_len(multIter)) {
            if (is.list(c1)) {
                combi1 <- matrix(unlist(lapply(c1, function(x)
                    which(colnames(data) %in% sample(x, size = r, replace = FALSE)))),
                ncol = r, byrow = TRUE
                )
                combi2 <- matrix(unlist(lapply(c2, function(x)
                    which(colnames(data) %in% sample(x, size = r, replace = FALSE)))),
                ncol = r, byrow = TRUE
                )
            } else {
                combi1 <- t(apply(c1, 2, function(x)
                    which(colnames(data) %in% sample(x, size = r, replace = FALSE))))
                combi2 <- t(apply(c2, 2, function(x)
                    which(colnames(data) %in% sample(x, size = r, replace = FALSE))))
            }
            combi <- rbind(combi, cbind(combi1, combi2))
        }
        results <- cbind(combi[2:dim(combi)[1], ],
            nFeatures = c(rep(0, multIter * dim(combMulti)[1]))
        )
        colnames(results) <- c(seq_len(2 * r), "nFeatures")
        n1 <- dim(data)[2]
        n2 <- min(table(classes))
        cl1 <- NA
        cl2 <- NA
        categories.control <- seq_len(dim(data)[2])
        categories.case <- seq_len(dim(data)[2])
        multi <- TRUE
    } else {
        # Binary contrast. 'Control' label defines
        # '0' class for eBayes algorithm.
        if (is.na(control)) {
            cl1 <- levels(classes)[1]
            cl2 <- levels(classes)[2]
            control <- cl1
        } else {
            if (!(control %in% levels(classes))) {
                stop("Control label not found in vector classes provided.")
            }
            cl1 <- control
            cl2 <- levels(classes)[which(levels(classes) != control)]
        }
        n1 <- length(which(classes == cl1))
        n2 <- length(which(classes == cl2))
        classes <- classes[c(
            which(classes == cl1),
            which(classes == cl2)
        )]

        # Maximum number of iterations
        maxSub <- prod(
            factorial(c(n1, n2)) /
                (factorial(r) * factorial(c(n1, n2) - r))
        )

        if (iterations > maxSub) {
            iterations <- maxSub
            message(
                .timestamp(), " -- 'iterations' is higher than
              maximum number of combinations of samples.
              Number of 'iterations' will be changed to ", maxSub,
                "."
            )
        }
        data <- data[, names(classes)]
        categories.control <- seq_len(n1)
        categories.case <- (n1 + 1):(n1 + n2)
        message("\n Classes vector defined as input.
            SUPERVISED analysis will be carry out.\n")
    }

    return(list(
        categories.control = categories.control, data = data,
        categories.case = categories.case, results = results,
        multi = multi, classes = classes, cl1 = cl1, cl2 = cl2,
        n1 = n1, n2 = n2
    ))
}


#################################
#################################
### Combinatorial calculation

.combCalc <- function(iterations, multi, classes, r, results,
                      categories.control, categories.case,
                      n1, n2) {
    cont <- 0
    user <- "n"
    # Asking for number of iterations if it was not defined.
    if (iterations == 0 || is.null(iterations) & !multi) {
        n.all <- subsamplingProb(x = r, iter = 1000, n1 = n1, n2 = n2)$n.all
        user <- readline(prompt = message(
            "\n All possible combinations:", n.all,
            "\n Run ALL combinations? (y/n):"
        ))
        if (user == "n") {
            iterations <- as.integer(readline(
                prompt = message("\n Define number of random combinations to calculate:")
            ))
            if (!(is.integer(iterations))) {
                stop("Non-integer argument passed as number of combinations.")
            }
        }
    }
    # Calculate random or all possible combinations between samples if 'binary',
    # 'multiclass' or 'unsupervised' has been proposed.
    if (user == "y" & !multi) {
        allCombControl <- combn(n1, r)
        allCombCase <- combn(n2, r)
        nSamples <- nrow(allCombControl) + nrow(allCombCase)
        ncomb <- ncol(allCombControl) * ncol(allCombCase)
        results <- matrix(nrow = ncomb, ncol = nSamples)
        colnames(results) <- seq_len(nSamples)
        cont <- 1

        # Building combination matrix.
        for (iCombCONTROL in seq_len(ncol(allCombControl)))
        {
            for (iCombCASE in seq_len(ncol(allCombCase)))
            {
                samples <- c(
                    categories.control[allCombControl[, iCombCONTROL]],
                    categories.case[allCombCase[, iCombCASE]]
                )
                if (length(which(duplicated(samples))) > 0) {
                    samples[duplicated(samples)] <-
                        sample(which(!(seq_len(n1 + n2) %in% samples)), 1)
                }
                results[cont, seq_len(nSamples)] <- samples
                cont <- cont + 1
            }
        }
        results <- cbind(results, nFeatures = c(rep(0, ncomb)))
        msg <- paste(.timestamp(), "-- Number of total iterations:", ncomb)
    } else if (!multi) {
        results <- matrix(nrow = iterations, ncol = 2 * r + 1)
        colnames(results) <- c(seq(from = 1, to = 2 * r), "nFeatures")
        # Building combination matrix.
        for (i in seq_len(iterations))
        {
            if (all(!(is.na(classes)))) {
                results[i, seq_len(r)] <- sample(categories.control, r, replace = FALSE)
                results[i, (r + 1):(2 * r)] <- sample(categories.case, r, replace = FALSE)
            }
            else {
                results[i, seq_len(r * 2)] <- sample(categories.control, 2 * r, replace = FALSE)
            }
        }
        results[, c("nFeatures")] <- 0
        msg <- paste(.timestamp(), "-- Randomly selected", iterations, "iterations.")
    } else {
        msg <- paste(.timestamp(), "-- Number of total iterations:", dim(results)[1])
    }

    message(msg)

    return(results)
}

#################################
#################################
### Internal LIMMA calculations

.LIMMAcalc <- function(data, results, classes, control, q.val,
                       call, r, cpus, parallel, temp.path,
                       multi, n1, n2) {
    # Creating final variable containing all subsampling results.
    res <- list(
        data = as.matrix(data), results = results,
        subStatFeature = NULL, incidenceMatrix = NULL,
        classes = classes, resampleSize = r, control = control,
        pos.iter = 0, q.val = q.val, call = call
    )
    # Some variables definition to parallel subsampling.
    top_eje <- matrix(ncol = 7)
    colnames(top_eje) <- c(
        "logFC", "AveExpr", "t", "P.Value",
        "adj.P.Val", "B", "ID"
    )
    top_eje[, 2:7] <- 0
    UP <- as.data.frame(data)
    if (all(is.na(classes)) | multi) {
        UP[, seq_len(n1)] <- 0
    } else {
        UP[, seq_len(n1 + n2)] <- 0
    }
    DOWN <- UP
    mx <- UP
    results <- cbind(results, counter = seq(seq_len(dim(results)[1])))

    suppressWarnings(limma1 <- bplapply(seq_len(dim(results)[1]),
        FUN = .limmaSubsamp,
        results = results[, seq_len(2 * r)], r = r,
        data = data, q.val = q.val,
        temp.path = temp.path
    ))

    pos <- unlist(lapply(limma1, function(x) x$pos))
    limma1 <- unlist(lapply(limma1, function(x) x$p))
    limma1 <- limma1[!(is.na(limma1))]

    results[, "nFeatures"] <- pos
    results <- results[, seq_len(dim(results)[2] - 1)]

    # Stop subsampling for non-differential signal.
    if (length(limma1) <= 1) {
        stop("NO-RESULT: Any combination shows DE features with these
          resampling size and p-value threshold.
          Try again with another parameters.")
    }

    return(list(
        res = res, limma1 = limma1,
        top_eje = top_eje, mx = mx,
        UP = UP, DOWN = DOWN,
        results = results
    ))
}


#################################
#################################

.statCalc <- function(counter, tab, top_eje, ncomb, j) {
    res_ <- vector(length = 12)
    names(res_) <- colnames(tab)

    lines <- top_eje[as.character(top_eje[, c("ID")]) %in%
        as.character(tab[counter, c("ID")]), ]

    # Statistics from this RDA step.
    if (length(lines[, c("ID")]) > 0) {
        res_[c("Avrg.logFC")] <- .colMeans(lines[, c("logFC")],
            m = length(lines[, c("ID")]), n = 1
        )
        res_[c("Best.adj.P.Val")] <- lines[1, c("P.Value")]
        res_[c("Repeats")] <- length(lines[, 1])
        res_[c("FR.Repeats")] <- res_[c("Repeats")] / ncomb
        res_[c("RF.Positive.Repeats")] <- res_[c("Repeats")] / j
        res_[c("Chi.Square")] <- (-2) * sum(log(lines[, c("P.Value")]))
        res_[c("P.Val.ChiSq")] <- 1 - pchisq(as.numeric(res_[c("Chi.Square")]),
            df = 2 * length(lines[, 1])
        )
        res_[c("ID")] <- as.character(tab[counter, c("ID")])
    }

    return(res_)
}


#################################
#################################

.limmaSubsamp <- function(counter, results, r, data, q.val,
                          temp.path) {
    # Taking different combinations from previous matrix and running LIMMA.
    CASE <- CONTROL <- NULL
    samples <- results[counter, seq_len(r * 2)]
    design <- cbind(CONTROL = rep(c(1, 0), each = r), CASE = rep(c(0, 1), each = r))
    fit <- limma::lmFit(data[, samples], design)
    cont.mat <- limma::makeContrasts(CASEvsCONTROL = CASE - CONTROL,
                                     levels = design)
    fit2 <- limma::contrasts.fit(fit, cont.mat)
    fit3 <- limma::eBayes(fit2)
    top <- limma::topTable(fit3,
        adjust.method = "none",
        number = nrow(data),
        p.value = q.val
    )
    pos <- dim(top)[1]

    # Generating intermediate files in temporary dir. If any DE feature is found for
    # any combination, no file would be created.
    resultfile <- paste(temp.path, "/", c("diff", "incid"),
        counter, ".dta",
        sep = ""
    )
    if (pos > 0) {
        top <- data.frame(top, ID = rownames(top))
        top <- top[seq_len(pos), ]
        # Differential features and incidence matrix files are separated.
        foreign::write.dta(top, file = resultfile[1])
        foreign::write.dta(as.data.frame(samples), file = resultfile[2])
    }
    else {
        resultfile <- c(NA, NA)
    }

    return(list(p = resultfile, pos = pos))
}


#################################
#################################
### Frequency matrix and intermediate results

.freqMatrix <- function(limmaRes, data, n1, n2, r,
                        multi, unsup) {
    mx0 <- limmaRes$mx
    mx20 <- limmaRes$mx
    pb <- txtProgressBar(style = 3, min = 1, max = length(limmaRes$limma1) - 1)
    UP <- limmaRes$UP
    DOWN <- limmaRes$DOWN
    top_eje <- limmaRes$top_eje

    # Reading and summarizing intermediate results.
    for (i in seq(from = 1, to = length(limmaRes$limma1) - 1, by = 2))
    {
        mx <- mx0
        mx2 <- mx20
        resultfile1 <- limmaRes$limma1[i]
        resultfile2 <- limmaRes$limma1[i + 1]
        if (is.matrix(limmaRes$limma1)) {
            resultfile1 <- limmaRes$limma1[1, i]
            resultfile2 <- limmaRes$limma1[2, i]
        }

        # Joining all differential event statistics.
        top_eje1 <- read.dta(file = resultfile1)
        colnames(top_eje1) <- colnames(top_eje)
        rownames(top_eje1) <- make.names(rep(LETTERS, dim(data)[1])[
            seq_len(dim(top_eje1)[1])
        ], unique = TRUE)
        samples <- as.vector(t(read.dta(file = resultfile2)))
        counter <- as.numeric(unlist(strsplit(unlist(
            strsplit(resultfile1, split = "diff")
        )[2], split = ".", fixed = TRUE))[1])

        # Generating global incidence matrix with UP and DOWN events.
        mx[rownames(mx) %in% top_eje1[top_eje1[, c("logFC")] < 0, c("ID")], samples[seq_len(r)]] <- mx[
            rownames(mx) %in% top_eje1[top_eje1[, c("logFC")] < 0, c("ID")], samples[seq_len(r)]
        ] + 1
        mx[rownames(mx) %in% top_eje1[top_eje1[, c("logFC")] > 0, c("ID")], samples[(r + 1):(r * 2)]] <- mx[
            rownames(mx) %in% top_eje1[top_eje1[, c("logFC")] > 0, c("ID")], samples[(r + 1):(r * 2)]
        ] + 1
        mx2[rownames(mx2) %in% top_eje1[top_eje1[, c("logFC")] < 0, c("ID")], samples[(r + 1):(r * 2)]] <- mx2[
            rownames(mx2) %in% top_eje1[top_eje1[, c("logFC")] < 0, c("ID")], samples[(r + 1):(r * 2)]
        ] + 1
        mx2[rownames(mx2) %in% top_eje1[top_eje1[, c("logFC")] > 0, c("ID")], samples[seq_len(r)]] <- mx2[
            rownames(mx2) %in% top_eje1[top_eje1[, c("logFC")] > 0, c("ID")], samples[seq_len(r)]
        ] + 1

        UP <- UP + mx
        DOWN <- DOWN + mx2
        top_eje <- rbind(top_eje, top_eje1)
        rownames(top_eje1) <- NULL
        setTxtProgressBar(pb, i)
    }

    close(pb)

    UP <- UP[rowSums(UP) > 0, ]
    DOWN <- DOWN[rowSums(DOWN) > 0, ]
    MULTI <- UP + DOWN

    # Making up incidence matrix with separated UP, DOWN and MIXED events.
    if (!unsup & !multi) {
        both.feat <- names(which(rowSums(UP[, seq_len(n1)]) > 0 &&
            rowSums(UP[, (n1 + 1):(n1 + n2)]) > 0))
        if (length(both.feat) > 0) {
            both.feat.mx <- UP[rowSums(UP) > 0, ][rownames(UP) %in% both.feat, ]
            rownames(both.feat.mx) <- paste(rownames(both.feat.mx), sep = "deco", "UP")
            UP <- UP[!(rownames(UP) %in% both.feat), ]
            DOWN <- DOWN[!(rownames(DOWN) %in% both.feat), ]
        } else {
            both.feat.mx <- matrix(0, nrow = 1, ncol = ncol(UP))
        }
    } else {
        both.feat <- NA
        both.feat.mx <- NA
    }

    # Vector containing direction of differential event if 'binary'
    UpDw <- c(
        rep("UP", length(rownames(UP)[rowSums(UP[, seq_len(n1)]) == 0])),
        rep("DOWN", length(rownames(UP)[rowSums(UP[, seq_len(n1)]) > 0])),
        rep("MIXED", length(both.feat))
    )
    names(UpDw) <- c(
        rownames(UP)[rowSums(UP[, seq_len(n1)]) == 0],
        rownames(UP)[rowSums(UP[, seq_len(n1)]) > 0], both.feat
    )
    rownames(UP) <- paste(rownames(UP), sep = "deco", "UP")
    rownames(MULTI) <- paste(rownames(MULTI), sep = "deco", "UP")
    rownames(DOWN) <- paste(rownames(DOWN), sep = "deco", "DOWN")

    return(list(
        top_eje = top_eje, UP = UP, DOWN = DOWN,
        MULTI = MULTI, UpDw = UpDw,
        both.feat = both.feat,
        both.feat.mx = both.feat.mx
    ))
}


#################################
#################################
### Function for filtering features
### by Repeats

.repThr <- function(sub, rep.thr, samp.perc, cpus) {
    suppressWarnings(g.names <- unlist(bplapply(rownames(sub$incidenceMatrix),
        FUN = function(x) unlist(strsplit(x, split = "deco", fixed = TRUE))[1]
    )))
    names(g.names) <- rownames(sub$incidenceMatrix)

    if (all(is.na(sub$classes)) || length(levels(sub$classes)) > 2) {
        f <- apply(sub$incidenceMatrix, 1, function(x) length(which(x >= rep.thr)))
        names(f) <- g.names[names(g.names) %in% names(f)]
    }
    else {
        f <- vapply(unique(g.names), function(x)
            sum(apply(sub$incidenceMatrix[grepl(
                rownames(sub$incidenceMatrix),
                pattern = x, fixed = TRUE
            ), ], 1, function(x)
                length(which(x >= rep.thr)))), numeric(1))
    }

    # 'Repeats' threshold.
    # All features under those thresholds will be removed.
    if (length(which(f > samp.perc * dim(sub$incidenceMatrix)[2])) > 0) {
        message(paste(
            " NOTE: Repeats index threshold have been applied:\n",
            length(which(f <= samp.perc * dim(sub$incidenceMatrix)[2])),
            "features removed from", dim(sub$subStatFeature)[1], "total."
        ))
    }
    sub$subStatFeature <- sub$subStatFeature[names(
        f[which(f > samp.perc * dim(sub$incidenceMatrix)[2])]
    ), ]
    sub$subStatFeature <- cbind(
        sub$subStatFeature,
        Repeats.index = f[names(
            f
        ) %in% rownames(sub$subStatFeature)] / dim(sub$incidenceMatrix)[2] * 100
    )
    sub$subStatFeature <- sub$subStatFeature[order(
        rownames(sub$subStatFeature),
        decreasing = TRUE
    ), ]

    # Removing features from incidenceMatrix with lower 'Repeats' than threshold.
    sub$incidenceMatrix <- sub$incidenceMatrix[rownames(
        sub$incidenceMatrix
    ) %in% names(
        g.names[g.names %in% as.character(sub$subStatFeature$ID)]
    ), ]
    if (length(which(f > samp.perc * dim(sub$incidenceMatrix)[2])) < 10) {
        stop("After applying 'Repeats' filter, there are not enough features
         (10 features minimum) to input NSCA.")
    }

    return(list(sub = sub, g.names = g.names))
}


#################################
#################################
### Function to call overlapFeature

overlapC <- function(sub, cl1, cpus) {
    message(.timestamp(), " -- Calculating overlapping signal per feature...")

    ## Calculating overlap...
    suppressWarnings(overlap <- bplapply(
        X = rownames(sub$subStatFeature),
        FUN = overlapFeature,
        data = sub$data, classes = sub$classes,
        control = cl1, analysis = "Binary"
    ))
    overlap <- unlist(overlap)
    names(overlap) <- rownames(sub$subStatFeature)

    ## Type of profile...
    prof <- vapply(overlap, function(x)
        max(which(c(0, 0.2, 0.4, 0.75) <= x)),
    FUN.VALUE = numeric(1)
    )

    profile <- prof
    profile[prof %in% c(3, 4)] <- "Minority"
    profile[prof == 2] <- "Majority"
    profile[prof == 1] <- "Complete"
    profile[sub$subStatFeature[names(profile), c("UpDw")] == "MIXED"] <- "Mixed"
    names(profile) <- names(overlap)
    profile <- profile[order(names(profile))]

    return(list(overlap = overlap, profile = profile))
}


#################################
#################################
### Intern R function of decoReport

.formattingObj <- function(deco, deco0,
                           analysis, info.size) {
    ## Formatting R objects to print.
    if (analysis == "Binary") {
        # p.value from NSCA analysis.
        p.val <- c(
            deco@NSCAcluster$Control$NSCA$P.Value,
            deco@NSCAcluster$Case$NSCA$P.Value
        )

        # General information about subclasses found.
        infoSubclass <- rbind(
            deco@NSCAcluster$Control$infoSubclass,
            deco@NSCAcluster$Case$infoSubclass
        )
        rel <- c(
            rep(
                p.val[1] <= 0.05,
                dim(deco@NSCAcluster$Control$infoSubclass)[1]
            ),
            rep(p.val[2] <= 0.05, dim(deco@NSCAcluster$Case$infoSubclass)[1])
        )
        infoSubclass <- data.frame(infoSubclass,
            Binary = rep(0, dim(infoSubclass)[1]),
            isRelevant = as.character(rel)
        )
        infoSubclass[vapply(
            rownames(infoSubclass),
            function(x) unlist(strsplit(x, split = " Subclass"))[1],
            character(1)
        ) != deco@control, "Binary"] <- 1

        # Sample membership to subclasses.
        samplesSubclass <- rbind(
            cbind(deco@NSCAcluster$Control$samplesSubclass[order(
                deco@NSCAcluster$Control$samplesSubclass[, 1]
            ), ]),
            cbind(deco@NSCAcluster$Case$samplesSubclass[order(
                deco@NSCAcluster$Case$samplesSubclass[, 1]
            ), ])
        )

        samplesSubclass <- cbind(
            Samples = rownames(samplesSubclass),
            Subclass = samplesSubclass[, 1]
        )
        rownames(samplesSubclass) <- seq_len(length(rownames(samplesSubclass)))

        # NSCA and hierarchical clustering information.
        nsc <- matrix(round(c(
            deco@NSCAcluster$Control$var, deco@NSCAcluster$Case$var,
            p.val, deco@NSCAcluster$Control$hclustSamp$huber,
            deco@NSCAcluster$Case$hclustSamp$huber
        ), digits = 3),
        nrow = 3, byrow = TRUE
        )
        colnames(nsc) <- c("Control samples", "Case samples")
        rownames(nsc) <- c(
            "Variability explained by NSCA",
            "NSCA C-statistic p.value",
            "Huber's gamma"
        )

        # Assignment of colors to all subclasses. It will be conserved along all the report.
        count <- table(vapply(
            rownames(infoSubclass),
            function(x) unlist(strsplit(x, split = " Subclass"))[1],
            character(1)
        ) != deco@control)
        color.cluster <- c(
            colorRampPalette(c("navyblue", "lightblue"))(count[1]),
            colorRampPalette(c("orange", "darkred"))(count[2])
        )
        names(color.cluster) <- rownames(infoSubclass)

        # Columns of 'featureTable' to be printed.
        names.col <- c(
            "ID", "SYMBOL", "UpDw", "Profile", "overlap.Ctrl.Case",
            "Standard.Chi.Square", "Repeats",
            "Repeats.index", "delta.signal",
            "sd.Ctrl", "Dendrogram.group.Ctrl",
            "h.Range.Ctrl", "sd.Case",
            "Dendrogram.group.Case", "h.Range.Case"
        )

        deco@featureTable[, names.col[!names.col %in% c("ID", "SYMBOL", "UpDw", "Profile")]] <- apply(
            deco@featureTable[, names.col[!names.col %in% c("ID", "SYMBOL", "UpDw", "Profile")]], 2, as.numeric
        )

        textF <- deco@featureTable[, colnames(deco@featureTable) %in% names.col]
        deco0@featureTable <- deco@featureTable
        textF <- textF[, na.omit(match(names.col, colnames(textF)))]

        # Chi, Delta , Repeats & SD(inverse)
        ordG <- order(rowMedians(apply(cbind(
            textF$Standard.Chi.Square, -textF$overlap.Ctrl.Case,
            -textF$overlap.Ctrl.Case, textF$Repeats, abs(textF$delta.signal),
            abs(textF$delta.signal), -textF$sd.Ctrl, -textF$sd.Case
        ), 2, function(x)
            rank(-x, ties.method = "max"))))
        # Chi, Delta & SD
        ordS <- rank(-rowMedians(apply(cbind(
            textF$Standard.Chi.Square, textF$Standard.Chi.Square,
            abs(textF$delta.signal), textF$Repeats,
            abs(textF$delta.signal), textF$h.Range.Ctrl, textF$h.Range.Case
        ), 2, function(x)
            rank(-x, ties.method = "max"))))

        # Saving initial modified 'featureTable' as R data.frame called 'textF0'.
        deco@featureTable <- deco@featureTable[ordG, ]
        textF0 <- textF
        textF <- textF[ordG, ]

        # Object defining samples membership to classes previously compared.
        classes <- unlist(lapply(deco@NSCAcluster, function(x)
            unlist(strsplit(x$samplesSubclass[1], split = " Subclass"))[1]))
    } else {
        # p.value from NSCA analysis.
        p.val <- deco@NSCAcluster$All$NSCA$P.Value

        # General information about subclasses found.
        infoSubclass <- deco@NSCAcluster$All$infoSubclass
        rel <- rep(p.val[1] <= 0.05, dim(
            deco@NSCAcluster$All$infoSubclass
        )[1])
        infoSubclass <- data.frame(infoSubclass,
            isRelevant = as.character(rel)
        )

        # Sample membership to subclasses.
        samplesSubclass <- cbind(deco@NSCAcluster$All$samplesSubclass[order(
            deco@NSCAcluster$All$samplesSubclass[, 1]
        ), ])

        samplesSubclass <- cbind(
            Samples = rownames(samplesSubclass),
            Subclass = samplesSubclass[, 1]
        )
        rownames(samplesSubclass) <- seq_len(length(rownames(samplesSubclass)))

        # NSCA and hierarchical clustering information.
        nsc <- round(rbind(
            deco@NSCAcluster$All$var,
            p.val, deco@NSCAcluster$All$hclustSamp$huber
        ), 3)
        colnames(nsc) <- c("All samples")
        rownames(nsc) <- c(
            "Variability explained by NSCA",
            "p.value", "Huber's gamma"
        )

        # Assignment of colors to all subclasses. It will be conserved along all the report.
        color.cluster <- colorRampPalette(c("greenyellow", "lightblue", "darkgreen", "darkred"))(
            dim(infoSubclass)[1])
        names(color.cluster) <- rownames(infoSubclass)

        # Columns of 'featureTable' to be printed.
        names.col <- c(
            "ID", "SYMBOL", "Standard.Chi.Square", "Repeats",
            "Repeats.index", "Avrg.logFC", "Dendrogram.group",
            "h", "sd", "h.Range"
        )


        ind <- colnames(deco@featureTable) %in%
            names.col[!names.col %in% c("ID", "SYMBOL")]
        deco@featureTable[, ind] <- apply(deco@featureTable[, ind], 2, as.numeric)

        textF <- deco@featureTable[, colnames(deco@featureTable) %in% names.col]
        textF <- textF[, na.omit(match(names.col, colnames(textF)))]

        ordG <- order(rowMedians(apply(
            cbind(
                textF$Standard.Chi.Square, textF$Standard.Chi.Square,
                textF$Repeats, textF$Repeats, textF$h.Range
            ), 2,
            function(x) rank(-x, ties.method = "max")
        )))
        ord <- ordG

        # No sense to make two different ranking for "Unsupervised" analysis.
        deco@featureTable <- deco@featureTable[ord, ]
        textF0 <- textF
        textF <- textF[ord, ][seq_len(info.size), ]
        classes <- sub$classes
    }
    return(list(
        deco = deco, p.val = p.val, color.cluster = color.cluster,
        infoSubclass = infoSubclass, ordG = ordG, ordS = ordS,
        rel = rel, classes = classes, nsc = nsc,
        textF = textF, textF0 = textF0, names.col = names.col,
        samplesSubclass = samplesSubclass
    ))
}


#################################
#################################
###

NSCACluster <- function(mx, data = NULL, id.names = NULL, k = NULL,
                        label = NA, v = 80, method.dend = NULL,
                        raw = NULL, UpDw = NULL, dir = NA) {
    if (is.null(id.names))
        id.names <- rownames(mx)

    ca_res <- NSCA(as.matrix(mx[rowSums(mx) > 0, colSums(mx) > 0]),
        v = v
    )
    nd <- min(which(ca_res$Inertia[, 3] >= v))
    if (nd == 1) {
        nd <- 2
        message(.timestamp(), "-- 1D is enough to reach explained variability from input.
            2D will be considered")
    }
    ca_ev <- rbind(ca_res$fbip[, seq_len(nd)], ca_res$g[, seq_len(nd)])

    # Renaming feature IDs
    if (all(!is.null(id.names))) {
        rownames(ca_ev)[rownames(ca_ev) %in% names(id.names)] <- as.character(id.names[
            rownames(ca_ev)[rownames(ca_ev) %in% names(id.names)]
        ])
    }

    # Calculating h statistic
    if (all(!is.null(raw)) & all(!is.null(data))) {
        dispersion <- data[id.names[rownames(ca_res$inner.prod)], colnames(ca_res$inner.prod)] -
            raw[id.names[rownames(ca_res$inner.prod)]]

        # Dispersion inversion for DOWN regulated within a category of samples
        if (all(!is.null(UpDw))) {
            dispersion[rownames(dispersion) %in% id.names[
              id.names %in% names(UpDw[UpDw == dir])], ] <- -dispersion[
                rownames(dispersion) %in% id.names[id.names %in% names(UpDw[UpDw == dir])],
            ]
        }

        inner <- as.matrix(ca_res$inner.prod)
        inner[] <- scale(as.vector(inner))
        dispersion[] <- scale(as.vector(dispersion))

        # Make positive ranges
        inner <- inner + abs(min(inner)) + 1
        rownames(dispersion) <- rownames(inner)
        ca_res$inner.prod <- as.matrix(inner * dispersion - mean(as.matrix(inner * dispersion)))

        # 'h' inversion for DOWN regulated within a category of samples
        if (all(!is.null(UpDw))) {
            ca_res$inner.prod[rownames(ca_res$inner.prod) %in% names(
                id.names
            )[id.names %in% names(UpDw[UpDw == dir])], ] <- -ca_res$inner.prod[
                rownames(ca_res$inner.prod) %in% names(id.names)[
                    id.names %in% names(UpDw[UpDw == dir])
                ],
            ]
        }
    } else if (all(!is.null(data))) {
        inner <- as.matrix(ca_res$inner.prod)
        inner[] <- scale(as.vector(inner))
        dispersion <- data[id.names[rownames(ca_res$inner.prod)], colnames(ca_res$inner.prod)] -
            rowMeans(data[id.names[rownames(ca_res$inner.prod)], colnames(ca_res$inner.prod)])
        dispersion[] <- scale(as.vector(dispersion))

        # Make positive ranges
        inner <- inner + abs(min(inner)) + 1

        # dispersion <- dispersion + abs(min(dispersion)) + 1
        ca_res$inner.prod <- as.matrix(inner * dispersion - mean(as.matrix(inner * dispersion)))
    } else {
        message("Data has not been provided, then 'h statistic' will be not computed.
        Inner product from NSCA\n
        would be used in biclustering of features and samples.")
    }
    # Make sure assignment
    rownames(ca_res$inner.prod) <- rownames(inner)

    # Biclustering
    info.dend.samp <- cophDECO(
        data = as.matrix(ca_res$inner.prod), method.heatmap = method.dend,
        k = k, scale = FALSE, coph = TRUE
    )
    info.dend.feat <- cophDECO(
        data = as.matrix(t(ca_res$inner.prod)), method.heatmap = "ward.D",
        k = 10, scale = FALSE, verbose = FALSE, coph = FALSE
    )
    k <- max(info.dend.samp$cluster)

    # Calculating some statistics...
    kk <- gg1 <- gg2 <- ss <- vector(length = k)
    message(.timestamp(), "-- Optimized for ", k, " subclasses.")

    info <- vapply(seq_len(k), function(y) rowMeans(cbind(
            as.matrix(ca_res$inner.prod)[, info.dend.samp$cluster == y]
        )),
    FUN.VALUE = numeric(dim(ca_res$inner.prod)[1])
    )
    info <- cbind(info, ca_res$di[rownames(info)])

    if (k > 1)
        info <- as.data.frame(cbind(info, matrix(nrow = dim(info)[1], ncol = 3)))
    else
        info <- as.data.frame(t(rbind(info, matrix(nrow = dim(info)[1], ncol = 3))))

    colnames(info) <- c(
        paste("Scl", seq_len(k), sep = ""),
        "Tau.feature", "Closer.subclass", "h.Best", "ID"
    )
    info <- info[order(rownames(info), decreasing = TRUE), ]

    if (k > 1) {
        info[, c("Closer.subclass")] <- apply(
            info[, seq_len(k)], 1,
            function(x) which(abs(x) == max(abs(x)))
        )
    } else {
        info[, c("Closer.subclass")] <- rep(1, dim(info)[1])
    }

    for (i in seq_len(k))
    {
        kk[i] <- length(which(colnames(mx) %in% names(
            info.dend.samp$cluster[info.dend.samp$cluster == i]
        )))
        gg1[i] <- length(which(info[which(vapply(rownames(info), function(x) unlist(
                strsplit(x, split = "deco", fixed = TRUE)
            )[2], character(1)) == "UP"), "Closer.subclass"] == i))
        gg2[i] <- length(which(info[which(vapply(rownames(info), function(x) unlist(
                strsplit(x, split = "deco", fixed = TRUE)
            )[2], character(1)) == "DOWN"), "Closer.subclass"] == i))
        p <- abs(info[info[, "Closer.subclass"] == i, i])
        ss[i] <- mean(p[order(p, decreasing = TRUE)][seq_len((length(p) * 0.05))])
    }
    if (all(!is.na(label))) {
        infoSubclass <- data.frame(
            Samples = kk, FeaturesUP = gg1, FeaturesDOWN = gg2, SpeValue = ss,
            row.names = paste(label, "Subclass", seq(from = 1, to = k), sep = " ")
        )
    } else {
        infoSubclass <- data.frame(
            Samples = kk, FeaturesUP = gg1, FeaturesDOWN = gg2, SpeValue = ss,
            row.names = paste("Subclass", seq(from = 1, to = k), sep = " ")
        )
    }
    if (k > 1) {
        info[, "h.Best"] <- apply(
            info[, c(paste("Scl", seq_len(k), sep = ""), "Closer.subclass")],
            1, function(x) x[x[k + 1]]
        )
    } else {
        info[, "h.Best"] <- info[, 1]
    }

    # Renaming some feature IDs and removing duplicated features.
    info[, c("ID")] <- as.character(id.names[rownames(info)])
    rownames(ca_ev)[rownames(ca_ev) %in% rownames(mx)][!(
        duplicated(id.names[rownames(ca_ev[rownames(ca_ev) %in% rownames(mx), ])]) | duplicated(id.names[
            rownames(ca_ev[rownames(ca_ev) %in% rownames(mx), ])
        ], fromLast = TRUE))] <- id.names[rownames(ca_ev[rownames(ca_ev) %in% rownames(mx), ])][!(
        duplicated(id.names[rownames(ca_ev[rownames(ca_ev) %in% rownames(mx), ])]) | duplicated(id.names[rownames(ca_ev[rownames(ca_ev) %in% rownames(mx), ])], fromLast = TRUE))]
    info <- info[!rownames(info) %in% rownames(info[info[, c("ID")] %in% info[duplicated(info[, c("ID")]), c("ID")], ])[
        which(unlist(strsplit(rownames(info[info[, c("ID")] %in% info[duplicated(info[, c("ID")]), c("ID")], ]),
            split = "deco", fixed = TRUE
        )) == "DOWN") / 2
    ], ]
    rownames(ca_res$inner.prod)[!(duplicated(
        id.names[rownames(ca_res$inner.prod)]
    ) | duplicated(id.names[rownames(ca_res$inner.prod)], fromLast = TRUE))] <- id.names[
        rownames(ca_res$inner.prod)
    ][!(duplicated(id.names[rownames(ca_res$inner.prod)]) | duplicated(id.names[
        rownames(ca_res$inner.prod)
    ], fromLast = TRUE))]

    # Calculating h statistic per feature per subclass of samples.
    ord <- seq_len(length(info.dend.feat$cluster))
    names(ord) <- rownames(ca_res$inner.prod)[rev(info.dend.feat$dend$order)]
    patt <- info.dend.feat$cluster
    for (i in seq_len(max(info.dend.feat$cluster)))
        patt[info.dend.feat$cluster == i] <- which(order(vapply(seq_len(max(info.dend.feat$cluster)), function(x)
            mean(which(info.dend.feat$cluster[rev(info.dend.feat$dend$order)] == x)), numeric(1))) == i)
    info.dend.feat$cluster <- patt
    info <- data.frame(info,
        h.Range = apply(info[, seq_len(k)], 1, function(x) diff(range(x))),
        Dendrogram.group = info.dend.feat$cluster[rownames(info)], Dendrogram.order = ord[rownames(info)]
    )
    rownames(info) <- as.character(info[, c("ID")])
    rankingH <- cbind(
        t(interleave(
            t(apply(info[, seq_len(k)], 2, function(x) rank(-abs(x), ties.method = "max"))),
            t(info[, seq_len(k)])
        )), info[, "h.Range"],
        info[, "Dendrogram.group"], ord[rownames(info)]
    )
    colnames(rankingH) <- c(
        paste(rep(c("Ranking", "h"), k), rep(paste("Scl", seq_len(k), sep = ""), each = 2), sep = "."),
        "h.Range", "Dendrogram.group", "Dendrogram.order"
    )

    # Find out sample membership to any subclass.
    samplesSubclass <- cbind(info.dend.samp$cluster)
    colnames(samplesSubclass) <- "Subclass"
    for (i in seq_len(dim(samplesSubclass)[1]))
        samplesSubclass[i, 1] <- rownames(infoSubclass)[as.numeric(samplesSubclass[i, 1])]

    # Returning results.
    if (all(!is.null(data))) {
        names(ca_res)[names(ca_res) == "inner.prod"] <- "h"
        results <- list(
            ca_res, info, rankingH, infoSubclass, mx, ca_ev, samplesSubclass, ca_res$Inertia[nd, 3],
            method.dend, info.dend.samp, info.dend.feat, inner, dispersion
        )
        names(results) <- c(
            "NSCA", "infoFeature", "rankingFeature.h", "infoSubclass", "incidenceMatrix",
            "Biplot.coordinates", "samplesSubclass", "var",
            "clus.method", "hclustSamp", "hclustFeat",
            "inner", "dispersion"
        )
    } else {
        results <- list(
            ca_res, info, rankingH, infoSubclass, mx, ca_ev, samplesSubclass, ca_res$Inertia[nd, 3],
            method.dend, info.dend.samp, info.dend.feat
        )
        names(results) <- c(
            "NSCA", "infoFeature", "rankingFeature.h", "infoSubclass", "incidenceMatrix",
            "Biplot.coordinates", "samplesSubclass", "var",
            "clus.method", "hclustSamp", "hclustFeat"
        )
    }

    names(results$var) <- "Variability explained by NSCA"
    message(.timestamp(), "-- Summary of clustering analysis:")

    print(infoSubclass)
    print(results$var)

    return(results)
}


##############################
##############################
###

.cluster.stat <- function(d = NULL, clustering, alt.clustering = NULL, noisecluster = FALSE,
                          silhouette = TRUE, G2 = FALSE, G3 = FALSE, wgap = TRUE, sepindex = TRUE,
                          sepprob = 0.1, sepwithnoise = TRUE, compareonly = FALSE,
                          aggregateonly = FALSE) {
    if (!is.null(d)) {
        d <- as.dist(d)
    }
    cn <- max(clustering)
    clusteringf <- as.factor(clustering)
    clusteringl <- levels(clusteringf)
    cnn <- length(clusteringl)
    if (cn != cnn) {
        warning("clustering renumbered because maximum != number of clusters")
        for (i in seq_len(cnn)) clustering[clusteringf == clusteringl[i]] <- i
        cn <- cnn
    }
    n <- length(clustering)
    noisen <- 0
    cwn <- cn
    if (noisecluster) {
        noisen <- sum(clustering == cn)
        cwn <- cn - 1
    }
    diameter <- average.distance <- median.distance <- separation <- average.toother <- cluster.size <-
        within.dist <- between.dist <- numeric(0)
    for (i in seq_len(cn)) cluster.size[i] <- sum(clustering == i)
    pk1 <- cluster.size / n
    pk10 <- pk1[pk1 > 0]
    h1 <- -sum(pk10 * log(pk10))
    corrected.rand <- vi <- NULL
    if (!is.null(alt.clustering)) {
        choose2 <- function(v) {
            out <- numeric(0)
            for (i in seq_len(length(v))) out[i] <- ifelse(v[i] >= 2,
                    choose(v[i], 2), 0
                )
            out
        }
        cn2 <- max(alt.clustering)
        clusteringf <- as.factor(alt.clustering)
        clusteringl <- levels(clusteringf)
        cnn2 <- length(clusteringl)
        if (cn2 != cnn2) {
            warning("alt.clustering renumbered because maximum != number of clusters")
            for (i in seq_len(cnn2)) alt.clustering[clusteringf == clusteringl[i]] <- i
            cn2 <- cnn2
        }
        nij <- table(clustering, alt.clustering)
        dsum <- sum(choose2(nij))
        cs2 <- numeric(0)
        for (i in seq_len(cn2)) cs2[i] <- sum(alt.clustering == i)
        sum1 <- sum(choose2(cluster.size))
        sum2 <- sum(choose2(cs2))
        pk2 <- cs2 / n
        pk12 <- nij / n
        corrected.rand <- (dsum - sum1 * sum2 / choose2(n)) / ((sum1 +
            sum2) / 2 - sum1 * sum2 / choose2(n))
        pk20 <- pk2[pk2 > 0]
        h2 <- -sum(pk20 * log(pk20))
        icc <- 0
        for (i in seq_len(cn)) for (j in seq_len(cn2)) if (pk12[i, j] > 0) {
                    icc <- icc + pk12[i, j] * log(pk12[i, j] / (pk1[i] *
                        pk2[j]))
                }
        vi <- h1 + h2 - 2 * icc
    }
    if (compareonly) {
        out <- list(corrected.rand = corrected.rand, vi = vi)
    }
    else {
        dmat <- as.matrix(d)
        within.cluster.ss <- 0
        overall.ss <- nonnoise.ss <- sum(d^2) / n
        if (noisecluster) {
            nonnoise.ss <- sum(as.dist(dmat[
                clustering <= cwn,
                clustering <= cwn
            ])^2) / sum(clustering <= cwn)
        }
        ave.between.matrix <- separation.matrix <- matrix(0,
            ncol = cn, nrow = cn
        )
        di <- list()
        for (i in seq_len(cn)) {
            cluster.size[i] <- sum(clustering == i)
            di <- as.dist(dmat[clustering == i, clustering ==
                i])
            if (i <= cwn) {
                within.cluster.ss <- within.cluster.ss + sum(di^2) / cluster.size[i]
                within.dist <- c(within.dist, di)
            }
            if (length(di) > 0) {
                diameter[i] <- max(di)
            } else {
                diameter[i] <- NA
            }
            average.distance[i] <- mean(di)
            median.distance[i] <- median(di)
            bv <- numeric(0)
            for (j in seq_len(cn)) {
                if (j != i) {
                    sij <- dmat[clustering == i, clustering ==
                        j]
                    bv <- c(bv, sij)
                    if (i < j) {
                        separation.matrix[i, j] <- separation.matrix[
                            j,
                            i
                        ] <- min(sij)
                        ave.between.matrix[i, j] <- ave.between.matrix[
                            j,
                            i
                        ] <- mean(sij)
                        if (i <= cwn & j <= cwn) {
                            between.dist <- c(between.dist, sij)
                        }
                    }
                }
            }
            separation[i] <- min(bv)
            average.toother[i] <- mean(bv)
        }
        average.between <- mean(between.dist)
        average.within <- mean(within.dist)
        nwithin <- length(within.dist)
        nbetween <- length(between.dist)
        between.cluster.ss <- nonnoise.ss - within.cluster.ss
        ch <- between.cluster.ss * (n - noisen - cwn) / (within.cluster.ss *
            (cwn - 1))
        clus.avg.widths <- avg.width <- NULL
        if (silhouette) {
            sii <- silhouette(clustering, dmatrix = dmat)
            sc <- summary(sii)
            clus.avg.widths <- sc$clus.avg.widths
            if (noisecluster) {
                avg.width <- mean(sii[clustering <= cwn, 3])
            } else {
                avg.width <- sc$avg.width
            }
        }
        g2 <- g3 <- cn2 <- cwidegap <- widestgap <- sindex <- NULL
        if (G2) {
            splus <- sminus <- 0
            for (i in seq_len(nwithin)) {
                splus <- splus + sum(within.dist[i] < between.dist)
                sminus <- sminus + sum(within.dist[i] > between.dist)
            }
            g2 <- (splus - sminus) / (splus + sminus)
        }
        if (G3) {
            sdist <- sort(c(within.dist, between.dist))
            sr <- nwithin + nbetween
            dmin <- sum(sdist[seq_len(nwithin)])
            dmax <- sum(sdist[(sr - nwithin + 1):sr])
            g3 <- (sum(within.dist) - dmin) / (dmax - dmin)
        }
        pearsongamma <- cor(c(within.dist, between.dist), c(rep(
            0,
            nwithin
        ), rep(1, nbetween)))
        dunn <- min(separation[seq_len(cwn)]) / max(diameter[seq_len(cwn)], na.rm = TRUE)
        acwn <- ave.between.matrix[seq_len(cwn), seq_len(cwn)]
        dunn2 <- min(acwn[upper.tri(acwn)]) / max(average.distance[seq_len(cwn)],
            na.rm = TRUE
        )
        if (wgap) {
            cwidegap <- rep(0, cwn)
            for (i in seq_len(cwn)) if (sum(clustering == i) > 1) {
                    cwidegap[i] <- max(hclust(as.dist(dmat[clustering ==
                        i, clustering == i]), method = "single")$height)
                }
            widestgap <- max(cwidegap)
        }
        if (sepindex) {
            psep <- rep(NA, n)
            if (sepwithnoise | !noisecluster) {
                for (i in seq_len(n)) psep[i] <- min(dmat[i, clustering !=
                        clustering[i]])
                minsep <- floor(n * sepprob)
            }
            else {
                dmatnn <- dmat[clustering <= cwn, clustering <=
                    cwn]
                clusteringnn <- clustering[clustering <= cwn]
                for (i in seq_len((n - noisen))) psep[i] <- min(dmatnn[
                        i,
                        clusteringnn != clusteringnn[i]
                    ])
                minsep <- floor((n - noisen) * sepprob)
            }
            sindex <- mean(sort(psep)[seq_len(minsep)])
        }
        if (!aggregateonly) {
            out <- list(
                n = n, cluster.number = cn, cluster.size = cluster.size,
                min.cluster.size = min(cluster.size[seq_len(cwn)]),
                noisen = noisen, diameter = diameter, average.distance = average.distance,
                median.distance = median.distance, separation = separation,
                average.toother = average.toother, separation.matrix = separation.matrix,
                ave.between.matrix = ave.between.matrix, average.between = average.between,
                average.within = average.within, n.between = nbetween,
                n.within = nwithin, max.diameter = max(diameter[seq_len(cwn)],
                    na.rm = TRUE
                ), min.separation = sepwithnoise *
                    min(separation) + (!sepwithnoise) * min(separation[seq_len(cwn)]),
                within.cluster.ss = within.cluster.ss, clus.avg.silwidths = clus.avg.widths,
                avg.silwidth = avg.width, g2 = g2, g3 = g3, pearsongamma = pearsongamma,
                dunn = dunn, dunn2 = dunn2, entropy = h1, wb.ratio = average.within / average.between,
                ch = ch, cwidegap = cwidegap, widestgap = widestgap,
                sindex = sindex, corrected.rand = corrected.rand,
                vi = vi
            )
        } else {
            out <- list(
                n = n, cluster.number = cn, min.cluster.size = min(cluster.size[seq_len(cwn)]),
                noisen = noisen, average.between = average.between,
                average.within = average.within, max.diameter = max(diameter[seq_len(cwn)],
                    na.rm = TRUE
                ), min.separation = sepwithnoise *
                    min(separation) + (!sepwithnoise) * min(separation[seq_len(cwn)]),
                ave.within.cluster.ss = within.cluster.ss / (n - noisen),
                avg.silwidth = avg.width, g2 = g2, g3 = g3, pearsongamma = pearsongamma,
                dunn = dunn, dunn2 = dunn2, entropy = h1, wb.ratio = average.within / average.between,
                ch = ch, widestgap = widestgap, sindex = sindex,
                corrected.rand = corrected.rand, vi = vi
            )
        }
    }
    out
}

##############################
##############################
###

NSCA <- function(N, v = 80) {
    I <- nrow(N) # Number of rows of table
    J <- ncol(N) # Number of columns of table
    Inames <- dimnames(N)[1] # Row category names
    Jnames <- dimnames(N)[2] # Column category names
    n <- sum(N) # Total number of classifications in the table
    p <- N * (1 / n) # Matrix of joint relative proportions
    Imass <- as.matrix(rowSums(p))
    Jmass <- as.matrix(colSums(p))
    ItJ <- Imass %*% t(Jmass)
    y <- p - ItJ
    dI <- diag(c(Imass), nrow = I, ncol = I)
    dJ <- diag(c(Jmass), nrow = J, ncol = J)
    Ih <- Imass^-0.5
    Jh <- Jmass^-0.5
    dIh <- diag(c(Ih), nrow = I, ncol = I)
    dJh <- diag(c(Jh), nrow = J, ncol = J)
    uni <- rep(1, times = J) # A unitary vector of length J
    x <- (p %*% solve(dJ) - Imass %*% t(uni)) %*% sqrt(dJ)
    # Column Weighted Matrix of Pi
    sva <- svd(x) # SVD of Pi
    dmu <- diag(sva$d)

    fbip <- sva$u
    dimnames(fbip) <- list(paste(Inames[[1]]), paste(seq_len(min(I, J))))


    f <- sva$u %*% dmu # Row Principal Coordinates
    g <- dJh %*% sva$v %*% dmu # Column Principal Coordinates
    dimnames(f) <- list(paste(Inames[[1]]), paste(seq_len(min(I, J))))
    dimnames(g) <- list(paste(Jnames[[1]]), paste(seq_len(min(I, J))))


    Principal.Inertia <- diag(t(f[, seq_len(min(I - 1, J - 1))]) %*% f[, seq_len(min(I - 1, J - 1))])
    Total.Inertia <- sum(Principal.Inertia)
    tau.num <- Total.Inertia # Numerator of Goodman-Kruskal tau
    tau.denom <- 1 - sum(Imass^2) # Denominator of Goodman-Kruskal tau
    tau <- tau.num / tau.denom # Goodman-Kruskal tau index
    Ctau <- (n - 1) * (I - 1) * tau # Light & Margolin's C-statistic
    Percentage.Inertia <- 100 * (Principal.Inertia / tau.num)
    Cumm.Inertia <- cumsum(Percentage.Inertia)
    Inertia <- cbind(Principal.Inertia, Percentage.Inertia, Cumm.Inertia)
    dimnames(Inertia)[1] <- list(paste("Axis", seq_len(min(I - 1, J - 1)), sep = " "))
    q.value <- 1 - pchisq(Ctau, df = (I - 1) * (J - 1))
    inner.prod <- fbip[, seq_len(which(Cumm.Inertia >= v)[1])] %*%
        t(g[, seq_len(which(Cumm.Inertia >= v)[1])])
    rownames(inner.prod) <- rownames(fbip)
    tau.num.j <- rowSums(apply(g^2, 2, function(x) Jmass * (x) / tau.num))
    names(tau.num.j) <- rownames(Jmass)

    list(
        N = N, f = f, fbip = fbip, g = g, tau = tau, tau.num.j = tau.num.j,
        di = apply(f, 1, function(x) sum(x^2)), dj = apply(g, 1, function(x) sum(x^2)),
        Cstat = Ctau, Total.Inertia = Total.Inertia, P.Value = q.value,
        Inertia = Inertia, inner.prod = inner.prod
    )
}

######################################################
## Short functions

range01 <- function(x) {
    (x - min(x)) / (max(x) - min(x))
}

is.even <- function(x) {
    x %% 2 == 0
}
