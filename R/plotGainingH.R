plotGainingH <- function(deco, data, ids, print.annot = FALSE,
                         orig.classes = TRUE) {
  
    if (all(!ids %in% rownames(deco@featureTable))) {
        stop("ERROR: Any provided IDs were found as significant by DECO.")
    }

    if (any(ids %in% rownames(deco@featureTable))) {
        msg1 <- c("NOTE: Following IDs have not been found by DECO, so they 
              will be not plotted:")
        msg2 <- paste(" ", ids[!(ids %in% rownames(deco@featureTable))], "\n ")
    }

    message(paste(msg1, msg2, sep = "\n"))

    if (print.annot & "SYMBOL" %in% colnames(featureTable(deco))) {
        msg <- "NOTE: IDs have been mapped to HGNC symbol."
        name <- deco@featureTable[ids[ids %in% rownames(deco@featureTable)], "SYMBOL"]
        names(name) <- ids[ids %in% rownames(deco@featureTable)]
    } else {
        msg <- "NOTE: Original feature IDs will be used."
        name <- ids[ids %in% rownames(deco@featureTable)]
        names(name) <- name
    }

    message(msg)

    ## Type of analysis run
    analysis <- .analysisType(deco)

    if (analysis == "Unsupervised" & orig.classes) {
        orig.classes <- FALSE
        message("NOTE: Unsupervised analysis is not associated to original classes,
            DECO subclasses will be represented.")
    }

    if (orig.classes) {
        classes <- deco@classes
    } else {
        if (analysis == "Binary") {
            classes <- c(
                deco@NSCAcluster[[1]]$samplesSubclass[, 1],
                deco@NSCAcluster[[2]]$samplesSubclass[, 1]
            )
        } else {
            classes <- deco@NSCAcluster[[1]]$samplesSubclass[, 1]
        }
    }

    ## Color
    color <- jColor(cbind(sort(classes)))

    ## Plotting
    for (id in ids[ids %in% rownames(deco@featureTable)]) {
        if (analysis == "Binary") {
            y <- c(
                deco@NSCAcluster[[1]]$NSCA$h[id, ],
                deco@NSCAcluster[[2]]$NSCA$h[id, ]
            )
        } else {
            y <- deco@NSCAcluster[[1]]$NSCA$h[id, ]
        }

        x <- data[id, names(y)]

        d <- data.frame(x, y)
        colnames(d) <- c("omic.data", "h.statistic")

        plot1 <- suppressWarnings(ggplot(
            data.frame(d, classes = classes[rownames(d)]),
            aes(y = omic.data, x = h.statistic, color = classes)
        ) + geom_point() +
            ggtitle(paste(name[id], "raw statistics")) +
            xlab("h-statistic") + ylab("omic data") +
            geom_smooth(method = "loess", span = 1, se = FALSE) + theme_minimal() +
            scale_color_manual(values = adjustcolor(names(color$ty), 0.7)))

        d2 <- data.frame(apply(d, 2, function(x)
            rank(x, ties.method = "random")), classes = classes[rownames(d)])

        plot2 <- suppressWarnings(ggplot(d2, aes(y = omic.data, x = h.statistic, color = classes)) +
            geom_point() + ggtitle(paste(name[id], "rankings")) +
            xlab("ranking h-statistic") + ylab("ranking omic data") +
            geom_smooth(method = "loess", span = 1, se = FALSE) + theme_minimal() +
            scale_color_manual(values = adjustcolor(names(color$ty), 0.7)))

        abc <- by(d2[, seq_len(2)], classes[rownames(d2)], data.frame)
        df <- lapply(names(abc), function(x) cbind(name = x, abc[[x]]))
        df <- reshape2::melt(do.call(rbind, df), id = "name")
        colnames(df) <- c("classes", "stat", "ranking")

        plot3 <- ggplot(df, aes(x = stat, y = ranking, color = classes)) +
            geom_boxplot(outlier.alpha = 0.6, col = rep(names(color$ty), each = 2), lwd = 0.5) +
            facet_grid(~classes) + theme_minimal()

        gridExtra::grid.arrange(plot1, plot2, plot3,
            layout_matrix = matrix(c(1, 2, 3, 3), byrow = TRUE, nrow = 2)
        )
    }
}
