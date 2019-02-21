##########################################################################################
##########################################################################################

plotAssociationH <- function(deco, info.sample) {
    analysis <- .analysisType(deco)
    nsc <- NSCAcluster(deco)

    if (analysis == "Binary") {
        hMatrix <- cbind(
            nsc[[1]]$NSCA$h,
            nsc[[2]]$NSCA$h
        )
        hRanking <- cbind(
            nsc[[1]]$rankingFeature.h,
            nsc[[2]]$rankingFeature.h
        )[rownames(hMatrix), ]
        subcl <- c(
            nsc[[1]]$samplesSubclass[na.omit(match(names(info.sample), 
                                                   rownames(nsc[[1]]$samplesSubclass))), 1],
            nsc[[2]]$samplesSubclass[na.omit(match(names(info.sample), 
                                                   rownames(nsc[[2]]$samplesSubclass))), 1]
        )
        subclU <- c(
            unique(nsc[[1]]$samplesSubclass),
            unique(nsc[[2]]$samplesSubclass)
        )
    } else {
        hMatrix <- nsc[[1]]$NSCA$h
        hRanking <- nsc[[1]]$rankingFeature.h[rownames(hMatrix), ]
        subcl <- nsc[[1]]$samplesSubclass[names(info.sample), 1]
    }

    if (all(!is.na(info.sample)) & (is(info.sample, "character") | is(info.sample, "factor"))) {
        msg1 <- c("Plotting h-statistic from 'info.sample' categories...")

        hRanking <- hRanking[, grepl("h.S", colnames(hRanking))]
        info.sample <- factor(info.sample[colnames(hMatrix)])
        hRanking2 <- t(apply(hMatrix, 1, function(x) by(x, info.sample, mean)))
    } else {
          stop("ERROR: 'info.sample' must be a named character vector or factor with information of samples.")
      }

    message(msg1)

    ## Calculating associations...
    rankings <- apply(hRanking, 1, function(x) which(abs(x) == max(abs(x))))
    tab <- table(rankings)[as.character(seq_len(length(colnames(hRanking))))]
    tab[is.na(tab)] <- 0
    if (analysis == "Binary") {
          names(tab) <- subclU
      }
    rankings <- names(tab)[rankings]

    rankings2 <- apply(hRanking2, 1, function(x) which(abs(x) == max(abs(x))))
    tab <- table(rankings2)[as.character(seq_len(length(colnames(hRanking2))))]
    tab[is.na(tab)] <- 0
    names(tab) <- colnames(hRanking2)
    rankings2 <- names(tab)[rankings2]

    # Assignment of colors to all subclasses.
    col <- colorRampPalette(RColorBrewer::brewer.pal("Spectral", n = 8))(dim(hRanking)[2])

    D <- reshape2::melt(hRanking2)
    colnames(D) <- c("feature", "Category", "h.statistic")
    D <- data.frame(D, Subclasses = rankings[D$feature])

    bp <- ggplot(D) + aes(x = Subclasses, y = h.statistic, fill = Subclasses) +
        xlab("DECO subclasses") +
        geom_abline(slope = 0, intercept = 0) +
        geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), color = adjustcolor("black", 0.5)) +
        # geom_jitter(height = 0, width = 0.1) +
        scale_fill_manual(values = adjustcolor(col, 0.7)) +
        ylab("h-statistic") +
        ggtitle("Association of h-statistic among DECO subclasses and input categories") +
        theme(panel.background = element_blank(), axis.text.x = element_blank())

    bp <- bp + facet_grid(. ~ Category)

    tab <- table(subcl)
    tab[is.na(tab)] <- 0
    if (analysis != "Binary") {
        tab <- tab[paste("Subclass", seq_len(length(colnames(hRanking))))]
        names(tab) <- paste("Subclass", seq_len(length(colnames(hRanking))))
    }
    subcl <- paste(subcl, " (n=", tab[subcl], ")", sep = "")

    freq <- reshape2::melt(as.matrix(table(DECO.subclasses = subcl, input.info = info.sample)))
    colnames(freq)[3] <- "Samples"
    freq[freq == 0] <- NA

    bp2 <- ggplot(freq, aes(x = input.info, y = DECO.subclasses)) + geom_tile(aes(fill = Samples)) +
        scale_fill_gradient(low = "white", high = "red", na.value = "white") +
        theme_light() + xlab("") + ylab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        geom_point(aes(size = Samples), show.legend = FALSE, na.rm = TRUE)

    freq <- reshape2::melt(as.matrix(table(DECO.subclasses = rankings, input.info = rankings2)))
    colnames(freq)[3] <- "Features"
    freq[freq == 0] <- NA

    bp3 <- ggplot(freq, aes(x = input.info, y = DECO.subclasses)) + geom_tile(aes(fill = Features)) +
        scale_fill_gradient(low = "white", high = "navyblue", na.value = "white") +
        theme_light() + xlab("") + ylab("") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        geom_point(aes(size = Features), show.legend = FALSE, na.rm = TRUE)


    gridExtra::grid.arrange(bp, bp2, bp3,
        layout_matrix = matrix(c(1, 1, 2, 3), byrow = TRUE, nrow = 2)
    )
}
