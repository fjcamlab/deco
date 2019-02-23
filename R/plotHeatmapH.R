####################################################################################################
####################################################################################################

plotHeatmapH <- function(deco,
                         info.sample = NA, info.feature = NA,
                         print.annot = FALSE,
                         cex.legend = 1, cex.names = 1) { 
    ## Errors
    if (all(!is.na(info.sample))) {
          if (!is(info.sample, "data.frame") & !is(info.sample, "matrix")) {
                stop("ERROR: 'info.sample' must be a matrix or data.frame object.")
            }
      }

    infoS <- lapply(list(sample = info.sample, feature = info.feature), function(x) {
        if (all(!is.na(x) & !is.null(x) & x != "")) {
            suppressWarnings(infoS <- jColor(x))
            pos <- TRUE
        }
        else if (length(x) > 1 & !pos) {
              stop(paste("ERROR:", names(x), "must avoid NA, NULL or empty information."))
          } else {
              infoS <- NA
          }
        return(infoS)
    })

    # Setting initial variables
    symbol <- FALSE
    if (print.annot & "SYMBOL" %in% colnames(featureTable(deco))) {
        msg1 <- "NOTE: IDs have been mapped to HGNC symbol."
        symbol <- TRUE
    } else {
          msg1 <- "NOTE: Original feature IDs will be used."
      }

    ## Type of analysis run
    analysis <- .analysisType(deco)

    if (analysis == "Binary") {
          msg2 <- "Binary analysis, two h-statistic matrices include control or case samples."
      } else {
          msg2 <- paste(analysis, "analysis, h-statistic matrix includes all samples.")
      }

    ## Prompt messages
    message(paste(msg1, msg2, sep = "\n"))

    if (analysis == "Binary") {
          user <- as.numeric(readline(
              prompt = message("\n Plot control (0) or case (1) matrix?:\n")
          ))
      }

    ## Intern objects depending on "analysis"
    if (analysis == "Binary") {
        hMatrix <- NSCAcluster(deco)[[user + 1]]$NSCA$h
        samplesSubclass <- NSCAcluster(deco)[user]$samplesSubclass[colnames(hMatrix), ]
        col <- list(c("navyblue", "lightblue"), c("orange", "darkred"))[[user + 1]]
        ##
        y1 <- max(deco@NSCAcluster[user + 1]$hclustSamp$cluster)[1]
        y2 <- max(deco@NSCAcluster[user + 1]$hclustSamp$cluster)[2]
        main <- paste(c("CONTROL (0)", "CASE (1)")[user + 1], "samples")
        info.dend <- deco@NSCAcluster[user + 1]$hclustSamp
        info.dendF <- deco@NSCAcluster[user + 1]$hclustFeat
        p.val <- deco@NSCAcluster[user + 1]$NSCA$P.Value
        ##
        feature.label <- as.vector(deco@featureTable[rownames(hMatrix), c("UpDw")])
        feature.label[feature.label == c("UP", "DOWN")[user + 1]] <- "chartreuse4"
        feature.label[feature.label == c("DOWN", "UP")[user + 1]] <- "firebrick2"
        feature.label[feature.label == "MIXED"] <- "black"
    } else {
        hMatrix <- NSCAcluster(deco)$All$NSCA$h
        samplesSubclass <- NSCAcluster(deco)$All$samplesSubclass[colnames(hMatrix), ]
        col <- c("greenyellow", "lightblue", "darkgreen", "darkred")
        ##
        y1 <- range(deco@NSCAcluster$All$hclustSamp$cluster)[1]
        y2 <- range(deco@NSCAcluster$All$hclustSamp$cluster)[2]
        main <- "ALL samples"
        info.dend <- deco@NSCAcluster$All$hclustSamp
        info.dendF <- deco@NSCAcluster$All$hclustFeat
        p.val <- deco@NSCAcluster$All$NSCA$P.Value
        feature.label <- rep(NA, dim(hMatrix)[1])
    }

    # Assignment of colors to all subclasses.
    color.cluster <- colorRampPalette(col)(length(table(samplesSubclass)))
    names(color.cluster) <- unique(samplesSubclass)

    # Labels to plot (original IDs or annotated ones).
    if (symbol) {
        labels <- as.character(featureTable(deco)[rownames(hMatrix), "SYMBOL"])
        names(labels) <- as.character(featureTable(deco)[rownames(hMatrix), "ID"])
        labels[is.na(labels)] <- names(labels)[is.na(labels)]
    } else {
        labels <- rownames(hMatrix)
        names(labels) <- rownames(hMatrix)
    }

    if (all(!is.na(infoS$sample))) {
        info.sample.color <- infoS$sample$col[colnames(hMatrix), ]
        cc <- color.cluster[samplesSubclass[colnames(hMatrix)]]
        # Final label matrix to place under samples dendrogram.
        info.sample.color <- cbind(info.sample.color, c(rep(NA, dim(hMatrix)[2])),
            DECO.Subclass = cc
        )
        colnames(info.sample.color)[seq_len(dim(info.sample)[2])] <- colnames(infoS$sample$col)
        m <- 0.15
        n <- -0.05
    } else {
        info.sample.color <- c(rep(NA, dim(hMatrix)[2]))
        info.sample.color <- color.cluster
        info.sample.color <- cbind(
            DECO.Subclass = info.sample.color,
            rep(NA, dim(hMatrix)[2])
        )
    }
    # Same procedure should be done for feature information.
    if (all(is.na(info.feature))) {
        infoFpat <- colorRampPalette(c("grey", "black"))(max(info.dendF$cluster))[info.dendF$cluster]
        if (analysis == "Binary") {
            feature.color <- cbind(Feature.pattern = infoFpat, rep(
                "white",
                length(feature.label)
            ), Exprs = feature.label)
        } else {
            feature.color <- cbind(Feature.pattern = infoFpat, rep(
                "white",
                length(feature.label)
            ))
        }
    } else {
        info.feature <- cbind(info.feature)
        info.feature.color <- info.feature
        infoFpat <- colorRampPalette(c("grey", "black"))(max(info.dendF$cluster))[info.dendF$cluster]
        for (z in seq_len(length(unlist(apply(info.feature, 2, unique))))) {
            info.feature.color[info.feature == unlist(apply(
                info.feature,
                2, unique
            ))[z]] <- as.character(rep(
                colorRampPalette(RColorBrewer::brewer.pal(
                    name = "Paired",
                    n = 10
                ))(length(unlist(apply(info.feature, 2, unique))))[z],
                length(info.feature.color[info.feature == unlist(apply(
                    info.feature,
                    2, unique
                ))[z]])
            ))
        }
        # Final label matrix to place under samples dendrogram.
        feature.color <- cbind(
            Feature.pattern = infoFpat, Exprs = feature.label,
            rep("white", length(feature.label)), Info.gene = info.feature.color
        )
    }
    rownames(feature.color) <- rownames(hMatrix)
    names(infoFpat) <- rownames(hMatrix)

    # h-statistic heatmap including sample and feature information.
    .heatplot.2(
        dataset = as.matrix(hMatrix), Colv = as.dendrogram(info.dend$dend),
        Rowv = as.dendrogram(info.dendF$dend), margins = c(10, 8), cexRow = cex.names,
        cexCol = cex.names + m, lmat = rbind(c(6, 0, 5, 0), c(
            0, 0, 5,
            0
        ), c(0, 0, 2, 0), c(4, 1, 3, 0)), lhei = c(
            0.25, 0.1, m,
            1
        ), lwid = c(0.5, 0.08, 1.5, 0.15), method = info.dend$dend$method,
        labRow = unname(labels), labRowCol = rep("black", dim(hMatrix)[1]),
        ColSideColorsSize = dim(info.sample.color)[2], RowSideColorsSize = dim(feature.color)[2],
        main = "", dualScale = FALSE, KeyValueName = "h-statistic", notecex = cex.names,
        ColSideColors = info.sample.color, RowSideColors = t(feature.color),
        col.breaks = 256, rm.out = TRUE
    )

    # Adding more information to plot.
    title(main = list(paste("h-statistic.", main), font = 2), line = 4.5)
    mtext(
        text = "Higher (absolute) values of h-statistic indicates more\nfeature-relevance for these samples.",
        side = 3, outer = TRUE, line = -3, las = 1
    )
    mtext(
        text = paste("f =", dim(hMatrix)[1], "features for heatmap"), side = 4,
        las = 3, font = 2, cex = 1.3
    )

    # Legends including feature and sample information.
    if (all(!is.na(infoS$sample))) {
        legend("topleft",
            legend = infoS$sample$ty, col = names(infoS$sample$ty), pch = 15,
            cex = cex.legend - 0.3, title = "Sample info provided by user",
            bty = "n", xpd = TRUE, inset = c(0, 0.15)
        )
    }
    legend(
        title = "DECO sample subclasses", "topright", legend = paste(
            "Subclass",
            seq_len(max(info.dend$cluster))
        ), col = color.cluster, pch = 15, cex = cex.legend - 0.1, bty = "n", xpd = TRUE,
        inset = c(0, 0)
    )
}
