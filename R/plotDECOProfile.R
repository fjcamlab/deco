##################################################################################################################################
##################################################################################################################################

plotDECOProfile <- function(deco, id, data, pdf.file = NA, plot.h = FALSE,
                            info.sample = NA, print.annot = TRUE, 
                            cex.legend = 1.1, cex.names = 1,
                            cex.samples = 1) {
  # Open pdf file connection.
  if (is.na(pdf.file)) {
    stop("ERROR: No pdf name has been provided.")
  }
  pdf(file = pdf.file, width = 18, height = 12, family = "Helvetica-Narrow")
  
  # Type of deco analysis
  if (all(!id %in% as.character(deco@featureTable[, c("ID")]))) {
    stop("ERROR: Input IDs do not match with any DE feature.")
  }
  if (all(is.na(deco@classes))) {
    analysis <- "Unsupervised"
  } else if (length(deco@NSCAcluster) == 1) {
    analysis <- "Multiclass"
  } else {
    analysis <- "Binary"
  }
  
  # Formatting objects...
  if (print.annot) {
    symbol <- deco@featureTable[, "SYMBOL"]
    names(symbol) <- deco@featureTable[, "ID"]
  }
  layout(mat = matrix(c(
    1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4, 5, 5,
    6, 6, 6, 0
  ), nrow = 2, byrow = TRUE))
  id_sus <- rownames(deco@NSCAcluster[[1]]$NSCA$h)
  limy <- range(na.omit(data))
  id0 <- as.character(id[id %in% id_sus])
  main <- paste("Raw data profile of:", id0, sep = " ")
  if (any(!(id %in% id_sus))) {
    message("NOTE: Following IDs have not been found by DECO or they have not been included in NSCA
            , so 'Profile' will be not plotted for them:")
    message(paste(id[!(id %in% id_sus)], "\n"))
  }
  data0 <- data
  # This plot will be repeated along all IDs input.
  for (a in seq_len(length(id0))) {
    par(mar = c(10, 6, 10, 6))
    data <- data0
    # Taking 'h' statistic for each ID and then clustering samples using this
    # information.
    if (analysis == "Binary") {
      # Sample membership to subclasses.
      samplesSubclass <- rbind(cbind(deco@NSCAcluster$Control$samplesSubclass[order(deco@NSCAcluster$Control$samplesSubclass[
        ,
        1
        ]), ]), cbind(deco@NSCAcluster$Case$samplesSubclass[order(deco@NSCAcluster$Case$samplesSubclass[
          ,
          1
          ]), ]))
      samplesSubclass <- cbind(
        Samples = rownames(samplesSubclass),
        Subclass = samplesSubclass[, 1]
      )
      rownames(samplesSubclass) <- seq_len(length(rownames(samplesSubclass)))
      # Subclasses
      infoSubclass <- rbind(deco@NSCAcluster$Control$infoSubclass, 
                            deco@NSCAcluster$Case$infoSubclass)
      # Assignment of colors to all subclasses. It will be conserved along all
      # the report.
      count <- table(vapply(rownames(infoSubclass), function(x) unlist(strsplit(x,
                                                                                split = " Subclass"
      ))[1], character(1)) != deco@control)
      color.cluster <- c(
        colorRampPalette(c("navyblue", "lightblue"))(count[1]),
        colorRampPalette(c("orange", "darkred"))(count[2])
      )
      names(color.cluster) <- rownames(infoSubclass)
      # h-statistic vector
      rangeH <- range(c(deco@NSCAcluster$Control$NSCA$h), c(deco@NSCAcluster$Case$NSCA$h))
      tab.info <- c(deco@NSCAcluster$Control$NSCA$h[as.character(id0)[a], ], deco@NSCAcluster$Case$NSCA$h[as.character(id0)[a], ])
      tab.info[is.na(tab.info)] <- 0
      # Clustering NSCA inner product for one feature.
      patt <- innerProductAssign(
        inner = tab.info, samples = deco@classes,
        control = deco@control, analysis
      )$cl
      count <- table(patt)
      thr <- which(cumsum(count) >= table(deco@classes)[names(table(deco@classes)) ==
                                                          deco@control])[1]
      colors <- c(colorRampPalette(c("blue", "lightblue"))(thr), colorRampPalette(c(
        "red",
        "darkorange"
      ))(length(count) - thr))
      names(count) <- c(
        paste("Pattern", seq_len(thr), names(table(deco@classes))[1]),
        paste("Pattern", seq_len(length(count) - thr), names(table(deco@classes))[2])
      )
      ord <- names(patt)
      bg.col <- adjustcolor(c(rep("navyblue", table(deco@classes)[names(table(deco@classes)) ==
                                                                    deco@control]), rep("darkred", table(deco@classes)[names(table(deco@classes)) !=
                                                                                                                         deco@control])), 0.6)
      limb <- list(range(deco@NSCAcluster$Control$rankingFeature.h[
        ,
        seq(
          2, dim(deco@NSCAcluster$Control$infoSubclass)[1] * 2,
          2
        )
        ]), range(deco@NSCAcluster$Case$rankingFeature.h[, seq(
          2,
          dim(deco@NSCAcluster$Case$infoSubclass)[1] * 2, 2
        )]))
      limb <- limb[[which(unlist(lapply(limb, diff)) == max(unlist(lapply(
        limb,
        diff
      ))))]]
      pch <- c(rep(15, table(bg.col)[1]), rep(16, table(bg.col)[2]))
    } else {
      # Sample membership to subclasses.
      samplesSubclass <- as.matrix(deco@NSCAcluster$All$samplesSubclass[order(deco@NSCAcluster$All$samplesSubclass[
        ,
        1
        ]), ])
      samplesSubclass <- cbind(
        Samples = rownames(samplesSubclass),
        Subclass = samplesSubclass[, 1]
      )
      rownames(samplesSubclass) <- seq_len(length(rownames(samplesSubclass)))
      # Subclasses
      infoSubclass <- deco@NSCAcluster$All$infoSubclass
      # Assignment of colors to all subclasses. It will be conserved along all
      # the report.
      color.cluster <- colorRampPalette(c(
        "greenyellow", "lightblue",
        "darkgreen", "darkred"
      ))(dim(infoSubclass)[1])
      names(color.cluster) <- rownames(infoSubclass)
      # h-statistic vector
      rangeH <- range(deco@NSCAcluster$All$NSCA$h)
      tab.info <- c(deco@NSCAcluster$All$NSCA$h[as.character(id0[a]), ])
      tab.info[is.na(tab.info)] <- 0
      # Clustering NSCA inner product for one feature.
      patt <- innerProductAssign(
        inner = tab.info, samples = deco@classes,
        analysis = analysis
      )$cl
      count <- table(patt)
      names(count) <- paste("Pattern", seq_len(length(count)))
      colors <- colorRampPalette(c("black", "grey"))(length(count))
      ord <- names(patt)
      bg.col <- adjustcolor(rep("darkred", length(names(patt))), 0.6)
      limb <- range(deco@NSCAcluster$All$rankingFeature.h[, seq(2, length(color.cluster) *
                                                                  2, 2)])
      pch <- 16
    }
    # Ordering raw data using clustering information.
    data <- data[, ord]
    names(colors) <- names(count)
    count.col <- rep(colors, count[names(colors)])
    
    # Plot feature profile.
    if (!plot.h) {
      plot(data[as.character(id0[a]), ],
           ylim = limy, pch = pch, xaxt = "n",
           lwd = 0.8, xlim = c(1, ncol(data)), main = paste(
             main[a],
             "\nsorted by h-statistic ranking."
           ), type = "p", ylab = "Raw Data Signal",
           col = bg.col, xlab = "", cex = cex.samples + 0.5, axes = FALSE
      )
    } else {
      interval <- unique(c(seq(rangeH[1], 0, length.out = 10), seq(0,
                                                                   rangeH[2],
                                                                   length.out = 11
      )))
      colorH <- colorRampPalette(c(
        "blue", "lightblue", "grey", "tomato",
        "red"
      ))(20)
      bg.col <- adjustcolor(colorH[vapply(tab.info[ord], function(x) max(which(interval <=
                                                                                 x)), numeric(1))], 0.8)
      plot(data[as.character(id0[a]), ],
           ylim = limy, pch = 22, xaxt = "n",
           lwd = 1.2, xlim = c(1, ncol(data)), main = main[a], type = "p",
           ylab = "Raw Data Signal", col = adjustcolor("white", 0.6),
           bg = bg.col, xlab = "", cex = cex.samples + 0.5, axes = FALSE
      )
    }
    
    axis(2, las = 2, cex.axis = 1.5)
    axis(side = 1, at = seq_len(dim(data)[2]), labels = rownames(samplesSubclass[match(
      ord,
      samplesSubclass[, c("Samples")]
    ), ]), las = 2, cex.axis = cex.names +
      0.5)
    m <- sum(abs(limy))
    if (analysis == "Binary") {
      thr <- dim(deco@NSCAcluster$Control$samplesSubclass)[1] + 0.5
      abline(v = thr, col = "grey", lwd = 2, lty = 2)
      segments(
        x0 = c(1, thr + 0.5), x1 = c(thr - 0.5, dim(data)[2]),
        y0 = c(mean(data[as.character(id0[a]), seq_len(thr)]), mean(data[
          as.character(id0[a]),
          (thr + 1):dim(data)[2]
          ])), col = adjustcolor("black", 0.4),
        lwd = 2
      )
    } else {
      segments(x0 = 1, x1 = dim(data)[2], y0 = mean(data[as.character(id0[a]), ]), col = adjustcolor("black", 0.4), lwd = 2)
    }
    if (print.annot) {
      id_hgnc <- symbol[names(symbol) %in% as.character(id0)[a]]
      mtext(text = paste("HGNC symbol:", id_hgnc), line = 1)
    }
    par(xpd = NA)
    # Including sample information under the expression plot.
    if (all(!is.na(info.sample))) {
      infoS <- jColor(info.sample)
      cc <- infoS$col[, sort(colnames(infoS$col))]
      for (i in seq_len(dim(cc)[2])) {
        rect(0.5:(length(data[as.character(id0[a]), ]) - 0.5), xright = 1.5:(length(data[as.character(id0[a]), ]) + 0.5), border = NA, limy[1] - (m * 0.05 * i + 0.1 *
                                                                                                                                                    m), ytop = limy[1] - (m * 0.05 * (i + 1) + 0.1 * m), col = cc[
                                                                                                                                                      colnames(data),
                                                                                                                                                      i
                                                                                                                                                      ])
        text(
          x = length(data[as.character(id0[a]), ]) + 1, labels = colnames(cc)[i],
          y = limy[1] - (m * 0.05 * (i + 0.5) + 0.1 * m), adj = c(
            0,
            0.5
          )
        )
      }
    }
    par(xpd = FALSE)
    # Sample information and feature statistics will be plotted on the right
    # side.
    title(main = list("RDA & NSCA: Single feature pattern",
                      cex = 1.8,
                      font = 2
    ), outer = TRUE, line = -2)
    par(mar = c(7, 5, 4, 5))
    # 'h' statistic per Subclass of samples.
    d <- data.frame(samplesSubclass, h = tab.info[samplesSubclass[, "Samples"]])
    m <- vapply(names(color.cluster), function(x) mean(d[d[, "Subclass"] ==
                                                           x, "h"]), numeric(1))
    s <- vapply(names(color.cluster), function(x) sd(d[d[, "Subclass"] ==
                                                         x, "h"]) / sqrt(length(d[d[, "Subclass"] == x, "h"])), numeric(1))
    par(xpd = NA)
    plot(0,
         xlim = c(0, length(m)), ylim = limb, type = "n", axes = FALSE,
         xlab = "", ylab = "h mean per DECO subclass"
    )
    rect(
      xleft = 0:(length(m) - 1), xright = seq_len(length(m)) - 0.1,
      ybottom = 0, ytop = m, col = adjustcolor(color.cluster, 0.6),
      border = NA
    )
    segments(
      x0 = seq(0.45, length(m), 1), x1 = seq(0.45, length(m), 1),
      y0 = m - s, y1 = m + s
    )
    segments(x0 = seq(0.45, length(m), 1) - 0.1, x1 = seq(
      0.45, length(m),
      1
    ) + 0.1, y0 = m - s, y1 = m - s)
    segments(x0 = seq(0.45, length(m), 1) - 0.1, x1 = seq(
      0.45, length(m),
      1
    ) + 0.1, y0 = m + s, y1 = m + s)
    axis(1, labels = rownames(infoSubclass), at = seq(
      0.45, length(m),
      1
    ), lty = 0, las = 2, cex.axis = 1.3)
    axis(2, las = 2, cex.axis = 1.5)
    
    par(xpd = FALSE)
    abline(h = 0)
    mtext(
      text = "h-statistic mean per subclass", side = 1, line = 10,
      font = 2
    )
    # Color code of 'h' statistic per sample.
    if (plot.h) {
      par(mar = c(7, 10, 6, 10))
      z1 <- matrix(seq_len(20), nrow = 1)
      y1 <- interval
      image(
        x = 0:1, y1, z1, col = adjustcolor(colorH, 0.7), cex = 1.2,
        axes = FALSE, xlab = "Feature's density", ylab = "", main = "'h' statistic\npoint's color code"
      )
      lines(y = density(tab.info, bw = 1)$x, range01(density(tab.info,
                                                             bw = 1
      )$y), lwd = 4, col = "cyan")
      axis(2, at = interval, lwd = 0, lwd.ticks = 1, labels = round(
        interval,
        2
      ), font = 2, las = 2)
    } else {
      plot(0, type = "n", axes = FALSE, xlab = "", ylab = "")
    }
    # Color labels and summary of statistics per feature.
    par(xpd = TRUE)
    if (analysis == "Binary") {
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
      legend("center",
             legend = paste(c("CONTROL", "CASE"), "samples"),
             col = "grey", pt.bg = adjustcolor("black", 0.7), pch = c(
               21,
               22
             ), bty = "n", pt.cex = 2, title = "Constrast design:\npoint's shape",
             xpd = TRUE, cex = cex.legend + 0.5
      )
    } else if (analysis == "Multiclass") {
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
      legend("center",
             legend = levels(deco@classes), col = colorRampPalette(c(
               "navyblue",
               "darkorange", "darkred"
             ))(length(levels(deco@classes))), pch = 16,
             bty = "n", pt.cex = 2, title = "Constrast design:\ncolor of point",
             xpd = TRUE, cex = cex.legend + 0.5
      )
    } else {
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
      legend("center",
             legend = c("No classes were defined\nfor RDA subsampling"),
             bty = "n", pt.cex = 2, xpd = TRUE, xjust = 0.5, cex = cex.legend +
               0.5
      )
    }
    
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
    if (all(!is.na(info.sample))) {
      legend("center",
             legend = infoS$ty, col = names(infoS$ty), xpd = TRUE,
             pch = 15, cex = cex.legend + 1 / (length(infoS$ty)) + 0.5, title = "Sample info provided by user",
             bty = "n"
      )
    }
    if (analysis == "Binary") {
      textplot(rbind(round(
        t(deco@featureTable[
          as.character(id0)[a],
          colnames(deco@featureTable) %in% c(
            "Standard.Chi.Square",
            "overlap.Ctrl.Case", "Repeats", "delta.signal", "sd.Ctrl",
            "sd.Case", "h.Range.Ctrl", "h.Range.Case", "Dendrogram.group"
          )
          ]),
        3
      ), Profile = as.character(deco@featureTable[
        as.character(id0)[a],
        c("Profile")
        ]), UpDw = as.character(deco@featureTable[
          as.character(id0)[a],
          c("UpDw")
          ])), halign = "right")
    } else {
      textplot(round(t(deco@featureTable[as.character(id0)[a], c(
        "Standard.Chi.Square",
        "Repeats", "sd", "h.Range", "Dendrogram.group"
      )]), 3), halign = "right")
    }
    par(xpd = FALSE)
  }
  suppressWarnings(remove(m))
  
  dev.off()
}
