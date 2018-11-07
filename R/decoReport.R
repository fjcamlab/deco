################################################
################################################
# PDF report function

decoReport <- function(deco, sub, id = NA, pdf.file = "decoReport.pdf",
                       info.sample = NA, info.feature = NA,
                       cex.samples = 1.2, cex.legend = 0.9,
                       cex.names = 0.8, print.annot = FALSE)
{
  if(is.na(pdf.file))
    stop("PDF file name is empty.")

  # Graphical and warning options.
  options(warn=-1)
  default.par <- par()
  default.layout <- layout(mat = c(1))

  # Saving initial R objects.
  deco0 <- deco
  data <- sub$data
  data0 <- data

  if(!is(deco, "deco"))
    stop("ERROR: Input R object is not 'deco-class'.")
  if(!all(rownames(info.sample) %in% colnames(data)))
    stop("ERROR: Rownames from sample information table do not macth with colnames from data.")
  # Defining design of analysis previously run.
  if(all(is.na(deco@classes)))
    analysis <- "Unsupervised"
  else if(length(levels(sub$classes)) == 2)
    analysis <- "Binary"
  else
    analysis <- "Multiclass"
  message(paste(analysis,"DECO analysis will be generated."))
  # Setting up maximum number of features shown in tables.
  if(length(id) >= 50)
  {
    message(paste("NOTE: ID's input restringed to 50 first IDs."))
    id <- id[1:50]
  }
  if(dim(deco@featureTable)[1] >= 50)
    info.size <- 50
  else
    info.size <- dim(deco@featureTable)[1]

  # Open pdf file.
  pdf(file = pdf.file,width = 18,height = 12,family = "Helvetica-Narrow")

  # Variable used to annotate if print.annot == T and biological context is supplied.
  if("SYMBOL"%in%colnames(deco@featureTable) & print.annot == TRUE)
  {
    symbol <- as.character(deco@featureTable$SYMBOL)
    symbol[is.na(as.character(deco@featureTable$SYMBOL))] <-
      as.character(deco@featureTable$ID)[is.na(as.character(deco@featureTable$SYMBOL))]
    names(symbol) <- as.character(deco@featureTable$ID)
    symbol <- symbol[sort(names(symbol))]
  }else
    symbol <- NA

  ## Formatting R objects to print.
  if(analysis == "Binary"){
    # p.value from NSCA analysis.
    p.val <- c(deco@NSCAcluster$Control$NSCA$P.Value,deco@NSCAcluster$Case$NSCA$P.Value)
    # General information about subclasses found.
    infoSubclass <- rbind(deco@NSCAcluster$Control$infoSubclass,deco@NSCAcluster$Case$infoSubclass)
    rel <- c(rep(p.val[1] <= 0.05, dim(deco@NSCAcluster$Control$infoSubclass)[1]),rep(p.val[2] <= 0.05, dim(
      deco@NSCAcluster$Case$infoSubclass)[1]))
    infoSubclass <- data.frame(infoSubclass, Binary = rep(0,dim(infoSubclass)[1]), isRelevant = as.character(rel))
    infoSubclass[sapply(rownames(infoSubclass),
                        function(x)unlist(strsplit(x,split=" Subclass"))[1]) != deco@control, "Binary"] <- 1
    # Sample membership to subclasses.
    samplesSubclass <- rbind(cbind(deco@NSCAcluster$Control$samplesSubclass[order(deco@NSCAcluster$Control$samplesSubclass[,1]),]),
                             cbind(deco@NSCAcluster$Case$samplesSubclass[order(deco@NSCAcluster$Case$samplesSubclass[,1]),]))
    samplesSubclass <- cbind(Samples = rownames(samplesSubclass), Subclass = samplesSubclass[,1])
    rownames(samplesSubclass) <- 1:length(rownames(samplesSubclass))
    # NSCA and hierarchical clustering information.
    nsc <- matrix(round(c(deco@NSCAcluster$Control$var,deco@NSCAcluster$Case$var,p.val,deco@NSCAcluster$Control$hclustSamp$huber,
                        deco@NSCAcluster$Case$hclustSamp$huber),digits = 3),
                nrow = 3, byrow = TRUE)
    colnames(nsc) <- c("Control samples","Case samples")
    rownames(nsc) <- c("Variability explained by NSCA", "NSCA C-statistic p.value", "Huber's gamma")
    # Assignment of colors to all subclasses. It will be conserved along all the report.
    count <- table(sapply(rownames(infoSubclass),
                          function(x)unlist(strsplit(x,split=" Subclass"))[1]) != deco@control)
    color.cluster <- c(colorRampPalette(c("navyblue","lightblue"))(count[1]),colorRampPalette(c("orange","darkred"))(count[2]))
    names(color.cluster) <- rownames(infoSubclass)
    # Columns of 'featureTable' to be printed.
    names.col <- c("ID","SYMBOL","UpDw","Profile","overlap.Ctrl.Case","Standard.Chi.Square","Repeats",
                   "Repeats.index","delta.signal","sd.Ctrl","Dendrogram.group.Ctrl",
                   "h.Range.Ctrl","sd.Case","Dendrogram.group.Case","h.Range.Case")
    deco@featureTable[,names.col[!names.col %in% c("ID","SYMBOL","UpDw","Profile")]] <- apply(
      deco@featureTable[,names.col[!names.col %in% c("ID","SYMBOL","UpDw","Profile")]],2,as.numeric)
    textF <- deco@featureTable[,colnames(deco@featureTable)%in%names.col]
    deco0@featureTable <- deco@featureTable
    textF <- textF[,na.omit(match(names.col,colnames(textF)))]
    # Chi, Delta , Repeats & SD(inverse)
    ordG <- order(apply(apply(cbind(textF$Standard.Chi.Square, -textF$overlap.Ctrl.Case,
                              -textF$overlap.Ctrl.Case, textF$Repeats,abs(textF$delta.signal),
                              abs(textF$delta.signal),-textF$sd.Ctrl,-textF$sd.Case),2,function(x)
                                rank(-x, ties.method = "max")), 1, median))
    # Chi, Delta & SD
    ordS <- rank(-apply(apply(cbind(textF$Standard.Chi.Square, textF$Standard.Chi.Square, abs(textF$delta.signal), textF$Repeats,
                                    abs(textF$delta.signal),textF$h.Range.Ctrl,textF$h.Range.Case),2,function(x)
                                      rank(-x, ties.method = "max")), 1, median))
    # Saving initial modified 'featureTable' as R data.frame called 'textF0'.
    deco@featureTable <- deco@featureTable[ordG,]
    textF0 <- textF
    textF <- textF[ordG,]
    # Object defining samples membership to classes previously compared.
    classes <- c(unlist(strsplit(deco@NSCAcluster$Control$samplesSubclass[1],split=" Subclass"))[1],
                 unlist(strsplit(deco@NSCAcluster$Case$samplesSubclass[1],split=" Subclass"))[1])
  }else{
    # p.value from NSCA analysis.
    p.val <- deco@NSCAcluster$All$NSCA$P.Value
    # General information about subclasses found.
    infoSubclass <- deco@NSCAcluster$All$infoSubclass
    rel <- rep(p.val[1] <= 0.05, dim(deco@NSCAcluster$All$infoSubclass)[1])
    infoSubclass <- data.frame(infoSubclass, isRelevant = as.character(rel))
    # Sample membership to subclasses.
    samplesSubclass <- cbind(deco@NSCAcluster$All$samplesSubclass[order(
      deco@NSCAcluster$All$samplesSubclass[,1]),])
    samplesSubclass <- cbind(Samples = rownames(samplesSubclass), Subclass = samplesSubclass[,1])
    rownames(samplesSubclass) <- 1:length(rownames(samplesSubclass))
    # NSCA and hierarchical clustering information.
    nsc <- rbind(round(deco@NSCAcluster$All$var, 3), round(p.val, 3), round(deco@NSCAcluster$All$hclustSamp$huber, 3))
    colnames(nsc) <- c("All samples")
    rownames(nsc) <- c("Variability explained by NSCA", "p.value", "Huber's gamma")
    # Assignment of colors to all subclasses. It will be conserved along all the report.
    color.cluster <- colorRampPalette(c("greenyellow","lightblue","darkgreen","darkred"))(dim(infoSubclass)[1])
    names(color.cluster) <- rownames(infoSubclass)
    # Columns of 'featureTable' to be printed.
    names.col <- c("ID","SYMBOL","Standard.Chi.Square","Repeats","Repeats.index","Avrg.logFC","Dendrogram.group","h","sd","h.Range")
    deco@featureTable[,colnames(deco@featureTable)%in%names.col[!names.col %in% c("ID","SYMBOL")]] <- apply(
      deco@featureTable[,colnames(deco@featureTable)%in%names.col[!names.col %in% c("ID","SYMBOL")]],2,as.numeric)
    textF <- deco@featureTable[,colnames(deco@featureTable)%in%names.col]
    textF <- textF[,na.omit(match(names.col,colnames(textF)))]
    ordG <- order(apply(apply(cbind(textF$Standard.Chi.Square,textF$Standard.Chi.Square,textF$Repeats,
                                    textF$Repeats,textF$h.Range),2,function(x) rank(-x, ties.method = "max")), 1, median))
    ord <- ordG
    # No sense to make two different ranking for "Unsupervised" analysis.
    deco@featureTable <- deco@featureTable[ord,]
    textF0 <- textF
    textF <- textF[ord,][1:info.size,]
    classes <- sub$classes

  }

  ### FIRST SECTION

  # RDA information table
  j <- rbind(analysis, dim(deco@incidenceMatrix)[2], round(deco@featureTable[1,c("Repeats")]/deco@featureTable[1,c("FR.Repeats")],digits = 0),
             deco@pos.iter, round(deco@pos.iter/round(deco@featureTable[1,c("Repeats")]/deco@featureTable[1,c("FR.Repeats")],digits = 0)*100, 2),
             dim(deco@featureTable)[1], min(deco@featureTable$Repeats)-1, deco@q.val, sub$resampleSize)
  # Generating colors for all different sample information input by user.
  if(all(is.na(info.sample)))
    info.sample <- matrix(NA, nrow = ncol(deco@incidenceMatrix),
                          ncol = 2, dimnames = list(colnames(deco@incidenceMatrix),c("","")))
  if(all(!is.na(info.sample)))
    infoS <- jColor(info.sample)
  else
    infoS <- NA
  # Configuring device to print RDA and NSCA general information.
  layout(mat = matrix(c(1,2,3,4,5,5,5,5),nrow = 4, 2, byrow = FALSE), heights = c(2,0.75,1.5,2))
  # Information about RDA step.
  colnames(j) <- c("")
  rownames(j) <- c("Contrast design:", "Number of samples:", "Total iterations:",
                   "Positive DE iterations:","% positive DE iterations:",
                   "DE features:","Minimum repeats:","LIMMA q.value threshold:",
                   "RDA resampling size:")
  textplot(j, col.rownames = "darkblue", cex = 2, valign = "bottom")
  mtext(text = "RDA information", cex = 1.4, font = 2, line = -5)
  # Information about NSCA and subclasses clustering.
  color.table <- matrix(data = c(rep("black", dim(nsc)[2])), nrow = 3, ncol = dim(nsc)[2], byrow = FALSE)
  color.table[2, nsc[2,] > 0.05] <- "red"
  textplot(nsc, col.rownames = "darkblue", cex = 1.7, valign = "top", col.data = color.table)
  mtext(text = "NSCA information", cex = 1.4, font = 2, line = 2)
  # Top 10 feature ranking.
  y <- na.omit(cbind(Ranking = 1:10, textF[1:10,colnames(textF)%in%c("ID","UpDw","Profile","SYMBOL"),drop=FALSE]))
  textplot(object = y, show.rownames = FALSE, col.colnames = "darkblue", valign = "top", cex = 1.5)
  mtext(text = "Feature ranking information", cex = 1.4, line = 1, font = 2)
  # Subclasses found information.
  col <- color.cluster
  col[infoSubclass[,2] == 0] <- "grey"
  col[!as.logical(infoSubclass[,"isRelevant"])] <- "grey"
  textplot(valign = "top", object = infoSubclass,
           col.data = matrix(rep(col,dim(infoSubclass)[2]), ncol = dim(infoSubclass)[2], byrow = FALSE),
           col.rownames = col, title = "1", cex = 3/(0.25*(dim(infoSubclass)[1]+3)))
  mtext(text = "Subclass information", cex = 1.4, line = 1, font = 2)
  mtext(side = 1, text = "Non-relevant subclasses are grey-colored", font = 1, line = -2)
  # Samples membership to subclasses.
  color.table <- matrix(data = c(rep("black",2)), ncol = 2, nrow = dim(samplesSubclass)[1], byrow = FALSE)
  color.table[is.na(samplesSubclass[,c("Subclass")]), 1] <- "red"
  samplesSubclass[is.na(samplesSubclass[,c("Subclass")]),c("Subclass")] <- "NoSubclass"
  textplot(samplesSubclass, valign = "top", cmar = 2, mar = c(3,3,5,3), col.colnames = "darkblue",
           col.data = color.table, col.rownames = "darkblue")
  title(main = list(paste("Overview\nDECO analysis report"), cex = 2.5, font = 2), line = -6, outer = TRUE)

  ### SECOND SECTION

  # Guide SECTION showing color code used along the report.
  if(any(!is.na(info.sample))){
    layout(matrix(c(1,2,3), ncol = 3))
    par(xpd = TRUE)
    plot(0, xlab = "", ylab = "", axes = FALSE, type = "n")
    legend("center", legend = sort(colnames(infoS$orig)), cex = cex.legend + 0.7,
           title = "Sample categories provided by user", bty = 'n', y.intersp = 1.3)
    plot(0, xlab = "", ylab = "", axes = FALSE, type = "n")
    legend("center", legend = infoS$ty, col = names(infoS$ty),
           pch = 15, cex = cex.legend + 0.7, title = "Sample info provided by user", bty = 'n')
    plot(0, xlab = "", ylab = "", axes = FALSE, type = "n")
    legend("center", legend = names(color.cluster), col = color.cluster, y.intersp = 1.3,
           pch = 15, cex = cex.legend + 0.7, xpd = TRUE,
           bty = "n", title = "Subclasses of samples found by DECO")
    title(main = list(paste("Sample color code"), cex = 2.5, font = 2),
          line = -5, outer = TRUE)
    mtext(text = "Reference colors for samples along the PDF report:\ninfo provided and subclasses of samples found.", cex = 1.4,
          line = -10, outer = TRUE)
    par(xpd = FALSE)
  }

  ### THIRD SECTION

  layout(default.layout)
  .plotEffectSize(deco, infoS, samplesSubclass, analysis)

  ### FOURTH SECTION

  if(analysis == "Binary"){
    par(las = 1)
    .plotTau(x = list(deco@NSCAcluster$Control, deco@NSCAcluster$Case),
            samplesSubclass = samplesSubclass,
            color.cluster, labels = symbol, n = 20)}
  else
    .plotTau(x = list(deco@NSCAcluster$All),
            samplesSubclass = samplesSubclass,
            color.cluster, labels = symbol, n = 20)

  ### FIFTH SECTION: Heterogeneity.

  i <- 1
  par(mfrow = c(1,2), mar = c(10,6,10,6))

  while(i %in% 1:2)
  {
    if(analysis != "Binary"){par(mfrow = c(1,1),mar = c(12,8,12,8));v <- deco@NSCAcluster$All$NSCA$Inertia[,3];
    i <- 3; main = "All dataset"
    var <- deco@NSCAcluster$All$var}
    if(i == 1){v <- deco@NSCAcluster$Control$NSCA$Inertia[,3]
    main = "Control samples"
    var <- deco@NSCAcluster$Control$var}
    if(i == 2){v <- deco@NSCAcluster$Case$NSCA$Inertia[,3]
    main = "Case samples"
    var <- deco@NSCAcluster$Case$var}

    .plotHeterogeneityDECO(v = v, main = main, var = var)
    i = i+1
  }
  par(mfrow = c(1,1), xpd=NA)
  legend("bottom",legend = c("Cumulative"),col = c("darkgoldenrod"),xpd = TRUE, inset = c(0,-0.1),
         pch = 16, lty = 1, lwd = 3, cex = 1.3, bty = "n")
  title(main = list("NSCA: Smooth-curve of variability explained"), line = -2.5, outer = TRUE)

  ### SIXTH SECTION: Rankings

  layout(default.layout)
  par(mar = c(7,7,7,7), xpd = FALSE)

  if(analysis %in% c("Unsupervised","Multiclass"))
  {
    color.table <- matrix(data = c(rep("black", dim(textF)[2])),nrow = min(info.size,dim(textF)[1]),
                          ncol = dim(textF)[2]+1, byrow = FALSE)
    colnames(color.table) <- c("Ranking",colnames(textF))
    if(all(is.na(id)))
      id <- as.character(textF$ID[1:15])
    color.table[as.character(textF[1:info.size,"ID"]) == id,"ID"] <- "darkred"
    textF <- textF[textF$ID != "",]
    textF[1:info.size,!colnames(textF)%in%c("SYMBOL","ID","Tau.feature")] <-
      round(textF[1:info.size, !colnames(textF)%in%c("SYMBOL","ID","Tau.feature")],4)
    textplot(cbind(Ranking = 1:info.size, textF[1:info.size,na.omit(match(names.col,colnames(textF)))]),
             valign = "top", cmar = 2, cex = 0.95,
             mar = c(1,1,2,1), show.rownames = FALSE, col.colnames = c(rep("deepskyblue4",dim(textF)[2]-2),rep("darkred",3)),
             col.data = color.table[1:info.size,])
    title(list(paste("RDA: Top 50 feature signature."),
               cex = 1.5), outer = TRUE, line = -2.5)
    legend("bottomright",legend = c("RDA info","NSCA info"), title = "Colnames color",
           col = c("deepskyblue4","darkred"), xpd = NA, bty = "n",
           pch = 65, cex = 1.3, inset = 0.05, pt.cex = 1.5,
           pt.lwd = 3, text.font = 2)
  }else
  {
    IDS <- list()
    for(i in 1:2)
    {
      if(i == 1){
        rank <- c("Majority","Complete")
        ord <- ordG}
      else{
        rank <- "Minority"
        ord <- ordS}
      textF <- textF0[ord,]
      textF <- textF[which(textF$Profile %in% rank)[1:info.size],]
      n <- c("sd.Ctrl","Dendrogram.group.Ctrl","h.Range.Ctrl",
             "sd.Case","Dendrogram.group.Case","h.Range.Case")
      textF[,c("overlap.Ctrl.Case","Standard.Chi.Square","Repeats.index",n)] <- round(
        textF[,c("overlap.Ctrl.Case","Standard.Chi.Square","Repeats.index",n)],3)
      textF[1:info.size,colnames(textF)%in%names.col][is.na(textF[1:info.size,colnames(textF)%in%names.col])] <- "NotAssigned"
      textF <- na.omit(textF[as.character(textF$ID) != "",])
      if(dim(textF)[1] == 0){
        message("NOTE: Top50 of ",c("Complete-Majority","Minority-Mixed")[i]," was no printed, there are no feature showing this profile.")
        next
      }
      tex <- cbind(Ranking = 1:dim(textF)[1],textF[,na.omit(match(names.col,colnames(textF)))])
      # Saving both top-50 features.
      IDS[[i]] <- na.omit(as.character(textF[which(textF$Profile %in% rank)[1:dim(tex)[1]],c("ID")]))
      color.table <- matrix(data = c(rep("black",length(which(colnames(deco@featureTable)%in%names.col)))),
                            nrow = dim(tex)[1], ncol = length(which(colnames(deco@featureTable)%in%names.col)), byrow = FALSE)
      colnames(color.table) <- colnames(textF)
      if(all(is.na(id)))
        id <- na.omit(IDS[[i]][1:15])
      # Color of table.
      color.table[as.character(textF[1:dim(tex)[1],c("ID")]) %in% id,1] <- "red"
      color.table[textF[1:dim(tex)[1],colnames(textF)%in%names.col] == "MIXED"] <- "gray20"
      color.table[textF[1:dim(tex)[1],colnames(textF)%in%names.col] == "UP"] <- "darkred"
      color.table[textF[1:dim(tex)[1],colnames(textF)%in%names.col] == "DOWN"] <- "darkolivegreen"
      color.table[textF[1:dim(tex)[1],colnames(textF)%in%names.col] == "Minority"] <- "goldenrod1"
      color.table[textF[1:dim(tex)[1],colnames(textF)%in%names.col] == "Majority"] <- "chocolate2"
      color.table[textF[1:dim(tex)[1],colnames(textF)%in%names.col] == "Complete"] <- "brown4"
      color.table[,c("delta.signal")] <- color.table[,c("UpDw")]
      textplot(tex, valign = "top", cex = 0.95,
               cmar = 1.1, mar = c(7,0,4,0), show.rownames = FALSE,
               col.colnames = c(rep("black",dim(tex)[2]-6),rep("navyblue",3),rep("darkred",3)),
               col.data = cbind(rep("black",info.size),color.table))
      title(list(paste("RDA: Top",table(textF$ID != "")[1],"feature signature,",c("MAJORITIES","MINORITIES")[i],
                       "ranking method was applied"),
                 cex = 1.5), outer = TRUE, line = -2.5)
      legend("bottom",legend = c("RDA info","NSCA: CONTROL info","NSCA: CASE info"),
             col = c("black","navyblue","darkred"), xpd = NA, bty = "n",
             pch = 88, cex = 1.3, inset = -0.05, pt.cex = 1.5, horiz = TRUE, yjust = 0,
             pt.lwd = 4, text.font = 2)
    }
  }

  ### SEVENTH SECTION: Repeats threshold
  .plotRepThr(sub = sub, deco = deco, id = NA, print.annot = print.annot)

  ### EIGHTH SECTION: Overlap & dendrograms
  layout(matrix(1))
  if(analysis == "Binary"){
    .plotOD(deco0, id, ord = list(ordG, ordS), symbol = symbol, print.annot = print.annot)
    .plotDend(list(deco@NSCAcluster$Control$hclustSamp, deco@NSCAcluster$Case$hclustSamp), analysis,
             color.cluster, samplesSubclass, cex.samples+0.1)
  }else
    .plotDend(deco@NSCAcluster$All$hclustSamp, analysis, color.cluster, samplesSubclass,
             cex.names = cex.samples+0.1)

  ### NINTH SECTION: Most relevant features per subclass.
  .plotRankingH(deco = deco, analysis)

  ### TENTH SECTION: Heatmaps of 'h' statistic and original data.
  .plotHeatmapDECO(deco, p.val, id, analysis, symbol, infoS, color.cluster, data0,
                   cex.legend, cex.names, cex.samples, print.annot,
                   info.sample, info.feature, samplesSubclass)

  ### ELEVENTH SECTION: Feature profile with original data ranked by 'h' values.
  .plotProfileDECO(deco, id, data, analysis, infoS, color.cluster,
                   info.sample, infoSubclass, samplesSubclass,
                   cex.legend, cex.names, cex.samples, print.annot, symbol)

  ### TWELVETH SECTION: Biplot from NSCA coordinates.

  if(all(c("CHR","CHRLOC","SYMBOL") %in% colnames(deco@featureTable)))
    .plotDECOmanhattan(deco, id)
  else
    message(paste("There are not any CHR and CHRLOC information within 'featureTable' slot,\nManhattan plot will be not plotted."))


  ### THIRTENTH SECTION: Biplot from NSCA coordinates.

  .plotDECObiplot(deco, color.cluster, cex.samples)

  ### FOURTENTH SECTION: Feature gaining

  # if(analysis == "Binary"){
  #   par(mfrow = c(1,2), mar = c(11,6,9,6))
  #   for(i in 1:2)
  #     .plotFeatureGaining(deco@featureTable[,paste("h.Range",c("Ctrl","Case")[i],sep=".")],
  #                         deco@featureTable[,paste("sd",c("Ctrl","Case")[i],sep=".")],
  #                         deco, print.annot, main = c("CONTROL SAMPLES","CASE SAMPLES")[i], id)
  # }else{
  #   par(mfrow = c(1,1), mar = c(11,10,9,10))
  #   .plotFeatureGaining(deco@featureTable[,"h.Range"], deco@featureTable[,"sd"],
  #                       deco, print.annot, main = "ALL SAMPLES", id)
  # }

  ### Generating .txt tables
  if(analysis == "Binary")
    rankingH <- cbind(deco@NSCAcluster[[1]]$rankingFeature.h,
                      deco@NSCAcluster[[2]]$rankingFeature.h[rownames(deco@NSCAcluster[[1]]$rankingFeature.h),])
  else
    rankingH <- deco@NSCAcluster[[1]]$rankingFeature.h

  # txt.file <- rev(unlist(strsplit(pdf.file, split = "/", fixed = T)))[1]
  txt.fileF <- paste(unlist(strsplit(pdf.file, split = ".pdf", fixed = TRUE))[1], "_DECO_feature_statistics",".tsv", sep = "")
  txt.fileH <- paste(unlist(strsplit(pdf.file, split = ".pdf", fixed = TRUE))[1], "_DECO_feature_h", ".tsv", sep = "")
  txt.fileS <- paste(unlist(strsplit(pdf.file, split = ".pdf", fixed = TRUE))[1], "_DECO_samples_subclass_membership", ".tsv", sep = "")
  write.table(deco@featureTable, txt.fileF, sep = "\t", row.names = FALSE)
  write.table(data.frame(ID=rownames(rankingH),rankingH), txt.fileH, sep = "\t", row.names = FALSE)
  write.table(data.frame(RefNumber=rownames(samplesSubclass),samplesSubclass), txt.fileS, sep = "\t", row.names = FALSE)

  message(paste("Both PDF report and TXT table with feature statistics included in 'featureTable' slot was generated succesfully."))

  dev.off()
}

##################################################################################################################################
##################################################################################################################################
# Intern functions to PDF report

.plotHeatmapDECO <- function(deco, p.val, id, analysis, symbol, infoS, color.cluster, data0,
                             cex.legend, cex.names, cex.samples, print.annot,
                             info.sample, info.feature, samplesSubclass)
{
  # Setting initial variables
  if(print.annot)
    message("NOTE: IDs have been mapped to HGNC symbol for 'Heatmap' plot.")

  i <- 1
  main <- list(NA)
  m <- 0.05
  n <- -0.03
  p.val <- rev(p.val)
  if(length(p.val) == 1 & p.val[1] > 0.05)
    message("\nWARNING: Heatmap from 'h' statistic corresponding to ALL samples\n shows NSCA p-value > 0.05.")
  else if(length(p.val) == 2 & any(p.val > 0.05))
    message(paste("\nWARNING: Heatmap from 'h' statistic corresponding to",c("CASES","CONTROL")[p.val > 0.05],"samples
                  shows NSCA p-value > 0.05."))
  while(i %in% 1:2)
  {
    # Taking main data depending on RDA design.
    if(analysis != "Binary"){
      i <- 3; MX <- as.matrix(deco@NSCAcluster$All$NSCA$h)
      nsca <- deco@NSCAcluster$All
      y1 <- range(deco@NSCAcluster$All$hclustSamp$cluster)[1]
      y2 <- range(deco@NSCAcluster$All$hclustSamp$cluster)[2]
      if(is.na(main)){main[i]="ALL samples"}
      info.dend <- deco@NSCAcluster$All$hclustSamp
      info.dendF <- deco@NSCAcluster$All$hclustFeat}
    if(i == 1){
      MX <- as.matrix(deco@NSCAcluster$Case$NSCA$h)
      nsca <- deco@NSCAcluster$Case
      y1 <- max(deco@NSCAcluster$Control$hclustSamp$cluster) + 1
      y2 <- max(deco@NSCAcluster$Case$hclustSamp$cluster) + y1
      main[i]="CASE (1) samples"}
    if(i == 2){
      MX <- as.matrix(deco@NSCAcluster$Control$NSCA$h)
      nsca <- deco@NSCAcluster$Control
      y1 <- 1
      y2 <- max(deco@NSCAcluster$Control$hclustSamp$cluster)
      main[i]="CONTROL (0) samples"}

    count <- table(deco@classes)
    # Labels to plot (original IDs or annotated ones).
    if(print.annot)
    {
      labels <- as.character(symbol[na.omit(match(rownames(MX),names(symbol)))])
      id.row <- c(rep("black",dim(MX)[1]))
      id.row[names(symbol[na.omit(match(rownames(MX),names(symbol)))]) %in% id] <- "red"
      names(id.row) <- labels
    }else
    {
      labels <- rownames(MX)
      id.row <- c(rep("black",dim(MX)[1]))
      id.row[rownames(MX) %in% id] <- "red"
      names(id.row) <- labels
    }
    # UP or DOWN label color in case of 'Supervised' design.
    y <- 0
    if(analysis == "Binary"){
      feature.label <- as.vector(deco@featureTable[rownames(MX),c("UpDw")])
      if(i==1){feature.label[feature.label == "DOWN"] <- "chartreuse4"
      feature.label[feature.label == "UP"] <- "firebrick2"
      feature.label[feature.label == "MIXED"] <- "black"
      y <- max(deco@NSCAcluster$Control$hclustSamp$cluster)
      info.dend <- deco@NSCAcluster$Case$hclustSamp
      info.dendF <- deco@NSCAcluster$Case$hclustFeat}
      if(i==2){feature.label[feature.label == "UP"] <- "chartreuse4"
      feature.label[feature.label == "DOWN"] <- "firebrick2"
      feature.label[feature.label == "MIXED"] <- "black";
      y <- 0
      info.dend <- deco@NSCAcluster$Control$hclustSamp
      info.dendF <- deco@NSCAcluster$Case$hclustFeat}
    }else{
      feature.label <- rep(NA,dim(MX)[1])
    }
    # If any sample information was input, we need generating label colors for each one.
    info.flag <- FALSE
    if(all(!is.na(info.sample)))
    {
      info.sample.color <- infoS$col[colnames(MX),]
      cc <- c(rep("white",dim(MX)[2]))
      for(j in 1:max(info.dend$cluster))
        cc[info.dend$cluster == j] <- color.cluster[j+y]
      # Final label matrix to place under samples dendrogram.
      info.sample.color <- cbind(info.sample.color, c(rep(NA,dim(MX)[2])), DECO.Subclass = cc)
      colnames(info.sample.color)[1:dim(info.sample)[2]] <- colnames(infoS$col)
      m <- 0.15
      n <- -0.05
    }
    else
    {
      info.sample.color <- c(rep(NA,dim(MX)[2]))
      for(j in 1:max(info.dend$cluster))
        info.sample.color[info.dend$cluster == j] <- color.cluster[j+y]
      info.sample.color <- cbind(DECO.Subclass = info.sample.color, rep(NA, dim(MX)[2]))
    }
    # Same procedure should be done for feature information.
    if(all(is.na(info.feature)))
    {
      infoFpat <- colorRampPalette(c("grey","black"))(max(info.dendF$cluster))[info.dendF$cluster]
      names(infoFpat) <- rownames(MX)
      if(analysis == "Binary")
        feature.color <- cbind(Feature.pattern = infoFpat, rep("white",length(feature.label)), Exprs = feature.label)
      else
        feature.color <- cbind(Feature.pattern = infoFpat, rep("white",length(feature.label)))
      rownames(feature.color) <- rownames(MX)
    }else{
      info.feature <- cbind(info.feature)
      info.feature.color <- info.feature
      infoFpat <- colorRampPalette(c("grey","black"))(max(info.dendF$cluster))[info.dendF$cluster]
      names(infoFpat) <- rownames(MX)
      for(z in 1:length(unlist(apply(info.feature,2,unique)))){info.feature.color[
        info.feature == unlist(apply(info.feature,2,unique))[z]] <- as.character(
          rep(colorRampPalette(RColorBrewer::brewer.pal(name = "Paired",n = 10))(length(unlist(apply(info.feature,2,unique))))[z],
              length(info.feature.color[info.feature == unlist(apply(info.feature,2,unique))[z]])))}
      # Final label matrix to place under samples dendrogram.
      feature.color <- cbind(Feature.pattern = infoFpat, Exprs = feature.label, rep("white",length(feature.label)),
                             Info.gene = info.feature.color)
    }
    # h-statistic heatmap including sample and feature information.
    .heatplot.2(dataset = as.matrix(MX), Colv = as.dendrogram(info.dend$dend),
                Rowv = as.dendrogram(info.dendF$dend), margins = c(14,12), cexRow = cex.names,
                cexCol = cex.names+m, lmat=rbind(c(6,0,5,0),c(0,0,5,0),c(0,0,2,0),c(4,1,3,0)),
                lhei=c(0.25,0.1,m,1),lwid=c(0.5,0.08,1.5,0.15),
                method = info.dend$dend$method,
                labCol = rownames(samplesSubclass)[na.omit(match(colnames(MX),samplesSubclass[,c("Samples")]))],
                labRow = labels, labRowCol = id.row, ColSideColorsSize = dim(info.sample.color)[2],
                RowSideColorsSize = dim(feature.color)[2],
                main = "", dualScale = FALSE, KeyValueName = "h-statistic",
                notecex = cex.names, ColSideColors = info.sample.color,
                RowSideColors = t(feature.color), col.breaks = 256, rm.out = TRUE)
    # Adding more information to plot.
    title(main = list(paste("h-statistic.",main[i]), font = 2), line = 4.5)
    mtext(text = "Higher h-statistic indicates more feature-relevance for these samples.", side = 3, outer = TRUE, line = -3, las = 1)
    mtext(text = paste("f =",dim(MX)[1],"features for heatmap"), side = 4, las = 3, font = 2, cex = 1.3)
    # mtext(text = paste("Cophenetic correlation of samples dendogram = ",round(sort(info.dend$coph,decreasing = TRUE)[1], digits = 2),
    #                    "\nAgglomeration method: ", info.dend$dend$method, sep = ""), side = 1, font = 1, line = 5, cex = 1)
    # Legends including feature and sample information.
    if(all(!is.na(info.sample)))
    {
      legend("topleft", legend = infoS$ty, col = names(infoS$ty), pch = 15,cex = cex.legend-0.2,
             title = "Sample info provided by user", bty = 'n', xpd = TRUE, inset = c(0,0.1))
    }
    legend(title = "DECO feature-sample subclasses", "topright",
           legend = paste("Subclass",1:max(info.dend$cluster)),col = color.cluster[(y+1):(j+y)], pch = 15,
           cex = cex.legend-0.1, bty = 'n', xpd = TRUE, inset = c(0,0))
    if(i == 3)
      legend("bottomleft",legend = c("DE features","ID provided features"), title = "Features",
             col = c("black","red"), bty = 'n', cex = cex.legend-0.1, pch = 15, xpd = TRUE, inset = c(0,-0.05))
    else
      legend("bottomleft",legend = c("UP","DOWN","MIXED","ID provided"), title = "Class feature and IDs",
             col = c("firebrick2","chartreuse4","black","red"), bty = 'n', cex = cex.legend-0.1, pch = 15, xpd = TRUE, inset = c(0,n))
    i = i+1
  }
  ## Raw Data heatmap without classes division for 'Supervised' design.
  MX <- data0[rownames(data0)%in%as.character(deco@featureTable[,"ID"]),
              colnames(data0)%in%samplesSubclass[,c("Samples")]]
  # Annotation of feature names if it is required.
  if(print.annot)
  {
    labels <- as.character(symbol[na.omit(match(rownames(MX),names(symbol)))])
    id.row <- c(rep("black",dim(MX)[1]))
    id.row[names(symbol[na.omit(match(rownames(MX),names(symbol)))]) %in% id] <- "red"
    names(id.row) <- labels
  }else
  {
    labels <- rownames(MX)
    id.row <- c(rep("black",dim(MX)[1]))
    id.row[rownames(MX) %in% id] <- "red"
    names(id.row) <- labels
  }
  # Generating label colors from sample information.
  if(all(!is.na(info.sample))){
    info.sample.color <- infoS$col[colnames(MX),]
    cc <- as.character(deco@classes)
    names(cc) <- names(deco@classes)

    if(analysis == "Binary"){
      cc[deco@classes == deco@control] <- "navyblue"
      cc[deco@classes != deco@control] <- "darkred"
      info.sample.color <- cbind(info.sample.color, c(rep(NA,dim(MX)[2])), DECO.Contrast.design = cc[colnames(MX)])
    }else if(analysis == "Unsupervised"){
      info.sample.color <- cbind(info.sample.color); n <- 0.07
    }else{
      col <- colorRampPalette(c("red","blue"))(length(table(deco@classes)))
      names(col) <- levels(deco@classes)
      cc <- col[cc]
      names(cc) <- names(deco@classes)
      info.sample.color <- cbind(info.sample.color, c(rep(NA,dim(MX)[2])), DECO.Contrast.design = cc[colnames(MX)])
      n <- 0.07
    }
    if(analysis == "Multiclass"){
      cols <- col
      col <- c(paste("Contrast design:",names(cols)),infoS$ty)
      names(col)[1:length(levels(deco@classes))] <- cols
    }else
      col <- infoS$ty
    colnames(info.sample.color)[1:dim(info.sample)[2]] <- colnames(infoS$col)
    n <- 0.1
  }
  else{
    cc <- as.character(deco@classes)
    names(cc) <- names(deco@classes)
    cc[deco@classes == deco@control] <- "navyblue"
    cc[deco@classes != deco@control] <- "darkred"
    if(analysis == "Binary")
      info.sample.color <- cbind(DECO.Contrast.design = cc[colnames(MX)], rep(NA, dim(MX)[2]))
    else
      info.sample.color <- cbind(rep(NA, dim(MX)[2]))
    n <- 0.05
  }
  # Hierarchical clustering of samples based on raw data and just for top-50 features.
  Colv <- cophDECO(data = MX, k = 2, verbose = FALSE, method.heatmap = "ward.D", coph = TRUE)
  Rowv <- cophDECO(data = t(MX), k = 2, verbose = FALSE, method.heatmap = "ward.D", coph = FALSE)
  # RAW data including NA's should be substituted.
  if(any(is.na(MX))){
    MX[is.na(MX)] <- min(na.omit(MX))
    message("WARNING: RAW data includes NA's, so they will be substituted by MINIMUM value to represent heatmap.")
  }
  par(las = 1)
  # HEATMAP representation of raw data.
  .heatplot.2(dataset = MX, Colv = as.dendrogram(Colv$dend), Rowv = as.dendrogram(Rowv$dend), margins = c(14,12), cexRow = cex.names, cexCol = cex.names+m,
             lmat=rbind(c(6,0,5,0),c(0,0,5,0),c(0,0,2,0),c(4,1,3,0)), lhei=c(0.25,n,m,1),lwid=c(0.5,0.02,1.5,0.15),
             method = Colv$dend$method, labCol = rownames(samplesSubclass), labRow = labels, labRowCol = id.row,
             ColSideColorsSize = dim(info.sample.color)[2]-1, RowSideColorsSize = 1, main = "", KeyValueName = "Raw data",
             notecex = cex.names, ColSideColors = info.sample.color,
             RowSideColors = c(rep(NA,dim(MX)[1])), dualScale = FALSE, cols.default = FALSE)
  # Adding some information, like legends for sample and feature information.
  title(main = list(paste("Raw data for differential features. All samples are represented."), font = 2), line = 4.5)
  mtext(text = "Classic raw data heatmap for features selected by DECO.", side = 3, outer = TRUE, line = -3, las = 1)
  mtext(text = paste("f =",dim(MX)[1]," features for heatmap"), side = 4, las = 3, font = 2, cex = 1.3)
  # mtext(text = paste("Cophenetic correlation of samples dendogram = ",round(sort(Colv$coph,decreasing = T)[1], digits = 2),
  #                    "\nAgglomeration method: ", Colv$dend$method, sep = ""), side = 1, font = 1, line = 5, cex = 1)
  if(all(!is.na(info.sample)))
  {
    legend("topleft", legend = col, col = names(col), pch = 15,cex = cex.legend-0.2,title = "Sample info provided by user",
           bty = 'n', xpd = TRUE, inset = c(0,0.1))
  }
  if(analysis == "Unsupervised")
    legend("bottomleft",legend = c("Differential features","ID provided features"), title = "Profile",col = c("black","red"),
           bty = 'n', cex = cex.legend-0.1, pch = 15, xpd = TRUE, inset = c(0,-0.05))
  else
    legend("bottomleft",legend = "ID provided", title = "Highlighted IDs",col = "red", bty = 'n', cex = cex.legend+0.1,
           pch = 15, xpd = TRUE, inset = c(0,n))

}


##################################################################################################################################


.plotProfileDECO <- function(deco, id, data, analysis, infoS, color.cluster, info.sample,
                             infoSubclass, samplesSubclass, cex.legend, cex.names, cex.samples,
                             print.annot, symbol)
{
  if(any(id %in% as.character(deco@featureTable[,c("ID")])))
  {
    flag = TRUE
    if(analysis != "Unsupervised"){
      classes <- deco@classes
      layout(mat = matrix(c(1,1,1,6,6,3,4,5,2,2),
                          ncol = 5, byrow = TRUE), heights = c(1.5,1))
    }else{
      if(all(!is.na(info.sample)))
        layout(mat = matrix(c(1:8), ncol = 4, byrow = TRUE), widths = c(5,3,1.5,1.5))
      else{
        layout(mat = matrix(c(1:6), ncol = 3, byrow = TRUE), widths = c(3,2,1))
        flag = FALSE
      }
    }
    id_sus <- rownames(deco@NSCAcluster[[1]]$NSCA$h)
    limy <- range(na.omit(data))
    id0 <- as.character(id[id %in% id_sus])
    main <- paste("Raw data profile of:",id0,sep=" ")
    if(any(!(id %in% id_sus))){
      message("NOTE: Following IDs have not been found by DECO or they have not been included in NSCA
              (if COMPLETE profiles were removed), so 'Profile' will be not plotted for them:")
      message(paste(id[!(id %in% id_sus)], "\n"))}
    data0 <- data

    # This plot will be repeated along all IDs input.
    for(a in 1:length(id0))
    {
      par(mar=c(10, 6, 10, 6))
      data <- data0
      # Taking 'h' statistic for each ID and then clustering samples using this information.
      if(analysis == "Binary")
      {
        tab.info <- c(deco@NSCAcluster$Control$NSCA$h[as.character(id0)[a],],
                      deco@NSCAcluster$Case$NSCA$h[as.character(id0)[a],])
        tab.info[is.na(tab.info)] <- 0
        # Clustering NSCA inner product for one feature.
        patt <- innerProductAssign(inner = tab.info, samples = deco@classes, control = deco@control, analysis)$cl
        count <- table(patt)
        thr <- which(cumsum(count) >= table(deco@classes)[names(table(deco@classes))==deco@control])[1]
        colors <- c(colorRampPalette(c("blue","lightblue"))(thr),colorRampPalette(c("red","darkorange"))(length(count)-thr))
        names(count) <- c(paste("Pattern",1:thr,names(table(deco@classes))[1]),
                          paste("Pattern",1:(length(count)-thr),names(table(deco@classes))[2]))
        ord <- names(patt)
        bg.col <- adjustcolor(c(rep("navyblue",table(deco@classes)[names(table(deco@classes)) == deco@control]),
                                rep("darkred",table(deco@classes)[names(table(deco@classes)) != deco@control])), 0.6)
        limb <- list(range(deco@NSCAcluster$Control$rankingFeature.h[,seq(2,dim(deco@NSCAcluster$Control$infoSubclass)[1]*2,2)]),
                     range(deco@NSCAcluster$Case$rankingFeature.h[,seq(2,dim(deco@NSCAcluster$Case$infoSubclass)[1]*2,2)]))
        limb <- limb[[which(unlist(lapply(limb, diff)) == max(unlist(lapply(limb, diff))))]]
      }else
      {
        tab.info <- c(deco@NSCAcluster$All$NSCA$h[as.character(id0[a]),])
        tab.info[is.na(tab.info)] <- 0
        # Clustering NSCA inner product for one feature.
        patt <- innerProductAssign(inner = tab.info, samples = deco@classes, analysis = analysis)$cl
        count <- table(patt)
        names(count) <- paste("Pattern",1:length(count))
        colors <- colorRampPalette(c("black","grey"))(length(count))
        ord <- names(patt)
        bg.col <- adjustcolor(rep("darkred",length(names(patt))),0.6)
        limb <- range(deco@NSCAcluster$All$rankingFeature.h[,seq(2,length(color.cluster)*2,2)])
      }
      # Ordering raw data using clustering information.
      data <- data[,ord]
      names(colors) <- names(count)
      count.col <- unlist(sapply(colors,function(x) rep(x,count[which(colors == x)])))

      # Plot feature profile.
      plot(x = 1:dim(data)[2], y = data[as.character(id0[a]),], ylim = limy, pch=21, xaxt="n",
           xlim=c(1,ncol(data)), main= paste(main[a],"\nsorted by h-statistic ranking"),
           type = "p", ylab = "Raw Data Signal", col="white",
           bg = bg.col, xlab="", cex = cex.samples+0.5)
      axis(side = 1,at = 1:dim(data)[2], labels = rownames(samplesSubclass[match(ord,samplesSubclass[,c("Samples")]),]),
           las=2, cex.axis = cex.names)
      m <- sum(abs(limy))
      if(print.annot){
        id_hgnc <- symbol[names(symbol)%in%as.character(id0)[a]]
        mtext(text = paste("HGNC symbol:",id_hgnc), line = 1)
        note <- TRUE}
      par(xpd=NA)
      # Including sample information under the expression plot.
      if(all(!is.na(info.sample))){
        cc <- rep("white",dim(samplesSubclass)[1])
        for(j in 1:length(color.cluster))
          cc[samplesSubclass[,"Subclass"] == names(color.cluster)[j]] <- color.cluster[j]
        names(cc) <- as.character(samplesSubclass[,"Samples"])
        cc <- cbind(infoS$col, DECO.Subclasses = cc[rownames(infoS$col)])
        for(i in 1:dim(cc)[2]){
          rect(0.5:(length(data[as.character(id0[a]),])-0.5),xright = 1.5:(length(data[as.character(id0[a]),])+0.5), border = NA,
               limy[1]-(m*0.05*i+0.1*m), ytop = limy[1]-(m*0.05*(i+1)+0.1*m), col = cc[colnames(data), i])
          text(x = length(data[as.character(id0[a]),])+1, labels = colnames(cc)[i],
               y = limy[1]-(m*0.05*(i+0.5)+0.1*m), adj = c(0,0.5))
        }
      }
      par(xpd=FALSE)
      # Adding dashed lines of mean values for each feature pattern found.
      for(i in 1:(length(count)))
      {
        if((max(which(patt == names(count)[i]))+0.5) < length(count.col) & analysis == "Binary"){
          thr <- dim(deco@NSCAcluster$Control$samplesSubclass)[1]+0.5
          abline(v = thr, col = "grey", lwd = 2, lty = 2)
          }
        if(i == 1){x0=0.75;x1=count[i]+0.25;y0 <- mean(data[id0[a],1:cumsum(count)[i]])}
        else{x0 = sum(count[1:(i-1)])+0.75; x1 = cumsum(count)[i]+0.25;y0 <- mean(data[id0[a],(cumsum(count)[(i-1)]+1):cumsum(count)[i]])}
        if((x1 - x0) > 1) segments(x0 = x0,x1 = x1,y0 = y0,y1 = y0,col = "black",lty = 1, lwd = 3)
      }
      # Sample information and feature statistics will be plotted on the right side.
      title(main = list("RDA & NSCA: Single feature patterns", cex = 1.5,
                        font = 2), outer = TRUE, line = -2)
      par(mar = c(7,5,4,5))
      # 'h' statistic per Subclass of samples.
      d <- data.frame(samplesSubclass, h = tab.info[samplesSubclass[,"Samples"]])
      m <- sapply(names(color.cluster), function(x) mean(d[d[,"Subclass"] == x,"h"]))
      s <- sapply(names(color.cluster), function(x)
        sd(d[d[,"Subclass"] == x,"h"])/sqrt(length(d[d[,"Subclass"] == x,"h"])))
      par(xpd=NA)
      plot(0,xlim = c(0,length(m)), ylim = limb, type = "n", axes = FALSE,
           xlab = "", ylab = "h mean per DECO subclass")
      h <- sapply(rownames(infoSubclass), function(x) d[d[,"Subclass"] == x,"h"])
      rect(xleft = 0:(length(m)-1), xright = 1:length(m)-0.1, ybottom = 0, ytop = m,
           col = adjustcolor(color.cluster, 0.6), border = NA)
      segments(x0 = seq(0.45,length(m),1), x1 = seq(0.45,length(m),1), y0 = m-s, y1 = m+s)
      segments(x0 = seq(0.45,length(m),1)-0.1, x1 = seq(0.45,length(m),1)+0.1, y0 = m-s, y1 = m-s)
      segments(x0 = seq(0.45,length(m),1)-0.1, x1 = seq(0.45,length(m),1)+0.1, y0 = m+s, y1 = m+s)
      lab <- unlist(lapply(sapply(names(color.cluster), function(x) strsplit(x, split = " Subclass ")),
                           function(x) paste(x,collapse = ".")))
      axis(1, labels = lab, at = seq(0.45,length(m),1),
           lty = 0, las = 2, cex.axis = 1.3)
      axis(2, las = 2, cex.axis = 1.5)
      par(xpd=FALSE)
      abline(h = 0)
      mtext(text = "h-statistic mean per subclass", side = 1, line = 10, font = 2)

      # Color labels and summary of statistics per feature.
      par(xpd = TRUE)
      if(analysis == "Binary"){
        plot(1, type = "n", axes=FALSE, xlab="", ylab="")
        legend("center",legend = paste(c("CONTROL","CASE"),"samples"),
               col = c("navyblue","darkred"), pch = 16, bty = "n", pt.cex = 2, title = "Constrast design:\ncolor of point",
               xpd = TRUE, cex = cex.legend+0.5)}
      else if(analysis == "Multiclass"){
        plot(1, type = "n", axes=FALSE, xlab="", ylab="")
        legend("center",legend = levels(classes),
               col = colorRampPalette(c("navyblue","darkorange","darkred"))(length(levels(classes))),
               pch = 16, bty = "n", pt.cex = 2, title = "Constrast design:\ncolor of point",
               xpd = TRUE, cex = cex.legend+0.5)
      }

      if(flag)
        plot(1, type = "n", axes=FALSE, xlab="", ylab="")
      if(all(!is.na(info.sample))){
        legend("center", legend = infoS$ty, col = names(infoS$ty), xpd = TRUE,
             pch = 15, cex = cex.legend+1/(length(infoS$ty))+0.5, title = "Sample info provided by user", bty = 'n')
      }
      if(analysis == "Binary"){
        textplot(rbind(round(t(deco@featureTable[as.character(id0)[a],colnames(deco@featureTable)%in%
                                                  c("Standard.Chi.Square","overlap.Ctrl.Case","Repeats","delta.signal",
                                                    "sd.Ctrl","sd.Case","h.Range.Ctrl","h.Range.Case","Dendrogram.group")]),3),
                       Profile = as.character(deco@featureTable[as.character(id0)[a],c("Profile")]),
                       UpDw = as.character(deco@featureTable[as.character(id0)[a],c("UpDw")])), halign = "right")}

      else{
        textplot(round(t(deco@featureTable[as.character(id0)[a],c("Standard.Chi.Square","Repeats","sd","h.Range","Dendrogram.group")]),3),
                 halign = "right")
      }
      if(analysis %in% c("Binary","Multiclass"))
        overlapFeature(as.character(id0[a]), data0, classes = deco@classes, analysis = analysis,
                       control = deco@control, infoS, plot = TRUE)
      par(xpd = FALSE)
    }
    remove(m)
  }else
    message("NOTE: Input IDs do not match with any DE feature. 'Expression' plot will be removed.")
  }

##################################################################################################################################

.plotDECObiplot <- function(deco, color.cluster, cex.samples)
{
  layout(matrix(1))
  color.cluster0 <- color.cluster
  par(mar = c(9,10,9,9))
  layout(matrix(c(1,3,2,2),nrow = 2, byrow = TRUE), heights = c(2,0.8))
  for(i in 1:length(deco@NSCAcluster)){
    color.cluster <- color.cluster0
    color.cluster <- color.cluster[names(color.cluster)%in%rownames(deco@NSCAcluster[[i]]$infoSubclass)]
    mx <- deco@NSCAcluster[[i]]$Biplot.coordinates[rownames(deco@NSCAcluster[[i]]$Biplot.coordinates)%in%colnames(deco@incidenceMatrix),]
    mxg <- deco@NSCAcluster[[i]]$Biplot.coordinates[!rownames(deco@NSCAcluster[[i]]$Biplot.coordinates)%in%colnames(deco@incidenceMatrix),]

    samp <- cbind(deco@NSCAcluster[[i]]$samplesSubclass, col = color.cluster[deco@NSCAcluster[[i]]$samplesSubclass])
    centroid <- matrix(0, nrow = length(color.cluster), ncol = 2, dimnames = list(names(color.cluster), 1:2))
    s <- c()
    for(j in 1:length(color.cluster)){
      centroid[j,] <- apply(mx[rownames(samp[samp[,1] %in% names(color.cluster)[j],]),1:2],2,mean)
      s[j] <- mean(apply(mx[rownames(samp[samp[,1] %in% names(color.cluster)[j],]),1:2],2,sd))
    }
    ###
    plot(mx, xlab = paste("Dim 1 NSCA, ",round(deco@NSCAcluster[[i]]$NSCA$Inertia[1,2],2),"% of total variability"),
         ylab = paste("Dim 2 NSCA, ",round(deco@NSCAcluster[[i]]$NSCA$Inertia[2,2],2),"% of total variability"),
         pch = 22, cex = cex.samples+1,
         bg = adjustcolor(samp[,"col"],0.5), col = adjustcolor("grey",0.3))
    abline(h = 0, v = 0, col = adjustcolor("black",0.3))
    points(centroid, cex = 0.7, pch = 21, lwd = 3)
    points(centroid, cex = s/max(s) * 16, pch = 21, col = adjustcolor(color.cluster,0.6), lwd = 4)
    text(centroid-sum(abs(range(mx[,1:2])))*0.01, labels = rownames(centroid), cex = 1.2, xpd = NA, font = 2)
    ###
    plot(0, type = "n", axes = FALSE, xlab = "", ylab = "")
    legend("center", legend = c(paste(rownames(centroid), ", n =",deco@NSCAcluster[[i]]$infoSubclass[,"Samples"],"samples"),
                                "Middle NSCA coordinates\nper subclass"), horiz = TRUE,
           pt.bg = c(color.cluster[rownames(centroid)],"white"), pch = c(rep(22,length(color.cluster)),21),
           xpd = NA, inset = c(-0.25,0), bty = "n", cex = 1.5, col = c(rep("grey",length(color.cluster)),"black"))
    if(dim(mx)[2] >= 3)
      scatterplot3d(mx[,1:3], bg = adjustcolor(samp[,"col"],0.6), color = adjustcolor("grey",0.3), mar = c(8,8,8,8),
                    pch = 22, cex.symbols = 2,
                    xlab = paste("Dim 1 NSCA, ",round(deco@NSCAcluster[[i]]$NSCA$Inertia[1,2],2),"% of total variability"),
                    ylab = paste("Dim 2 NSCA, ",round(deco@NSCAcluster[[i]]$NSCA$Inertia[2,2],2),"% of total variability"),
                    zlab = paste("Dim 3 NSCA, ",round(deco@NSCAcluster[[i]]$NSCA$Inertia[3,2],2),"% of total variability"))
    else{
      plot(0, type = "n", axes = TRUE, xlab = "", ylab = "")
      message(paste("NOTE: NSCA analysis on",names(deco@NSCAcluster)[[i]],
                    "shows less than 3 dimensions,\nthen 3D plot will be removed."))
    }
    title(list(paste("Biplot and three-dimensional plot of NSCA coordinates\n",names(deco@NSCAcluster)[i], "samples"),
               font = 2, cex = 1.4), outer = TRUE, line = -5)

  }
}

##################################################################################################################################
.plotRepThr <- function(sub, deco, id = NA, print.annot = FALSE)
{
  if(all(is.na(id)))
    id <- unique(c(as.character(deco@featureTable[1:5, "ID"]),as.character(deco@featureTable[order(
      deco@featureTable$Repeats.index, decreasing = TRUE), "ID"][1:5])))
  if(print.annot & "SYMBOL" %in% colnames(sub$subStatFeature))
    names <- as.character(deco@featureTable[id, "SYMBOL"])
  else
    names <- id

  g.names <- sapply(sort(rownames(sub$incidenceMatrix)), function(x) unlist(strsplit(x,split = "deco",fixed = TRUE))[1])

  if(all(is.na(sub$classes)) || length(levels(sub$classes)) > 2){
    x <- apply(sub$incidenceMatrix[sort(rownames(sub$incidenceMatrix)),], 1,
               function(x) length(which(x > 0)))/dim(sub$incidenceMatrix)[2]
    z <- apply(sub$incidenceMatrix[sort(rownames(sub$incidenceMatrix)),], 1,
          function(x) length(which(x >= deco@rep.thr)))/dim(sub$incidenceMatrix)[2]
    m <- max(sub$incidenceMatrix)
    names(z) <- g.names[names(g.names)%in%names(z)]
    names(x) <- names(z)
  }
  else{
    x <- sapply(unique(g.names), function(x)
      sum(apply(sub$incidenceMatrix[grepl(rownames(sub$incidenceMatrix), pattern = x, fixed = TRUE),],1,
                function(x) length(which(x > 0))))/dim(sub$incidenceMatrix)[2])
    z <- sapply(unique(g.names), function(x)
      sum(apply(sub$incidenceMatrix[grepl(rownames(sub$incidenceMatrix), pattern = x, fixed = TRUE),],1,
                function(x) length(which(x > deco@rep.thr))))/dim(sub$incidenceMatrix)[2])
  }

  y <- sub$subStatFeature[order(as.character(sub$subStatFeature$ID)),"Repeats"]
  names(y) <- sort(as.character(sub$subStatFeature$ID))
  y=y[names(x)]
  data <- data.frame(x,z)

  # Setting up Colors
  col <- sapply(y, function(j) max(which(seq(1,max(y),length.out = 16) <= j)))
  color <- colorRampPalette(c("orange","brown","green","blue"))(16)
  col <- color[col]
  col[!names(y) %in% deco@featureTable$ID] <- "red"
  n <- length(which(col == "red"))
  col <- adjustcolor(col, 0.3)

  z1=matrix(1:16,nrow=1)
  x1=1
  y1=seq(1,max(y),length.out = 16)

  # Plot
  layout(mat = matrix(c(1,2,3,3),ncol = 2, byrow = TRUE),
         widths = c(1,0.17), heights = c(1.3,0.5))
  par(mar = c(4,10,8,2))
  plot(data, type = "p", pch = 21, col = adjustcolor("white",0.1),
       bg = col, axes = FALSE, cex = 2,
       xlab = "% samples amounting at least 1 repeat", xlim = c(0,1),
       ylab = paste("% samples amounting threshold:",deco@rep.thr,"repeats"))
  axis(1, at = seq(0,1,0.1), font = 2)
  axis(2, font = 2, las = 2)
  for(a in 1:length(id))
    text(data[id[a],]-0.01, labels = names[a], font = 2, col = adjustcolor("black",0.8))
  legend("topleft", legend = paste(n,"features within <= ",round(deco@samp.perc,3)*100,"% samples with at least",
                                   deco@rep.thr,"repeats."),
         bty = "n", cex = 1.5, col = "red", pch = 16, lty = NA, lwd = 2)
  par(mar = c(3,3,3,8))
  image(x = 0:1,y1,z1,col=color, cex = 1.2,
        axes=FALSE,xlab="",ylab=paste("Repeats"))
  axis(2, at = seq(1,max(y),length.out = 16), labels = round(seq(1,max(y),length.out = 16),0), font = 2, las = 2)
  title(main = list("RDA: Repeats threshold based on differential events distribution per feature", font = 2, cex = 1.4),
        outer = TRUE, line = -5)

  layout(mat = 1)

}

################################################
################################################
# Effect Size plot

.plotEffectSize <- function(deco, infoS = NA, samplesSubclass, analysis)
{
  incid <- deco@incidenceMatrix
  if(analysis == "Binary"){
    s <- apply(incid[,names(deco@classes)], 2, sum)
    s <- c(sort(s[names(s)%in%names(deco@classes[deco@classes == deco@control])], decreasing = TRUE),
           sort(s[names(s)%in%names(deco@classes[deco@classes != deco@control])]))
    col.design <- c(rep("navyblue",table(deco@classes)[names(table(deco@classes)) == deco@control]),
                    rep("darkred",table(deco@classes)[names(table(deco@classes)) != deco@control]))
  }else{
    s <- sort(apply(incid, 2, sum))
    col.design <- "black"
  }
  s <- s[s > 0]
  limy <- c(0,range(s)[2])
  m <- sum(range(s))

  if(all(!is.na(infoS))){
    layout(matrix(c(1,2,3,3),ncol = 2,byrow = FALSE), widths = c(3,1), heights = c(2,0.7))
    par(mar = c(3,5,8,5), xpd = NA)
  }else
    par(mar = c(7,7,7,7), xpd = NA)
  plot(NA, ylim = limy, xlim = c(1,length(s)), xlab = "Reference number from samples",
       ylab = "Repeats per sample", axes = FALSE, type = "n")
  axis(1, at = 1:length(s), labels = rownames(samplesSubclass)[match(samplesSubclass[,1],names(s))],
       tick = FALSE, las = 2)
  axis(2, at = seq(0, limy[2], signif(m*0.1,0)), las = 2)
  abline(h = seq(0, limy[2], signif(m*0.1,0)), lty = 2, xpd = FALSE, col = "grey")
  rect(1:length(s)-0.35,xright = 1:length(s)+0.35, ybottom = 0, ytop = s,
       xlab = "", ylab = "",
       col = adjustcolor(col.design, alpha.f = 0.7), border = "grey")
  par(mar = c(5,5,5,5))
  if(analysis == "Binary"){
    abline(v = length(deco@classes[deco@classes == deco@control])+0.5, xpd = FALSE, lty = 2, lwd = 3)
    text(x = c(length(deco@classes[deco@classes == deco@control])/2, length(deco@classes[deco@classes != deco@control])/2+
                 length(deco@classes[deco@classes == deco@control])), font = 2, cex = 1.2,
         y = limy[2] + 0.05*sum(limy), labels = c("CONTROL","CASE"), col = c("navyblue","darkred"))}
  if(all(!is.na(infoS))){
    infoS$orig <- as.matrix(infoS$orig[names(s),])
    infoS$col <- as.matrix(infoS$col[names(s),])
    plot(0, ylim = c(0,dim(infoS$orig)[2]), xlim = c(1,dim(infoS$orig)[1]), xlab = "", ylab = "", axes = FALSE, type = "n")
    for(i in 1:dim(infoS$col)[2]){
      rect(0.5:(length(s)-0.5),xright = 1.5:(length(s)+0.5), border = NA,
           ybottom = i-1, ytop = i, col = infoS$col[names(s), i])
      text(x = length(s) + 2, labels = colnames(infoS$col)[i],
           y = i-0.5, adj = c(0,0.5))}
    plot(0, xlab = "", ylab = "", axes = FALSE, type = "n")
    legend("center", legend = infoS$ty, col = names(infoS$ty), pch = 15,
           cex = 1.2, title = "Sample info provided by user", bty = 'n')
  }
  par(xpd = FALSE)

  title(list("RDA: Differential events counted per sample", cex = 1.3), outer = TRUE, line = -3)
  mtext(text = "Samples with higher amounts of 'Repeats' resemble\nmore different profiles from other class of samples.",
        outer = TRUE, line = -7, cex = 1.1)
}


################################################
################################################
# Heterogeinity plot

.plotHeterogeneityDECO <- function(v, var = 80, main = NA)
{
  if(!is.vector(v)){stop("Function need a vector of variabilities explained by NSCA.")}
  if(is.na(main)){main <- "Heterogeinity dataset analysis"}

  thr <- which(v >= var)[1]
  co <- rep("darkgoldenrod",length(v))
  co[thr] <- "red"

  x <- 1:length(v)
  lo1 <- loess(v~x)
  xl <- seq(min(x),max(x), (max(x) - min(x))/1000)

  plot(xl, predict(lo1,xl),type = "l", ylim = c(0,100), xlab = "Number of NSCA dimensions",
       ylab = "Variability explained", lwd = 5, col = "darkgoldenrod",
       axes = FALSE)
  axis(side = 1, at = 1:length(v), labels = 1:length(v))
  axis(side = 2, at = seq(from = 0, to = 100, by = 10))
  points(v, x = x, lwd = 5, col = co)
  abline(h = var, lty = 2, lwd = 2, col = "red")

  title(main = main, cex = 1.3)
  mtext(text = paste("Minimum dimensions to explain ", signif(var,4), "% : ", thr, sep=""), las = 1)

}

##############################
##############################
## plotTau

.plotTau <- function(x, samplesSubclass, labels = NA, n = 10, color.cluster,
                    id = NA, main.classes = NA)
{
  if(all(is.na(labels)))
  {
    labels <- sapply(names(sort(x[[1]]$NSCA$di, decreasing = TRUE)),function(y) unlist(strsplit(y,"deco"))[1])
    names(labels) <- labels
  }
  for(j in 1:2)
  {
    if(length(x)==2){
      if(is.na(main.classes))
        main <- c("CONTROL","CASE")
      else
        main <- main.classes
      gene <- cbind(Control = x[[1]]$NSCA$di, Case = x[[2]]$NSCA$di)
      if(n > dim(gene)[1])
        n <- dim(gene)[1]
      if(j == 1)
        gene <- gene[names(sort(apply(apply(gene,2,rank),1,mean),decreasing = TRUE)[1:n]),]
      else
        gene <- gene[names(sort(apply(apply(gene,2,rank),1,mean),decreasing = FALSE)[1:n]),]
      rownames(gene) <- sapply(rownames(gene),function(x)unlist(strsplit(x,split="deco"))[1])
      layout(mat = rbind(c(1,2),c(3,4)), widths = c(2,1), heights = c(1,1))
    }
    else{
      main <- "ALL"
      gene <- cbind(All = x[[1]]$NSCA$di)
      if(n > dim(gene)[1])
        n <- dim(gene)[1]
      if(j == 1)
        gene <- cbind(gene[names(sort(apply(apply(gene,2,rank),1,mean),decreasing = TRUE)[1:n]),])
      else
        gene <- cbind(gene[names(sort(apply(apply(gene,2,rank),1,mean))[1:n]),])
      layout(mat = cbind(1,2), widths = c(2,1.2))
    }
  }
  par(mar = c(6,8,10,6))
  for(i in 1:length(x))
  {
    d <- data.frame(x[[i]]$samplesSubclass, tau = sort(x[[i]]$NSCA$tau.num.j))
    names(x[[i]]$NSCA$tau.num.j) <- as.character(1:length(x[[i]]$NSCA$tau.num.j))
    if(i == 1)
      names(x[[i]]$NSCA$tau.num.j) <- as.character(1:length(x[[i]]$NSCA$tau.num.j))
    else
      names(x[[i]]$NSCA$tau.num.j) <- as.character((length(x[[i-1]]$NSCA$tau.num.j)+1):
                                                     (length(x[[i]]$NSCA$tau.num.j)+length(x[[i-1]]$NSCA$tau.num.j)))
    ### Formatting tau contribution into a list
    dat <- by(d$tau, d$Subclass, c)
    dat[[length(dat)+1]] <- na.omit(d[as.character(samplesSubclass[,"Samples"]),"tau"])
    names(dat[[length(dat)]]) <- rownames(samplesSubclass)[as.character(samplesSubclass[,"Samples"]) %in% rownames(d)]
    names(dat)[length(dat)] <- paste(main[i],"samples")
    ## Drawing boxplot...
    .customBoxplot(dat, col = c(color.cluster[rownames(x[[i]]$infoSubclass)],"black"),
                  ylab = "% Tau contribution", ylim = c(0,max(d[,"tau"])),
                  pt.cex = 2, main = paste(main[i],"samples"), print.name = TRUE)

    text <- round(rbind(Tau = x[[i]]$NSCA$tau, C.Statistic = x[[i]]$NSCA$Cstat, p.value = x[[i]]$NSCA$P.Value),digits = 5)
    colnames(text) <- paste(main[i],"NSCA")
    text.color <- matrix("black",nrow = 3, ncol = 1)
    if(x[[i]]$NSCA$P.Value > 0.05)
      text.color[3,1] <- "orange"
    if(x[[i]]$NSCA$P.Value > 0.1)
      text.color[3,1] <- "red"
    if(x[[i]]$NSCA$P.Value <= 0.05)
      text.color[3,1] <- "green"
    textplot(halign = "center",valign = "center", object = text, cex = 2, col.data = text.color)
  }
  title(main = list("NSCA: Goodman and Kruskal's Tau contribution", font = 2, cex = 1.5), line = -3, outer = TRUE)
  mtext(text = "A higher Tau value per sample indicates more sample-specific signal.", outer = TRUE, line = -5)
}

##############################
##############################
## plot ranking 'h' statistic

.plotRankingH <- function(deco, analysis, id = NA)
{
  layout(mat = matrix(1))
  if(analysis != "Binary")
    j <- 2
  else
    j <- 0
  for(i in 1:length(deco@NSCAcluster)){

    d <- deco@NSCAcluster[[i]]$rankingFeature.h[order(deco@NSCAcluster[[i]]$rankingFeature.h[,"h.Range"],decreasing = TRUE),]
    d <- d[1:min(c(50,dim(d)[1])),]
    colnames(d) <- paste(colnames(d),c("Ctrl","Case","All")[i+j])
    col <- c(rep("black",dim(deco@featureTable[rownames(d),colnames(deco@featureTable)%in%c("ID","SYMBOL","Standard.Chi.Square")])[2]),
             rep("deepskyblue4",length(which(grepl(colnames(d), pattern = "Scl", fixed = TRUE)))),rep("darkred",3))
    d <- data.frame(deco@featureTable[rownames(d),colnames(deco@featureTable)%in%c("ID","SYMBOL","Standard.Chi.Square")], d)
    color.table <- matrix("black", nrow = min(c(50,dim(d)[1])), ncol = dim(d)[2])
    color.table[,which(grepl(colnames(d), pattern = "Ranking", fixed = TRUE))][
      apply(d[, which(grepl(colnames(d), pattern = "Ranking", fixed = TRUE))], 2, function(x) x %in% 1:10)] <- "red"
    color.table[as.character(d[,1]) %in% id,1] <- "darkred"
    textplot(d, show.rownames = FALSE, col.colnames = col, mar = c(2,2,5,2), col.data = color.table)
    title(list(paste("Mean 'h' statistic per subclass within",c("CONTROL","CASE","ALL")[i+j],"samples"), cex = 1.4, font = 2),
          outer = TRUE, line = -4)
    mtext(paste("Top",min(c(50,dim(d)[1])),"discriminant features among subclasses found by DECO algorithm."),
          outer = TRUE, side = 3, line = -5.5)
  }
}


##############################
##############################
## stackpoly.2

stackpoly.2 <- function (x, y = NULL, main = "", xlab = "", ylab = "", xat = NA, axis1 = TRUE,
                         xaxlab = NA, xlim = NA, ylim = NA, lty = 1, lwd = 1, border = NA,
                         col = NULL, stack = FALSE, axis2 = TRUE, axis4 = TRUE,
                         padj = 0, ...)
{
  ydim <- dim(y)
  if (is.null(y[1])) {
    y <- x
    ydim <- dim(y)
    if (is.null(ydim))
      x <- 1:length(y)
    else x <- matrix(rep(1:ydim[1], ydim[2]), ncol = ydim[2])
  }
  if (stack)
    y <- t(unlist(apply(as.matrix(y), 1, cumsum)))
  if (is.na(xlim[1]))
    xlim <- range(x)
  if (is.na(ylim[1]))
    ylim <- range(y)
  plot(0, main = main, xlab = xlab, ylab = ylab, xlim = xlim,
       ylim = ylim, type = "n", xaxs = "i", yaxs = "i", axes = FALSE,
       ...)
  if (is.matrix(y) || is.list(y)) {
    plotlim <- par("usr")
    if (is.na(xat[1]))
      xat <- x[, 1]
    if (is.na(xaxlab[1]))
      xaxlab <- xat
    if(axis1)
      axis(1, at = xat, labels = xaxlab, padj = padj, las = 2)
    if (axis2)
      axis(2)
    if (axis4)
      axis(4)
    if (is.null(col[1]))
      col = rainbow(ydim[2])
    else if (length(col) < ydim[2])
      col <- rep(col, length.out = ydim[2])
    if (length(border) < ydim[2])
      border <- rep(border, length.out = ydim[2])
    if (length(lty) < ydim[2])
      lty <- rep(lty, length.out = ydim[2])
    if (length(lwd) < ydim[2])
      lwd <- rep(lwd, length.out = ydim[2])
    for (pline in seq(ydim[2], 1, by = -1)) {
      polygon(c(x[1], x[, pline], x[ydim[1]]), c(plotlim[3],
                                                 y[, pline], plotlim[3]), border = border[pline],
              col = col[pline], lty = lty[pline], lwd = lwd[pline])
    }
  }
  else {
    polygon(c(min(x), x, max(x), 0), c(0, y, 0, 0), border = border,
            col = col, lty = lty, lwd = lwd)
    if (is.na(xat[1]))
      xat <- x
    if (is.na(xaxlab[1]))
      xaxlab <- xat
    axis(1, at = xat, labels = xaxlab)
    if (axis2)
      axis(2)
    if (axis4)
      axis(4)
  }
}


##############################
##############################
## plotOD

.plotOD <- function(deco, id, ord, symbol, print.annot = FALSE)
{
  par(mfrow = c(1,1), fig = c(0,1,0,1), las = 1)
  layout(matrix(c(2,1,4,3,5,5), ncol = 3), widths = c(5,1.5,1.5), heights = c(1.5,5))
  deco@featureTable <- deco@featureTable[order(deco@featureTable$Standard.Chi.Square, decreasing = FALSE),]
  id.color <- rep(adjustcolor("black",alpha.f = 0.3),dim(deco@featureTable)[1])
  id.color[deco@featureTable$ID %in% id] <- adjustcolor("red",alpha.f=0.5)
  bo.color <- rep("white",dim(deco@featureTable)[1])
  s <- data.frame(as.numeric(deco@featureTable$delta.signal), 1 - deco@featureTable$overlap.Ctrl.Case,
                  row.names = as.character(deco@featureTable$ID))
  sd.thr <- c(quantSD(s[,1]),quantSD(s[,2]))
  l <- "Thresholds among different\nfeature profiles"
  par(mar = c(10,10,1,1))
  plot(s, pch = 21, xlab = "delta.signal", lwd = 2,
       ylab = "Non-overlap = 1 - Overlap",
       cex = 15 * deco@featureTable$Standard.Chi.Square/max(deco@featureTable$Standard.Chi.Square),
       col = bo.color, bg = id.color, cex.lab = 1.5)
  abline(v = 0, h = c(0.2,0.4,0.75),
         lty = 2, xpd = FALSE, lwd = 1)
  title(main = list("RDA: overlap Signal VS delta Signal plot", cex = 1.5), outer = TRUE, line = -4)
  high <- as.character(deco@featureTable[order(deco@featureTable$Standard.Chi.Square, decreasing = TRUE),]$ID)[1:10]
  text(s[high, ], labels = symbol[high], cex = 0.9, col = "darkred", font = 2)
  if(any(!id %in% high))
  {if(print.annot)
    text(s[id[!id %in% high],], labels = symbol[id[!id %in% high]], cex = 0.9, col = "black", font = 2)
    else
      text(s[id[!id %in% high],], labels = id[!id %in% high], cex = 0.9, col = "black", font = 2)}
  par(mar = c(2,10,10,1))
  x1 <- min(which(density(s[,1])$x <= 0))
  x2 <- max(which(density(s[,1])$x <= 0))
  x3 <- max(which(density(s[,1])$x > 0))
  plot(x = density(s[,1])$x, y = density(s[,1])$y, xlim = range(s[,1]),
       type = "l", main = "", lwd = 3, ylab = "Density", col = "black")
  with(density(s[,1]), polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col= adjustcolor("green",alpha.f=0.3)))
  with(density(s[,1]), polygon(x=c(x[c(x2,x2:x3,x3)]), y= c(0, y[x2:x3], 0), col= adjustcolor("red3",alpha.f=0.3)))
  par(mar = c(10,2,1,10))
  x1 <- min(which(density(s[,2])$x <= 0.2))
  x2 <- max(which(density(s[,2])$x <= 0.2))
  x3 <- max(which(density(s[,2])$x <= 0.4))
  x4 <- max(which(density(s[,2])$x <= 0.75))
  x5 <- max(which(density(s[,2])$x <= 1))
  plot(y = density(s[,2])$x, x = density(s[,2])$y, ylim = range(s[,2]),
       type = "l", main = "", xlab = "Density", lwd = 4, col = "black")
  if(is.finite(x1) & is.finite(x2))
    with(density(s[,2]), polygon(y=c(x[c(x1,x1:x2,x2)]), x= c(0, y[x1:x2], 0), col= adjustcolor("gray20",alpha.f=0.3)))
  if(is.finite(x2) & is.finite(x3))
    with(density(s[,2]), polygon(y=c(x[c(x2,x2:x3,x3)]), x= c(0, y[x2:x3], 0), col= adjustcolor("goldenrod1",alpha.f=0.3)))
  if(is.finite(x3) & is.finite(x4))
    with(density(s[,2]), polygon(y=c(x[c(x3,x3:x4,x4)]), x= c(0, y[x3:x4], 0), col= adjustcolor("chocolate2",alpha.f=0.3)))
  if(is.finite(x4) & is.finite(x5))
    with(density(s[,2]), polygon(y=c(x[c(x4,x4:x5,x5)]), x= c(0, y[x4:x5], 0), col= adjustcolor("brown4",alpha.f=0.3)))

  plot(1, type="n", axes=FALSE, xlab="", ylab="")
  plot(1, type="n", axes=FALSE, xlab="", ylab="")
  par(xpd = TRUE, mar = c(1,1,1,1))
  legend("left", pch = c(16,65,NA,15,15,15,15), lty = c(NA,NA,2,NA,NA,NA,NA), pt.cex = 2.5, lwd = c(NA,NA,2,NA,NA,NA,NA),
         legend = c("ID provided by user", paste("Top-15 features based\non Standard.Chi.Square"), l,
                    "COMPLETE features","MAJORITY features","MINORITY features","MIXED features"),
         col = c(adjustcolor("red",alpha.f=0.5),"darkred","black",adjustcolor(c("brown4","chocolate2","goldenrod1","gray20"),alpha.f=0.3)),
         bty = "n", cex = 2, xpd = TRUE, y.intersp = 2)
  par(xpd = FALSE)
  mtext("Circle size corresponds to relative amount of 'Standard.Chi.Square' per feature.
        Higher circles indicate more DIFFERENTIAL SIGNAL between both classes.", outer = TRUE, line = -8)
}

##############################
##############################
## plotRepThr

plotRepThr <- function(sub, deco, id = NA, print.annot = FALSE)
{
  if(all(is.na(id)))
    id <- unique(c(as.character(deco@featureTable[1:5, "ID"]),as.character(
      deco@featureTable[order(deco@featureTable$Repeats.index, decreasing = TRUE), "ID"][1:5])))
  if(print.annot & "SYMBOL" %in% colnames(sub$subStatFeature))
    names <- as.character(deco@featureTable[id, "SYMBOL"])
  else
    names <- id

  g.names <- sapply(sort(rownames(sub$incidenceMatrix)), function(x)
    unlist(strsplit(x,split = "deco",fixed = TRUE))[1])

  if(all(is.na(sub$classes))){
    x <- apply(sub$incidenceMatrix[sort(rownames(sub$incidenceMatrix)),], 1,
               function(x) length(which(x > 0)))/dim(sub$incidenceMatrix)[2]
    z <- apply(sub$incidenceMatrix[sort(rownames(sub$incidenceMatrix)),], 1,
               function(x) length(which(x > deco@rep.thr)))/dim(sub$incidenceMatrix)[2]
    names(z) <- g.names[names(g.names)%in%names(z)]
  }
  else{
    x <- sapply(unique(g.names), function(x)
      sum(apply(sub$incidenceMatrix[grepl(rownames(sub$incidenceMatrix), pattern = x, fixed = TRUE),],1,
                function(x) length(which(x > 0))))/dim(sub$incidenceMatrix)[2])
    z <- sapply(unique(g.names), function(x)
      sum(apply(sub$incidenceMatrix[grepl(rownames(sub$incidenceMatrix), pattern = x, fixed = TRUE),],1,
                function(x) length(which(x > deco@rep.thr))))/dim(sub$incidenceMatrix)[2])
  }

  y <- sub$subStatFeature[order(sub$subStatFeature$ID),"Repeats"]
  data <- data.frame(x,y)

  # Setting up Colors
  col <- sapply(z, function(x) max(which(seq(0,1,deco@samp.perc) <= x)))
  color <- adjustcolor(c("red",colorRampPalette(c("orange","brown","green","blue"))(length(seq(0,1,deco@samp.perc))-1)),0.6)
  col <- color[col]
  n <- length(which(col == "red"))

  z1=matrix(1:length(seq(0,1,deco@samp.perc)),nrow=1)
  x1=1
  y1=seq(0,1,deco@samp.perc)

  # Building a non-linear model
  m <- lm(y ~ poly(x,3), data)

  # Plot
  layout(mat = matrix(c(1,2,3,3),ncol = 2, byrow = TRUE), widths = c(1,0.17), heights = c(1.3,0.5))
  par(mar = c(4,10,8,2))
  plot(x,y, type = "p", pch = 21, col = adjustcolor("white",0.3), bg = col, axes = FALSE, cex = 2,
       xlab = "% samples affected per feature", ylab = "Repeats per feature")
  axis(1, at = seq(0,1,0.1), font = 2)
  axis(2, font = 2)
  lines(smooth.spline(x,predict(m)), lwd = 3, lty = 5)
  for(a in 1:length(id))
    text(data[grep(x = rownames(data), pattern = id[a], fixed = TRUE),]-0.01, labels = names[a], font = 2)
  legend("topleft", legend = c("Non-linear model describing correlation",
                               paste(n,"features <= ",round(deco@samp.perc,2)*100,"% samples")),
         bty = "n", cex = 1.5, col = c("black","red"), pch = c(NA,16), lty = c(5,NA), lwd = 2)
  par(mar = c(3,3,3,8))
  image(x = 0:1,y1,z1,col=color, cex = 1.2,
        axes=FALSE,xlab="",ylab=paste("% samples affected:",3,"repeats at least"))
  axis(2,at = seq(0,max(y1),0.05)-deco@samp.perc/2, labels = seq(0,max(y1),0.05), font = 2, las = 2)
  par(mar = c(8,10,5,8))
  textplot(round(summary(m)$coefficients, 4), cex = 1.3)
  mtext(paste("Calculated non-linear model, adjusted r-squared =",round(summary(m)$adj.r.squared,3)), las = 1)
  title(main = list("RDA: Repeats threshold based on differential events distribution per feature", font = 2, cex = 1.4),
        outer = TRUE, line = -5)

  layout(mat = 1)

}

##############################
##############################
## Plot dendrogram from samples

.plotDend <- function(dend, analysis, color.cluster, samplesSubclass, cex.names)
{
  par(mar = c(10,6,6,6), xpd = TRUE)
  if(analysis == "Binary")
  {
    h <- c(dend[[1]]$huber, dend[[2]]$huber)
    dend <- list(dend[[1]]$dend, dend[[2]]$dend)
    layout(matrix(c(1,2,3,3), nrow = 2, byrow = FALSE), widths = c(6, 1.5))
    main <- c("CONTROL samples", "CASE samples")
    col <- color.cluster[samplesSubclass[,2]]
    names(col) <- rownames(samplesSubclass)
    for(i in 1:2){
      dend[[i]]$labels <- rownames(samplesSubclass)[match(dend[[i]]$labels,samplesSubclass[,1])]
      plot(dend[[i]], axes = FALSE, cex = cex.names, lwd = 2, hang = -1, main = "",
           xlab = "", ylab = "", sub = "")
      rect(0.5:(length(dend[[i]]$labels)-0.5),xright = 1.5:(length(dend[[i]]$labels)+0.5), border = NA,
           ybottom = -0.25*max(dend[[i]]$height), ytop = -0.15*max(dend[[i]]$height),
           col = col[dend[[i]]$labels[dend[[i]]$order]])
      mtext(main[i], font = 2, cex = 1.2)
      mtext(side = 1, text = paste("Hubber's gamma coefficient for cutting dendrogram:", round(h[i],3)),
            font = 2, cex = 1.2, line = 5)
    }
    plot(0, axes = FALSE, type = "n", xlab = "", ylab = "")
    pos <- "center"
  }else
  {
    h <- dend$huber
    dend <- dend$dend
    col <- color.cluster[samplesSubclass[,2]]
    names(col) <- rownames(samplesSubclass)
    dend$labels <- rownames(samplesSubclass)[match(dend$labels,samplesSubclass[,1])]
    plot(as.dendrogram(dend), axes = FALSE, cex = cex.names)
    rect(0.5:(length(dend$labels)-0.5),xright = 1.5:(length(dend$labels)+0.5), border = NA,
         ybottom = -0.25*max(dend$height), ytop = -0.15*max(dend$height),
         col = col[dend$labels[dend$order]])
    pos <- "topright"
  }
  legend(pos, legend = names(color.cluster), col = color.cluster, y.intersp = 1.3,
         pch = 15, cex = 1.3, xpd = TRUE,
         bty = "n", title = "Subclasses of samples found")
  title(list("NSCA: Subclasses of samples found based on 'h' statistic", font = 2, cex = 1.5),
        outer = TRUE, line = -3)
  par(xpd = FALSE)
}


##############################
##############################
## Modified function of heatplot (made4 require). Remove all printings of information.

.heatplot.2 <- function (dataset, Colv, Rowv, dend = c("both", "row", "column", "none"),
                        cols.default = TRUE, lowcol = "green", highcol = "red", scale = "none",
                        classvec = NULL, classvecCol = NULL, classvec2 = NULL, distfun = NULL, labRowCol = NULL,
                        returnSampleTree = FALSE, method = "ave", dualScale = TRUE, col.breaks = 256,
                        scaleKey = TRUE, rm.out = FALSE, ...)
{
  data <- made4::array2ade4(dataset)
  data <- as.matrix(data)
  if (dualScale)
    data <- scale(data)
  # Separate outliers
  if(rm.out){
    data[data <= quantile(data, prob = 0.025)] <- quantile(data, prob = 0.025)
    data[data >= quantile(data, prob = 0.975)] <- quantile(data, prob = 0.975)
  }
  cols <- function() {
    if(lowcol == "green" & highcol == "red")
      hmcol <- colorRampPalette(c(rev(RColorBrewer::brewer.pal(10, "Reds"))[1:6],"white",
                                  RColorBrewer::brewer.pal(10, "Greens")[4:9]))(col.breaks)
    else
      hmcol <- colorRampPalette(c(highcol,lowcol))(col.breaks)
    return(rev(hmcol))
  }
  cols.gentleman <- function() {
    hmcol <- colorRampPalette(c(rev(RColorBrewer::brewer.pal(10, "Reds"))[1:6],"white",
                                RColorBrewer::brewer.pal(10, "Blues")[4:9]))(col.breaks)
    return(rev(hmcol))
  }
  if (cols.default)
    plotcols = cols.gentleman()
  else plotcols = cols()
  if (is.null(distfun))
    distf = distFunc
  else distf = function(x) dist(t(x), method = distfun)

  # Rowv <- FALSE
  dend <- match.arg(dend)
  dend = tolower(dend)
  if (dend %in% c("row", "r")) {
    # Rowv = as.dendrogram(hclust(distf(t(data)), method = method))
    # Rowv = as.dendrogram(cophDECO(t(data), k = 2, method.heatmap = "ward.D", verbose = F)$dend)
  }
  if (dend %in% c("column", "col", "c")) {
    dend = "column"
  }
  if (dend %in% c("both", "TRUE")) {
    # Rowv = as.dendrogram(hclust(distf(t(data)), method = method))
    # Rowv = as.dendrogram(cophDECO(t(data), k = 2, method.heatmap = "ward.D", verbose = F)$dend)
    dend = "both"
  }
  RSideColors = CSideColors = NULL
  if (any(!is.null(classvec), !is.null(classvec2))) {
    proc.classvec <- function(classvec) {
      classvec = as.factor(classvec)
      if (is.null(classvecCol))
        classvecCol = getcol(length(levels(classvec)))
      SideCols = factor(classvec, labels = classvecCol)
      SideCols = as.character(SideCols)
      nSC = length(SideCols)
      return(list(nSC, SideCols))
    }
    if (!is.null(classvec)) {
      out = proc.classvec(classvec)
      nSC = out[[1]]
      SideCols = out[[2]]
      if (!nSC %in% dim(data))
        print("Error: classvec length not equal to nrow or ncol in data")
      if (nSC == nrow(data))
        RSideColors = SideCols
      if (nSC == ncol(data))
        CSideColors = SideCols
    }
    if (!is.null(classvec2)) {
      out = proc.classvec(classvec2)
      nSC = out[[1]]
      SideCols = out[[2]]
      if (!nSC %in% dim(data))
        print("Error: classvec2 length not equal to nrow or ncol in data")
      if (nSC == nrow(data))
        RSideColors = SideCols
      if (nSC == ncol(data))
        CSideColors = SideCols
    }
  }
  if (all(is.null(RSideColors), is.null(CSideColors)))
    .heatmap.3(data, Colv = Colv, Rowv = Rowv, col = plotcols, labRowCol = labRowCol,
              scale = scale, trace = "none", density.info = "density",
              dendrogram = dend, ...)
  if (all(!is.null(RSideColors), is.null(CSideColors)))
    .heatmap.3(data, Colv = Colv, Rowv = Rowv, col = plotcols, labRowCol = labRowCol,
              scale = scale, trace = "none", density.info = "density",
              RowSideColors = RSideColors, dendrogram = dend, ...)
  if (all(is.null(RSideColors), !is.null(CSideColors)))
    .heatmap.3(data, Colv = Colv, Rowv = Rowv, col = plotcols, labRowCol = labRowCol,
              scale = scale, trace = "none", density.info = "density",
              ColSideColors = CSideColors, dendrogram = dend, ...)
  if (all(!is.null(RSideColors), !is.null(CSideColors)))
    .heatmap.3(data, Colv = Colv, Rowv = Rowv, col = plotcols, labRowCol = labRowCol,
              scale = scale, trace = "none", density.info = "density",
              RowSideColors = RSideColors, ColSideColors = CSideColors,
              dendrogram = dend, ...)
  if (returnSampleTree)
    return(Colv)
}


##########################
### heatmap.3 function ###
##########################

.heatmap.3 <- function(x,
                      Rowv, Colv,
                      distfun = dist,
                      hclustfun = hclust,
                      dendrogram = c("both","row", "column", "none"),
                      symm = FALSE,
                      scale = c("none","row", "column"),
                      na.rm = TRUE,
                      revC = identical(Colv,Rowv),
                      add.expr,
                      breaks,
                      symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                      col = "heat.colors",
                      colsep,
                      rowsep,
                      sepcolor = "white",
                      sepwidth = c(0.05, 0.05),
                      cellnote,
                      notecex = 1,
                      notecol = "cyan",
                      na.color = par("bg"),
                      trace = c("none", "column","row", "both"),
                      tracecol = "black",
                      hline = median(breaks),
                      vline = median(breaks),
                      linecol = tracecol,
                      margins = c(5,5),
                      ColSideColors,
                      RowSideColors,
                      side.height.fraction=0.3,
                      cexRow = 0.2 + 1/log10(nr),
                      cexCol = 0.2 + 1/log10(nc),
                      labRow = NULL,
                      labCol = NULL,
                      labRowCol = NULL,
                      key = TRUE,
                      keysize = 1.5,
                      density.info = c("none", "histogram", "density"),
                      denscol = tracecol,
                      symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      lmat = NULL,
                      lhei = NULL,
                      lwid = NULL,
                      ColSideColorsSize = 1,
                      RowSideColorsSize = 1,
                      KeyValueName="Value",...){

  invalid <- function (x) {
    if (missing(x) || is.null(x) || length(x) == 0)
      return(TRUE)
    if (is.list(x))
      return(all(sapply(x, invalid)))
    else if (is.vector(x))
      return(all(is.na(x)))
    else return(FALSE)
  }
  x <- as.matrix(x)
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale))
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col))
    col <- get(col, mode = "function")
  if (!missing(breaks) && (scale != "none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  if (is.null(Rowv) || is.na(Rowv))
    Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv))
    Colv <- FALSE
  else if (Colv == "Rowv" && !isTRUE(Rowv))
    Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2)
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                 c("both", "row"))) {
      if (is.logical(Colv) && (Colv))
        dendrogram <- "column"
      else dedrogram <- "none"
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                 c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv))
        dendrogram <- "row"
      else dendrogram <- "none"
      warning("Discrepancy: Colv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  }
  else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc)
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
      (1:nr)[rowInd]
  else rownames(x)
  else
    labRow <- labRow[rowInd]
  labRowCol <- labRowCol[rowInd]
  if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
    if (missing(col) || is.function(col))
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks)
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function")
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei))
    lhei <- c(keysize, 4)
  if (missing(lwid) || is.null(lwid))
    lwid <- c(keysize, 4)
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)

    if (!missing(ColSideColors)) {
      if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
        stop("'ColSideColors' must be a matrix of nrow(x) rows")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      lhei=c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
    }

    if (!missing(RowSideColors)) {
      if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
        stop("'RowSideColors' must be a matrix of ncol(x) columns")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
      lwid <- c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }

  if (length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))

  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)

  if (!missing(RowSideColors)) {
    if (!is.matrix(RowSideColors)){
      par(mar = c(margins[1], 0, 0, 0.5))
      image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
      text(x = 0.5, y = seq(1,nr,nr/6), labels = 1:6)
    } else {
      par(mar = c(margins[1], 0, 0, 0.5))
      rsc = t(RowSideColors[,rowInd, drop=FALSE])
      rsc.colors = matrix()
      rsc.names = names(table(rsc))
      rsc.i = 1
      for (rsc.name in rsc.names) {
        rsc.colors[rsc.i] = rsc.name
        rsc[rsc == rsc.name] = rsc.i
        rsc.i = rsc.i + 1
      }
      rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
      image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
      text(y = c(cumsum(table(rsc[,1])) - table(rsc[,1])/2)/dim(rsc)[1], x = 0, labels = length(table(rsc[,1])):1,
           font = 2, col = "white", cex = 1.2)
      if (length(rownames(RowSideColors)) > 0) {
        axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE)
      }
    }
  }

  if (!missing(ColSideColors)) {

    if (!is.matrix(ColSideColors)){
      par(mar = c(0.5, 0, 0, margins[2]))
      image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    } else {
      par(mar = c(0.5, 0, 0, margins[2]))
      csc = ColSideColors[colInd, , drop=FALSE]
      csc.colors = matrix()
      csc.names = names(table(csc))
      csc.i = 1
      for (csc.name in csc.names) {
        csc.colors[csc.i] = csc.name
        csc[csc == csc.name] = csc.i
        csc.i = csc.i + 1
      }
      csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
      image(csc, col = as.vector(csc.colors), axes = FALSE)
      if (length(colnames(ColSideColors)) > 0) {
        axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
      }
    }
  }

  par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr"))
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
  retval$carpet <- x
  if (exists("ddr"))
    retval$rowDendrogram <- ddr
  if (exists("ddc"))
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  if (!invalid(na.color) & any(is.na(x))) {
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
          col = na.color, add = TRUE)
  }
  axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
       cex.axis = cexCol)
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  if(is.null(labRowCol))
    labRowCol <- rep("black",length(labRow))
  Map(function(x,y,z)
    axis(4, at=x, col.axis=y, labels=z, lwd=0, las=1, cex.axis = cexRow), 1:length(labRow), labRowCol, labRow)
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr))
    eval(substitute(add.expr))
  if (!missing(colsep))
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)),
                              xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep),
                              lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  if (!missing(rowsep))
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1,
                              ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1,
                              lwd = 1, col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol,
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        abline(h = i + hline, col = linecol, lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote))
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
         col = notecol, cex = notecex)
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }
  else plot.new()
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  }
  else plot.new()
  if (!is.null(main))
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    }
    else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }

    z <- seq(min.raw, max.raw, length = length(col))
    par(mar = c(3,3,3,3))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv)
    if (scale == "row")
      mtext(side = 1, "Row Z-Score", line = 2)
    else if (scale == "column")
      mtext(side = 1, "Column Z-Score", line = 2)
    else mtext(side = 1, KeyValueName, line = 2)
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
            lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
      # title("Color Key\nand Density Plot")
      par(cex = 0.5)
      mtext(side = 2, "Density", line = 3, las = 3)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
            col = denscol)
      axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
      # title("Color Key\nand Histogram")
      par(cex = 0.5)
      mtext(side = 2, "Count", line = 3, las = 3)
    }
    else title("Color Key")
  }
  else plot.new()
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                  high = retval$breaks[-1], color = retval$col)
  invisible(retval)
}


##############################
##############################
## customBoxplot

.customBoxplot <- function(data, col = NA, ylab = NA, xlab = NA, ylim = NULL,
                          paired.test = FALSE, p.val = 0.05, pt.cex = 1,
                          print.name = FALSE, ...){
  if(!is.list(data))
    stop("Data must be a 'list'")

  if(all(is.na(col)))
    col <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(data))

  if(is.null(ylim) & paired.test)
    ylim <- range(data)+c(-diff(range(data))*0.1,diff(range(data))*0.3)
  else if(is.null(ylim))
    ylim <- range(data)+c(-diff(range(data))*0.1,diff(range(data))*0.1)

  bdata <- lapply(data, function(x) boxplot(x, plot = FALSE))

  plot(0, type = "n", xlab = xlab, ylab = ylab, axes = FALSE,
       xlim = c(0.5,length(data)+0.5), ylim = ylim, ...)

  for(i in 1:length(bdata))
  {
    rect(ybottom = bdata[[i]]$stats[2,], ytop = bdata[[i]]$stats[4,], xleft = i-0.25, xright = i+0.25,
         col = adjustcolor(col[i],alpha.f = 0.4), border = "black")
    segments(y0 = bdata[[i]]$stats[1,], y1 = bdata[[i]]$stats[2,], x0 = i, x1 = i, lwd = 2)
    segments(y0 = bdata[[i]]$stats[4,], y1 = bdata[[i]]$stats[5,], x0 = i, x1 = i, lwd = 2)
    segments(y0 = bdata[[i]]$stats[3,], y1 = bdata[[i]]$stats[3,], x0 = i-0.25, x1 = i+0.25, lwd = 4)
    stripchart(data[[i]], at = i, method = "jitter", pch = 21, cex = pt.cex,
               vertical = TRUE, add = TRUE, bg = adjustcolor(col[i],alpha.f = 0.3),
               col = adjustcolor("white",0.6))
    if(print.name & length(bdata[[i]]$out) > 0)
      text(data[[i]][names(bdata[[i]]$out)], x = i+0.15, col = "red",
           labels = names(bdata[[i]]$out), font = 2)
  }
  axis(1, at = 1:length(data), labels = names(data), font = 2, lwd = 0, ...)
  axis(2, ...)

  if(paired.test){
    resp <- unlist(data)
    cl <- c()
    for(j in 1:length(data))
      cl <- c(cl, rep(names(data)[j],lapply(data, length)[[j]]))
    test <- pairwise.wilcox.test(resp, cl, p.adjust.method = "fdr")
    print(test)

    pos <- which(test$p.value < p.val)
    print(pos)
    interval <- seq(diff(range(data))*0.05,diff(range(data))*0.2,length.out = length(pos))
    par(xpd = TRUE)
    for(i in 1:length(pos)){
      x1 <- which(names(data) == colnames(test$p.value)[arrayInd(pos[i], dim(test$p.value))[,2]])
      x0 <- which(names(data) == rownames(test$p.value)[arrayInd(pos[i], dim(test$p.value))[,1]])
      x <- sort(c(x1,x0))
      segments(y0 = range(data)[2]+interval[i], y1 = range(data)[2]+interval[i],
               x0 = x[1], x1 = x[2], lwd = 2)
      segments(y0 = range(data)[2]+interval[i], y1 = range(data)[2]+interval[i]-diff(range(data))*0.01,
               x0 = x[1], x1 = x[1], lwd = 2)
      segments(y0 = range(data)[2]+interval[i], y1 = range(data)[2]+interval[i]-diff(range(data))*0.01,
               x0 = x[2], x1 = x[2], lwd = 2)
    }
  }
}


##################################################################################################################################

.plotFeatureGaining <- function(h, sd, deco, print.annot, main, id = NA){

  # h <- apply(deco@NSCAcluster$All$NSCA$h, 1, function(x) sum(abs(range(x))))
  # h <- log2((range01(h)+1)/(range01(sd)+1))
  # sd <- 0.5*log2((range01(h)+1)*(range01(sd)+1))

  names(h) <- as.character(deco@featureTable[,"ID"])
  names(sd) <- as.character(deco@featureTable[,"ID"])
  # sd <- sd[names(h)]
  m <- lm(h ~ sd)

  if(print.annot & "SYMBOL"%in%colnames(deco@featureTable))
    names <- as.character(deco@featureTable[,"SYMBOL"])
  else
    names <- as.character(deco@featureTable[,"ID"])

  names(names) <- as.character(deco@featureTable[,"ID"])

  col <- sapply(m$residuals, function(x) max(which(seq(min(m$residuals),max(m$residuals),length.out = 32) <= x)))
  color <- colorRampPalette(c("blue","green","orange","red"))(32)
  col <- color[col]

  plot(sd, h, bg = adjustcolor(col,0.7), col = adjustcolor("grey",0.5), xlab = "Standard deviation", ylab = "h Range",
       pch = 21, xlim = c(0,max(sd)), ylim = c(0,max(h)), cex = 1.2)
  abline(h = 0, col = adjustcolor("grey",0.6))
  lines(smooth.spline(sd, predict(m)), lwd = 3, lty = 5, col = adjustcolor("darkblue",0.7))
  text(sd[names(sort(abs(m$residuals), decreasing = TRUE))[1:10]], h[names(sort(abs(m$residuals), decreasing = TRUE))[1:10]],
       labels = names[names(sort(abs(m$residuals), decreasing = TRUE))[1:10]], font = 2, col = "black", cex = 0.9)
  text(sd[id], h[id], labels = names[id], font = 2, col = "red", cex = 0.9)
  mtext(side = 3, text = main, font = 2, col = "deepskyblue4")
  mtext(side = 1, text = paste("Adjusted R-squared =",round(summary(m)$adj.r.squared, 3),"\n h =",
                               round(summary(m)$coefficients[2,1],3), "* sd + (",round(summary(m)$coefficients[1,1],3),")"),
        line = -1)
  title(list("Feature gaining respect to original variance.\nHigher values of 'h' remarks specificity feature-sample.",
             font = 2, cex = 1.4), outer = TRUE, line = -4)

}


##################################################################################################################################
##################################################################################################################################

.plotDECOmanhattan <- function(deco, id)
{
  deco@featureTable[,c("CHR","CHRLOC")] <- apply(deco@featureTable[,c("CHR","CHRLOC")], 2, as.numeric)
  deco@featureTable <- deco@featureTable[with(deco@featureTable, order(CHR, CHRLOC)),]

  freq <- table(deco@featureTable[,c("CHR")])
  col <- c("black","grey")[as.numeric(is.odd(deco@featureTable[,c("CHR")]))+1]

  layout(mat = matrix(c(1)))
  par(mar = c(8,8,8,8))
  plot(deco@featureTable[,"Standard.Chi.Square"], xlab = "CHROMOSOME", xaxt = "n", ylab = "Standard Chi Square",
       bg = adjustcolor(col,0.7), pch = 21, col = adjustcolor("grey",0.5), las = 2, cex = 1.3)
  pos <- freq/2 +
    c(0,cumsum(freq[-length(freq)]))
  axis(1, at = pos, labels = names(freq))
  n <- deco@featureTable[,"Standard.Chi.Square"] >= quantile(deco@featureTable[,"Standard.Chi.Square"], 0.95)
  col <- rep("black",length(n))
  col[deco@featureTable[n,"ID"] %in% id] <- "red"
  text(x = which(n), font = 2, cex = 1.1,
       y = deco@featureTable[n,"Standard.Chi.Square"]+deco@featureTable[n,"Standard.Chi.Square"]*0.05,
       labels = deco@featureTable[n,"SYMBOL"], xpd = NA, col = col)
  title(main = list("Manhattan plot showing Standard Chi Square values\nper feature and their chromosome locations", cex = 1.4, font = 2),
        outer = TRUE, line = -4)
}


##################################################################################################################################

plotDECOProfile <- function(deco, id, data, pdf.file = NA, plot.h = FALSE,
                            info.sample = NA, print.annot = TRUE,
                            cex.legend = 1.1, cex.names = 1, cex.samples = 1)
{
  # Open pdf file connection.
  if(is.na(pdf.file))
    stop("ERROR: No pdf name has been provided.")
  pdf(file = pdf.file,width = 18,height = 12,family = "Helvetica-Narrow")

  # Type of deco analysis
  if(all(!id %in% as.character(deco@featureTable[,c("ID")])))
  {
    stop("ERROR: Input IDs do not match with any DE feature.")
  }
  if(all(is.na(deco@classes)))
    analysis <- "Unsupervised"
  else if(length(deco@NSCAcluster) == 1)
    analysis <- "Multiclass"
  else
    analysis <- "Binary"

  # Formatting objects...
  if(print.annot){
    symbol <- deco@featureTable[,"SYMBOL"]
    names(symbol) <- deco@featureTable[,"ID"]
  }
  layout(mat = matrix(c(1,1,1,1,1,1,2,2,2,2,3,3,4,4,5,5,6,6,6,0), nrow = 2, byrow = TRUE))
  id_sus <- rownames(deco@NSCAcluster[[1]]$NSCA$h)
  limy <- range(na.omit(data))
  id0 <- as.character(id[id %in% id_sus])
  main <- paste("Raw data profile of:",id0,sep=" ")
  if(any(!(id %in% id_sus))){
    message("NOTE: Following IDs have not been found by DECO or they have not been included in NSCA
            (if COMPLETE profiles were removed), so 'Profile' will be not plotted for them:")
    message(paste(id[!(id %in% id_sus)], "\n"))}
  data0 <- data
  # This plot will be repeated along all IDs input.
  for(a in 1:length(id0))
  {
    par(mar=c(10, 6, 10, 6))
    data <- data0
    # Taking 'h' statistic for each ID and then clustering samples using this information.
    if(analysis == "Binary")
    {
      # Sample membership to subclasses.
      samplesSubclass <- rbind(cbind(deco@NSCAcluster$Control$samplesSubclass[order(deco@NSCAcluster$Control$samplesSubclass[,1]),]),
                               cbind(deco@NSCAcluster$Case$samplesSubclass[order(deco@NSCAcluster$Case$samplesSubclass[,1]),]))
      samplesSubclass <- cbind(Samples = rownames(samplesSubclass), Subclass = samplesSubclass[,1])
      rownames(samplesSubclass) <- 1:length(rownames(samplesSubclass))
      # Subclasses
      infoSubclass <- rbind(deco@NSCAcluster$Control$infoSubclass,deco@NSCAcluster$Case$infoSubclass)
      # Assignment of colors to all subclasses. It will be conserved along all the report.
      count <- table(sapply(rownames(infoSubclass),
                            function(x)unlist(strsplit(x,split=" Subclass"))[1]) != deco@control)
      color.cluster <- c(colorRampPalette(c("navyblue","lightblue"))(count[1]),colorRampPalette(c("orange","darkred"))(count[2]))
      names(color.cluster) <- rownames(infoSubclass)
      # h-statistic vector
      rangeH <- range(c(deco@NSCAcluster$Control$NSCA$h),c(deco@NSCAcluster$Case$NSCA$h))
      tab.info <- c(deco@NSCAcluster$Control$NSCA$h[as.character(id0)[a],],
                    deco@NSCAcluster$Case$NSCA$h[as.character(id0)[a],])
      tab.info[is.na(tab.info)] <- 0
      # Clustering NSCA inner product for one feature.
      patt <- innerProductAssign(inner = tab.info, samples = deco@classes, control = deco@control, analysis)$cl
      count <- table(patt)
      thr <- which(cumsum(count) >= table(deco@classes)[names(table(deco@classes))==deco@control])[1]
      colors <- c(colorRampPalette(c("blue","lightblue"))(thr),colorRampPalette(c("red","darkorange"))(length(count)-thr))
      names(count) <- c(paste("Pattern",1:thr,names(table(deco@classes))[1]),
                        paste("Pattern",1:(length(count)-thr),names(table(deco@classes))[2]))
      ord <- names(patt)
      bg.col <- adjustcolor(c(rep("navyblue",table(deco@classes)[names(table(deco@classes)) == deco@control]),
                              rep("darkred",table(deco@classes)[names(table(deco@classes)) != deco@control])), 0.6)
      limb <- list(range(deco@NSCAcluster$Control$rankingFeature.h[,seq(2,dim(deco@NSCAcluster$Control$infoSubclass)[1]*2,2)]),
                   range(deco@NSCAcluster$Case$rankingFeature.h[,seq(2,dim(deco@NSCAcluster$Case$infoSubclass)[1]*2,2)]))
      limb <- limb[[which(unlist(lapply(limb, diff)) == max(unlist(lapply(limb, diff))))]]
      pch <- c(rep(15,table(bg.col)[1]),rep(16,table(bg.col)[2]))
    }else
    {
      # Sample membership to subclasses.
      samplesSubclass <- as.matrix(deco@NSCAcluster$All$samplesSubclass[order(deco@NSCAcluster$All$samplesSubclass[,1]),])
      samplesSubclass <- cbind(Samples = rownames(samplesSubclass), Subclass = samplesSubclass[,1])
      rownames(samplesSubclass) <- 1:length(rownames(samplesSubclass))
      # Subclasses
      infoSubclass <- deco@NSCAcluster$All$infoSubclass
      # Assignment of colors to all subclasses. It will be conserved along all the report.
      color.cluster <- colorRampPalette(c("greenyellow","lightblue","darkgreen","darkred"))(dim(infoSubclass)[1])
      names(color.cluster) <- rownames(infoSubclass)
      # h-statistic vector
      rangeH <- range(deco@NSCAcluster$All$NSCA$h)
      tab.info <- c(deco@NSCAcluster$All$NSCA$h[as.character(id0[a]),])
      tab.info[is.na(tab.info)] <- 0
      # Clustering NSCA inner product for one feature.
      patt <- innerProductAssign(inner = tab.info, samples = deco@classes, analysis = analysis)$cl
      count <- table(patt)
      names(count) <- paste("Pattern",1:length(count))
      colors <- colorRampPalette(c("black","grey"))(length(count))
      ord <- names(patt)
      bg.col <- adjustcolor(rep("darkred",length(names(patt))),0.6)
      limb <- range(deco@NSCAcluster$All$rankingFeature.h[,seq(2,length(color.cluster)*2,2)])
      pch <- 16
    }
    # Ordering raw data using clustering information.
    data <- data[,ord]
    names(colors) <- names(count)
    count.col <- unlist(sapply(colors,function(x) rep(x,count[which(colors == x)])))

    # Plot feature profile.
    if(!plot.h)
      plot(data[as.character(id0[a]),], ylim = limy, pch=pch, xaxt="n", lwd = 0.8,
           xlim=c(1,ncol(data)), main=paste(main[a],"\nsorted by h-statistic ranking."),
           type="p", ylab = "Raw Data Signal",
           col = bg.col, xlab="", cex = cex.samples+0.5, axes = FALSE)
    else{
      interval <- unique(c(seq(rangeH[1],0,length.out = 10),seq(0,rangeH[2],length.out = 11)))
      colorH <- colorRampPalette(c("blue","lightblue","grey","tomato","red"))(20)
      bg.col <- adjustcolor(colorH[sapply(tab.info[ord], function(x)
        max(which(interval<=x)))],0.8)
      plot(data[as.character(id0[a]),], ylim = limy, pch=22, xaxt="n", lwd = 1.2,
         xlim=c(1,ncol(data)), main=main[a], type="p", ylab = "Raw Data Signal", col = adjustcolor("white",0.6),
         bg = bg.col, xlab="", cex = cex.samples+0.5, axes = FALSE)
    }

    axis(2, las = 2, cex.axis = 1.5)
    axis(side = 1,at = 1:dim(data)[2], labels = rownames(samplesSubclass[match(ord,samplesSubclass[,c("Samples")]),]),
         las=2, cex.axis = cex.names+0.5)
    m <- sum(abs(limy))
    if(analysis == "Binary"){
      thr <- dim(deco@NSCAcluster$Control$samplesSubclass)[1]+0.5
      abline(v = thr, col = "grey", lwd = 2, lty = 2)
      segments(x0 = c(1,thr+0.5), x1 = c(thr-0.5, dim(data)[2]),
               y0 = c(mean(data[as.character(id0[a]),1:thr]), mean(data[as.character(id0[a]),(thr+1):dim(data)[2]])),
               col = adjustcolor("black", 0.4), lwd = 2)
    }else
      segments(x0 = 1, x1 = dim(data)[2], y0 = mean(data[as.character(id0[a]),]),
               col = adjustcolor("black", 0.4), lwd = 2)
    if(print.annot){
      id_hgnc <- symbol[names(symbol)%in%as.character(id0)[a]]
      mtext(text = paste("HGNC symbol:",id_hgnc), line = 1)
      note <- TRUE}
    par(xpd=NA)
    # Including sample information under the expression plot.
    if(all(!is.na(info.sample))){
      infoS <- jColor(info.sample)
      cc <- infoS$col[,sort(colnames(infoS$col))]
      for(i in 1:dim(cc)[2]){
        rect(0.5:(length(data[as.character(id0[a]),])-0.5),xright = 1.5:(length(data[as.character(id0[a]),])+0.5), border = NA,
             limy[1]-(m*0.05*i+0.1*m), ytop = limy[1]-(m*0.05*(i+1)+0.1*m), col = cc[colnames(data), i])
        text(x = length(data[as.character(id0[a]),])+1, labels = colnames(cc)[i],
             y = limy[1]-(m*0.05*(i+0.5)+0.1*m), adj = c(0,0.5))
      }
    }
    par(xpd=FALSE)
    # Sample information and feature statistics will be plotted on the right side.
    title(main = list("RDA & NSCA: Single feature pattern", cex = 1.8,
                      font = 2), outer = TRUE, line = -2)
    par(mar = c(7,5,4,5))
    # 'h' statistic per Subclass of samples.
    d <- data.frame(samplesSubclass, h = tab.info[samplesSubclass[,"Samples"]])
    m <- sapply(names(color.cluster), function(x) mean(d[d[,"Subclass"] == x,"h"]))
    s <- sapply(names(color.cluster), function(x)
      sd(d[d[,"Subclass"] == x,"h"])/sqrt(length(d[d[,"Subclass"] == x,"h"])))
    par(xpd=NA)
    plot(0,xlim = c(0,length(m)), ylim = limb, type = "n", axes = FALSE,
         xlab = "", ylab = "h mean per DECO subclass")
    h <- sapply(rownames(infoSubclass), function(x) d[d[,"Subclass"] == x,"h"])
    rect(xleft = 0:(length(m)-1), xright = 1:length(m)-0.1, ybottom = 0, ytop = m,
         col = adjustcolor(color.cluster, 0.6), border = NA)
    segments(x0 = seq(0.45,length(m),1), x1 = seq(0.45,length(m),1), y0 = m-s, y1 = m+s)
    segments(x0 = seq(0.45,length(m),1)-0.1, x1 = seq(0.45,length(m),1)+0.1, y0 = m-s, y1 = m-s)
    segments(x0 = seq(0.45,length(m),1)-0.1, x1 = seq(0.45,length(m),1)+0.1, y0 = m+s, y1 = m+s)
    axis(1, labels = rownames(infoSubclass), at = seq(0.45,length(m),1),
         lty = 0, las = 2, cex.axis = 1.3)
    axis(2, las = 2, cex.axis = 1.5)

    par(xpd=FALSE)
    abline(h = 0)
    mtext(text = "h-statistic mean per subclass", side = 1, line = 10, font = 2)
    # Color code of 'h' statistic per sample.
    if(plot.h){
      par(mar = c(7,10,6,10))
      z1=matrix(1:20,nrow=1)
      y1=interval
      image(x = 0:1,y1,z1,col=adjustcolor(colorH,0.7), cex = 1.2,
            axes=FALSE,xlab="Feature's density",ylab="",main="'h' statistic\npoint's color code")
      lines(y = density(tab.info, bw = 1)$x, range01(density(tab.info, bw = 1)$y), lwd = 4, col = "cyan")
      axis(2, at = interval,
           lwd = 0, lwd.ticks = 1,
           labels = round(interval,2),
           font = 2, las = 2)
    }else
      plot(0, type = "n", axes = FALSE, xlab = "", ylab = "")
    # Color labels and summary of statistics per feature.
    par(xpd = TRUE)
    if(analysis == "Binary"){
      plot(1, type = "n", axes=FALSE, xlab="", ylab="")
      legend("center",legend = paste(c("CONTROL","CASE"),"samples"), col = "grey",
             pt.bg = adjustcolor("black",0.7), pch = c(21,22), bty = "n",
             pt.cex = 2, title = "Constrast design:\npoint's shape", xpd = TRUE,
             cex = cex.legend+0.5)}
    else if(analysis == "Multiclass"){
      plot(1, type = "n", axes=FALSE, xlab="", ylab="")
      legend("center",legend = levels(deco@classes),
             col = colorRampPalette(c("navyblue","darkorange","darkred"))(length(levels(deco@classes))),
             pch = 16, bty = "n", pt.cex = 2, title = "Constrast design:\ncolor of point", xpd = TRUE,
             cex = cex.legend+0.5)}
    else{
      plot(1, type = "n", axes=FALSE, xlab="", ylab="")
      legend("center",legend = c("No classes were defined\nfor RDA subsampling"),
             bty = "n", pt.cex = 2, xpd = TRUE, xjust = 0.5,
             cex = cex.legend+0.5)
    }

    plot(1, type = "n", axes=FALSE, xlab="", ylab="")
    if(all(!is.na(info.sample)))
      legend("center", legend = infoS$ty, col = names(infoS$ty), xpd = TRUE,
             pch = 15, cex = cex.legend+1/(length(infoS$ty))+0.5, title = "Sample info provided by user", bty = 'n')
    if(analysis == "Binary"){
      textplot(rbind(round(t(deco@featureTable[as.character(id0)[a],colnames(deco@featureTable)%in%
                                                 c("Standard.Chi.Square","overlap.Ctrl.Case","Repeats","delta.signal",
                                                   "sd.Ctrl","sd.Case","h.Range.Ctrl","h.Range.Case","Dendrogram.group")]),3),
                     Profile = as.character(deco@featureTable[as.character(id0)[a],c("Profile")]),
                     UpDw = as.character(deco@featureTable[as.character(id0)[a],c("UpDw")])), halign = "right")
    }else{
      textplot(round(t(deco@featureTable[as.character(id0)[a],c("Standard.Chi.Square","Repeats","sd","h.Range","Dendrogram.group")]),3),
               halign = "right")
    }
    par(xpd = FALSE)
  }
  remove(m)

  dev.off()
}
