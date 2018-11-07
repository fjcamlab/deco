###########################################################################################################
########## Intern functions ##########

################################################
################################################
# R summary function for 'deco' object

summary.deco <- function(object, ...)
{
  cat("Decomposing Heterogeneous Cohorts from Omic profiling: DECO\nSummary:\n")
  cat("\nAnalysis design: ")
  if(length(object@NSCAcluster) == 2){
    cat("Binary\nClasses compared:");print(table(object@classes))}
  else if (all(is.na(object@classes))){
    cat("Unsupervised\n")}
  else{
    cat("Multiclass\nClasses compared:");print(table(object@classes))}

  cat("\n")
  thr <- data.frame("RDA q.value"= object@q.val, "Minimum repeats" = object@rep.thr,
                    "Percentage of affected samples" = object@samp.perc*100, "NSCA variability" = object@NSCAcluster[[1]]$var)
  rownames(thr) <- "Thresholds"
  printCoefmat(thr, digits = 3)
  cat("\nNumber of features out of thresholds:",dim(object@featureTable)[1],"\n")
  if("Profile"%in%colnames(object@featureTable)){
    cat("Feature profile table:")
    print(table(object@featureTable$Profile))
  }
  cat("Number of samples affected:",dim(object@incidenceMatrix)[2],"\n")
  cat("Number of positive RDA comparisons:",object@pos.iter,"\n")
  cat("Number of total RDA comparisons:", round(object@featureTable[1,c("Repeats")]/object@featureTable[1,c("FR.Repeats")],digits = 0), "\n\n")

}

################################################
################################################
# R print function for 'deco' object

print.deco <- function(x, ...)
{
  cat("\nDecomposing Heterogeneous Cohorts from Omic profiling: DECO\nSummary:\n")
  cat("\nAnalysis design: ")
  if(length(x@NSCAcluster) == 2){
    cat("Binary\nClasses compared:");print(table(x@classes))}
  else if (all(is.na(x@classes))){
    cat("Unsupervised\n")}
  else{
    cat("Multiclass\nClasses compared:");print(table(x@classes))}

  cat("\n")
  thr <- data.frame("RDA q.value"= x@q.val, "Minimum repeats" = x@rep.thr,
                    "Percentage of affected samples" = x@samp.perc*100, "NSCA variability" = x@NSCAcluster[[1]]$var)
  rownames(thr) <- "Thresholds"
  printCoefmat(thr, digits = 3)
  cat("\nNumber of features out of thresholds:",dim(x@featureTable)[1],"\n")
  cat("Number of samples affected:",dim(x@incidenceMatrix)[2],"\n")
  cat("Number of positive RDA comparisons:",x@pos.iter,"\n")
  cat("Number of total RDA comparisons:", round(x@featureTable[1,c("Repeats")]/x@featureTable[1,c("FR.Repeats")],digits = 0), "\n\n")
  cat("RDA call:\n")
  print(x@subsampling.call)
  cat("NSCA call:\n")
  print(x@deco.call)
}

################################################
################################################
# localMinima and localMaxima

localMinima <- function(x) {
  y <- diff(c(.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {y <- y[-1]}

  return(y)
}

localMaxima <- function(x) {
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {y <- y[-1]}

  return(y)
}

################################################
################################################
# Calculate quantile threshold per standard deviation

quantSD <- function(sdx, p = NULL)
{
  if(is.null(p))
    p <- 0.75
  if(quantile(sdx, probs = p) <= median(range(sdx)))
    thr <- quantile(sdx, probs = p)
  else
    thr <- median(sdx)
  return(thr)
}

################################################
################################################
# AnnotateDECO

AnnotateDECO <- function(ids, id.type, attributes = NA, pack.db = "org.Hs.eg.db",
                         rm.redundant = TRUE, verbose = FALSE)
{
  options(warn = -1)

  if(!exists(pack.db)) stop(paste("ERROR: Annotation package '",pack.db,"' is not in library."))
  else requireNamespace(pack.db, quietly = TRUE)
  if(verbose)
    message(paste(format(Sys.time(), " %H:%M:%S"),"-- Annotating IDs with Bioconductor..."))
  if(is.na(attributes))
    attributes <- c("SYMBOL","ENSEMBL","ENTREZID","CHR","CHRLOC","GENENAME")
  if(id.type == "PROBEID"){
    xx <- as.list(get(paste(unlist(strsplit(pack.db,split = ".db")),"ENSEMBL",sep="")))
    infogenes <- suppressMessages(AnnotationDbi::select(x = get(pack.db), obj = names(xx),
                                                        columns = c(attributes,id.type),
                                                        keys = keys(get(pack.db), keytype = "ENSEMBL"),
                                                        keytype = "ENSEMBL"))
  }else
    infogenes <- suppressMessages(AnnotationDbi::select(x = get(pack.db), obj = ids, columns = attributes,
                                                        keys = keys(get(pack.db), keytype = id.type),
                                                        keytype = id.type))
  infogenes <- infogenes[apply(infogenes[,attributes], 1, function(x) all(!is.na(x))),]
  if(all(is.na(infogenes)))
    stop("Annotation library does not found input IDs.")
  infogenes <- infogenes[order(infogenes[,id.type]),]
  infogenes <- infogenes[infogenes[,id.type]%in%ids,]
  if(rm.redundant){
    infogenes <- infogenes[!(duplicated(infogenes[,id.type])),]
    rownames(infogenes) <- infogenes[,id.type]
    if(verbose)
      message(paste(format(Sys.time(), " %H:%M:%S"),"-- Removed duplicated info from IDs."))
  }

  options(warn = 0)
  return(infogenes)
}

################################################
################################################
# Sample dependence

sampleDependence <- function(results, r, s)
{
  mm <- matrix(data = 0, ncol = s, nrow = s, dimnames = list(1:s,1:s))
  if(table(results[,"nFeatures"] > 0) > 0){
    for(i in 1:length(which(results[,"nFeatures"] > 0)))
    {
      mm[results[i,1:r],results[i,(r+1):(2*r)]] <- mm[results[i,1:r],results[i,(r+1):(2*r)]]+results[i,"nFeatures"]
      mm[results[i,1:r],results[i,1:r]] <- mm[results[i,1:r],results[i,1:r]]-results[i,"nFeatures"]
      mm[results[i,(r+1):(2*r)],results[i,(r+1):(2*r)]] <- mm[results[i,(r+1):(2*r)],results[i,(r+1):(2*r)]]-results[i,"nFeatures"]
    }
  }
  mm[lower.tri(mm, diag = TRUE)] <- 0
  return(mm)
}


##############################
##############################
## innerProductAssign

innerProductAssign <- function(inner, samples = NA, control = NA, analysis)
{
  samples <- samples[names(inner)]
  if(analysis == "Binary")
  {
    res <- 0
    inner0 <- inner
    for(j in 1:length(table(samples)))
    {
      if(j == 1)
        inner <- inner0[names(samples[samples==control])]
      else
        inner <- inner0[names(samples[samples!=control])]
      error=try(expr = {d <- density.lf(x = as.matrix(inner))}, silent = TRUE)
      if(inherits(error, "try-error"))
        d <- density(x = inner)
      x <- d$x[localMaxima(d$y)]
      error=try(expr = {p <- kmeans(inner, centers = x)}, silent = TRUE)
      if(inherits(error, "try-error"))
        error=try(expr = {p <- kmeans(inner, centers = length(x)-1, nstart = 10)}, silent = TRUE)
      res <- c(res, p$cluster+max(res))
    }
    res <- res[2:length(res)]
    res <- res[order(res, inner0)]
  }else
  {
    error=try(expr = {d <- density.lf(x = as.matrix(inner))}, silent = TRUE)
    if(inherits(error, "try-error"))
      d <- density(x = inner)
    x <- d$x[localMaxima(d$y)]
    error=try(expr = {p <- kmeans(inner, centers = x)}, silent = TRUE)
    if(inherits(error, "try-error"))
      error=try(expr = {p <- kmeans(inner, centers = length(x), nstart = 10)}, silent = TRUE)
    res <- p$cluster[order(p$cluster, inner)]
  }
  return(list(cl=res,centroid=x))
}

##############################
##############################
## Distance function with correlation

distFunc <- function(x, use = "pairwise.complete.obs", cor.method = "pearson") {
  co.x <- cor(x, use = use, method = cor.method)
  dist.co.x <- 1 - co.x
  return(as.dist(dist.co.x))
}


##############################
##############################
## Cophenetic correlation and clustering dendogram

cophDECO <- function(data, method.heatmap = "ward.D", k = NULL, scale = FALSE,
                     verbose = TRUE, coph = TRUE)
{
  if(is.null(method.heatmap))
    method.heatmap <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")

  if(scale)
    data <- scale(data)

  d <- distFunc(x = as.matrix(data), cor.method = "pearson")
  d[is.na(d)] <- 0

  if(coph)
    hmet <- sapply(method.heatmap, function(x){
      sampleTree <- as.dendrogram(hclust(d, method = x))
      coph <- cor(c(d), c(cophenetic(sampleTree)),
                  method = "pearson")
    })
  else{
    hmet <- rep(1,length(method.heatmap))
    names(hmet) <- method.heatmap
  }
  sampleTree <- hclust(d, method = names(hmet[order(hmet, decreasing = TRUE)])[1])
  sampleTree$height <- sampleTree$height[order(sampleTree$height)]

  if(is.null(k)){
    p <- c()
    for(i in 2:(dim(data)[2]-1))
    {
      cutt <- cutree(sampleTree, k = i)
      p[i-1] <- .cluster.stat(d = d, clustering = cutt)$pearsongamma
    }
    p <- unlist(na.omit(p))
    cutt <- cutree(sampleTree, k = which(p == max(p)) + 1)}
  else{
    cutt <- cutree(sampleTree, k = k)
    p <- .cluster.stat(d = d, clustering = cutt)$pearsongamma
    if(verbose)
      message(paste(format(Sys.time(), " %H:%M:%S"),"-- NOTE: Number of 'k' subclasses defined by user have been settle down."))
  }
  return(list(dend = sampleTree, coph = hmet, cluster = cutt, huber = max(p)))
}


##############################
##############################
## probability subsampling

subsamplingProb <- function(x, iter, n1, n2 = 0)
{
  if(n2 > 0)
  {
    p1o <- x/n1
    p2o <- x/n2

    m <- p1o * p2o

    r1 <- factorial(n1)/(factorial(x)*factorial(n1-x))
    r2 <- factorial(n2)/(factorial(x)*factorial(n2-x))

    m1 <- m * p2o * iter
    m2 <- m * p1o * iter

    m12 <- m * (p2o*p1o)^(x-1)

    all <- r1 * r2
  }else
  {
    p1o <- x/n1
    r1 <- factorial(n1)/(factorial(x)*factorial(n1-x))
    m = 0
    all = r1
  }
  return(list(prob.s = m, prob.comb = iter/all, n.all = all))
}


##############################
##############################
## Overlapping omic signal

overlapFeature <- function(id, data, classes, control, analysis, infoS = NA, plot = FALSE){

  if(analysis == "Binary"){
    da <- density.lf(as.matrix(data[id,names(classes[classes == control])]), from = range(data)[1],
                     to = range(data)[2], n = 1000, width = 1)
    db <- density.lf(as.matrix(data[id,names(classes[classes != control])]), from = range(data)[1],
                     to = range(data)[2], n = 1000, width = 1)
    d <- data.frame(x=da$x, a=da$y, b=db$y)
    # calculate intersection densities
    d$w <- pmin(d$a, d$b)
    # integrate areas under curves
    total <- integrate.xy(d$x, d$a) + integrate.xy(d$x, d$b)
    intersection <- integrate.xy(d$x, d$w)
    col <- adjustcolor(c("navyblue","darkred"),alpha.f = 0.3)
    pos <- 2:dim(d)[2]
    txt <- c("Control", "Case")
  }
  if(analysis == "Multiclass")
  {
    d <- sapply(levels(classes), function(x) {
      error=try(expr = {d <- density.lf(as.matrix(data[id,names(classes[classes == x])]),
                                        from = range(data)[1], to = range(data)[2], n = 1000, width = 1)$y}, silent = TRUE)
      if(inherits(error, "try-error"))
        d <- density(as.numeric(data[id,names(classes[classes == x])]))$y
      return(d)
      })
    error=try(expr = {x <- density.lf(as.matrix(data[id,names(classes[classes == levels(classes)[1]])]),
                                      from = range(data)[1], to = range(data)[2], n = 1000, width = 1)$x}, silent = TRUE)
    if(inherits(error, "try-error"))
      x <- density(as.numeric(data[id,names(classes[classes == levels(classes)[1]])]))$x
    w <- apply(apply(expand.grid(levels(classes), levels(classes)), 1, function(y) pmin(d[,y[1]], d[,y[2]])), 1, max)
    d <- data.frame(d, x = x, w = w)
    total <- sum(sapply(1:length(levels(classes)), function(y) integrate.xy(d$x, d[,y])))^2
    intersection <- integrate.xy(d$x, d$w)

    pos <- 1:length(levels(classes))
    col <- adjustcolor(colorRampPalette(c("navyblue","darkorange","darkred","darkgreen"))(length(levels(classes))),alpha.f = 0.3)
    txt <- levels(classes)
  }
  # compute overlap coefficient
  overlap <- length(levels(classes)) * intersection / total

  if(plot){
    par(mar = c(10,6,10,8))
    stackpoly.2(x = d[,pos], col = col, border = "grey", axis1 = TRUE,
                stack = FALSE, lwd = 1, lty = 3, xaxlab = round(d[seq(1,dim(d)[1],dim(d)[1]*0.05),"x"],1),
                xat = seq(1,dim(d)[1],dim(d)[1]*0.05), xlab = "Raw Data Signal", ylab = "Density",
                main = paste("Overlap:",round(overlap,3),"\nLikelihood density estimation, bandwidth = 1"))
    legend("topright", legend = txt, pch = 15, cex = 1.5, pt.cex = 1.5,
           col = col, bty = "n")
  }

  return(overlap)
}

##############################
##############################
## Colors assignment

jColor <- function(info)
{
  info <- as.matrix(info)
  nam <- rownames(info)
  double <- c("ivory2","gold","gray30","cornflowerblue","chocolate","honeydew1")

  colnames(info) <- paste(rev(LETTERS[1:dim(info)[2]]),": ",colnames(info), sep="")
  for(j in 1:dim(info)[2])
    info[,j] <- paste(rev(LETTERS[1:dim(info)[2]])[j],": ",as.character(info[,j]), sep = "")
  rownames(info) <- nam
  info.sample.color <- info
  info.sample.color <- apply(info.sample.color,2,as.character)
  if(dim(info)[2] > 1){
    pos <- apply(info,2,unique)
    if(is.matrix(pos))
      pos <- rep(dim(pos)[1] == 2,dim(pos)[2])
    else
      pos <- lapply(pos, length) == 2
  }
  else
    pos <- dim(apply(info,2,unique))[1] <= 2
  ty <- sort(unlist(apply(info[,!pos, drop = FALSE],2,unique)))
  tyy <- sort(unlist(apply(info[,pos, drop = FALSE],2,unique)))
  myPalette <- c(brewer.pal(name = "Set1",n = 9)[1:6],"white",brewer.pal(name = "Set1",n = 9)[8:9],"black")
  if(length(ty) > 10)
    myPalette <- colorRampPalette(myPalette)(length(unlist(apply(info,2,unique))))
  if(length(which(pos)) < dim(info)[2])
    for(z in 1:length(ty))
    {
      info.sample.color[info %in% ty[z]] <-
        as.character(rep(myPalette[z],length(info.sample.color[info == unlist(apply(info,2,unique))[z]])))
      names(ty)[z] <- unique(as.character(rep(myPalette[z],length(info.sample.color[info == unlist(apply(info,2,unique))[z]]))))
    }
  rownames(info.sample.color) <- nam
  ty <- c(tyy, ty)
  if(any(pos)){
    for(j in 1:(length(which(pos))*2))
      info.sample.color[info == ty[j]] <- rep(double,10)[j]
    names(ty)[1:j] <- rep(double,10)[1:j]
  }

  return(list(orig = info, col = info.sample.color, ty = ty))
}

##############################
##############################
## RNAseqFilter

RNAseqFilter <- function(data, q = 0.95, thr = 1)
{
  # thr <- density(data)$x[localMaxima(density(data)$y)][2]

  dataF <- data[apply(data, 1, function(x) length(which(x < thr))/length(x)) >= q,]
  dataNoF <- data[apply(data, 1, function(x) length(which(x < thr))/length(x)) < q,]

  message(paste(format(Sys.time(), " %H:%M:%S")," NOTE: Number of RNAseq features filtered is", dim(dataF)[1],
                ",\n",dim(dataNoF)[1],"features have been selected for subsampling."))

  return(dataNoF)
}

##############################
##############################
## calchDi

calchDi <- function(deco, samplesSubclass)
{
  n <- table(samplesSubclass[,c("Subclass")])
  sapply(1:n, function(x) deco@featureTable$Dendrogram.group)
}

###########################
### extractBestFeatures ###
###########################
# BestFeatures per subclass
# Returns only a character vector with IDs

bestFeatures <- function(nsca, data, f = 5){

  d <- t(apply(nsca$NSCA$h, 1, function(x) {
    p <- innerProductAssign(x)
    match(colnames(nsca$NSCA$h),names(p$cl))
  }))

  colnames(d) <- colnames(nsca$NSCA$h)

  infoSubclass <- nsca$infoSubclass
  samplesSubclass <- nsca$samplesSubclass

  r <- sapply(rownames(infoSubclass), function(y)
    apply(d[,rownames(samplesSubclass)[samplesSubclass[,c("Subclass")]%in%y]],1,function(x)
      mean(dist(x))))

  h <- sapply(rownames(infoSubclass), function(y)
    apply(nsca$NSCA$h[,rownames(samplesSubclass)[samplesSubclass[,c("Subclass")]%in%y]],1,function(x)
      mean(x)))

  for(i in 1:dim(infoSubclass)[1])
    r[,i] <- r[,i]/((infoSubclass[i,"Samples"]+1)/3)

  res <- h * r

  # l <- perc/100 * dim(res)[1]

  top1 <- apply(res, 2, function(x) names(sort(x, decreasing = TRUE)[c(1:f,(length(x)-f+1):length(x))]))
  top2 <- apply(res, 2, function(x) sort(x, decreasing = TRUE)[c(1:f,(length(x)-f+1):length(x))])

  return(list(topN = top1, topV = top2, h = h, r = r))
}

#################################
#################################
###

NSCAcluster <- function(mx, data = NULL, id.names = NULL, k = NULL,
                        label = NA, v = 80, method.dend = NULL,
                        raw = NULL, UpDw = NULL, dir = NA)
{
  options(warn = -1)

  if(is.null(id.names))
    id.names <- rownames(mx)

  ca_res <- NSCA(as.matrix(mx[apply(mx,1,sum) > 0, apply(mx,2,sum) > 0]), v = v)
  nd <- min(which(ca_res$Inertia[,3] >= v))
  if(nd == 1){nd <- 2; message(paste(format(Sys.time(), "\r %H:%M:%S"),"-- 1D is enough to reach explained variability from input. 2D will be considered"))}
  ca_ev <- rbind(ca_res$fbip[,1:nd],ca_res$g[,1:nd])
  # Renaming feature IDs
  if(all(!is.null(id.names)))
    rownames(ca_ev)[rownames(ca_ev) %in% names(id.names)] <- as.character(id.names[
      rownames(ca_ev)[rownames(ca_ev) %in% names(id.names)]])

  # Calculating h statistic
  if(all(!is.null(raw)) & all(!is.null(data))){
    dispersion <- data[id.names[rownames(ca_res$inner.prod)],colnames(ca_res$inner.prod)]-
      raw[id.names[rownames(ca_res$inner.prod)]]
    # Dispersion inversion for DOWN regulated within a category of samples
    if(all(!is.null(UpDw)))
      dispersion[rownames(dispersion) %in% id.names[id.names %in% names(UpDw[UpDw == dir])],] <- -dispersion[
        rownames(dispersion) %in% id.names[id.names %in% names(UpDw[UpDw == dir])],]
    inner <- as.matrix(ca_res$inner.prod)
    inner[] <- scale(as.vector(inner))
    dispersion[] <- scale(as.vector(dispersion))
    # Make positive ranges
    inner <- inner + abs(min(inner)) + 1
    rownames(dispersion) <- rownames(inner)
    ca_res$inner.prod <-  as.matrix(inner * dispersion - mean(as.matrix(inner * dispersion)))
    # 'h' inversion for DOWN regulated within a category of samples
    if(all(!is.null(UpDw)))
      ca_res$inner.prod[rownames(ca_res$inner.prod) %in% names(id.names)[id.names %in% names(UpDw[UpDw == dir])],] <- -ca_res$inner.prod[
        rownames(ca_res$inner.prod) %in% names(id.names)[id.names %in% names(UpDw[UpDw == dir])],]
  }else if(all(!is.null(data))){
    inner <- as.matrix(ca_res$inner.prod)
    inner[] <- scale(as.vector(inner))
    dispersion <- data[id.names[rownames(ca_res$inner.prod)],colnames(ca_res$inner.prod)]-
      apply(data[id.names[rownames(ca_res$inner.prod)],colnames(ca_res$inner.prod)],1,mean)
    dispersion[] <- scale(as.vector(dispersion))
    # Make positive ranges
    inner <- inner + abs(min(inner)) + 1
    # dispersion <- dispersion + abs(min(dispersion)) + 1
    ca_res$inner.prod <-  as.matrix(inner * dispersion - mean(as.matrix(inner * dispersion)))
  }
  else
    message("Data has not been provided, then 'h statistic' will be not computed. Inner product from NSCA\n
            would be used in biclustering of features and samples.")
  # Make sure assignment
  rownames(ca_res$inner.prod) <- rownames(inner)

  # Biclustering
  info.dend.samp <- cophDECO(data = as.matrix(ca_res$inner.prod), method.heatmap = method.dend,
                             k = k, scale = FALSE, coph = TRUE)
  info.dend.feat <- cophDECO(data = as.matrix(t(ca_res$inner.prod)), method.heatmap = "ward.D",
                             k = 10, scale = FALSE, verbose = FALSE, coph = FALSE)
  k <- max(info.dend.samp$cluster)

  # Calculating some statistics...
  if(dim(mx)[2] > 5)
    distf = mean
  else
    distf = median
  kk <- vector(length=k)
  gg1 <- vector(length=k)
  gg2 <- vector(length=k)
  ss <- vector(length=k)
  message(paste(format(Sys.time(), "\r %H:%M:%S"),"-- Optimized for",k,"subclasses."))
  info <- sapply(1:k, function(y) apply(cbind(as.matrix(ca_res$inner.prod)[, info.dend.samp$cluster==y]), 1, distf))
  info <- cbind(info, ca_res$di[rownames(info)])
  if(k > 1)
    info <- as.data.frame(cbind(info,matrix(nrow = dim(info)[1], ncol = 3)))
  else
    info <- as.data.frame(t(rbind(info,matrix(nrow = dim(info)[1], ncol = 3))))
  colnames(info) <- c(paste("Scl",1:k,sep=""),"Tau.feature","Closer.subclass","h.Best","ID")
  info <- info[order(rownames(info),decreasing = TRUE),]
  if(k>1)
    info[,c("Closer.subclass")] <- apply(info[,1:k],1,function(x) which(abs(x) == max(abs(x))))
  else
    info[,c("Closer.subclass")] <- rep(1, dim(info)[1])

  for(i in 1:k)
  {
    kk[i] <- length(which(colnames(mx) %in% names(info.dend.samp$cluster[which(info.dend.samp$cluster==i)])))
    gg1[i] <- length(which(info[which(sapply(rownames(info), function(x) unlist(
      strsplit(x,split = "deco",fixed = TRUE))[2])=="UP"),"Closer.subclass"] == i))
    gg2[i] <- length(which(info[which(sapply(rownames(info), function(x) unlist(
      strsplit(x,split = "deco",fixed = TRUE))[2])=="DOWN"),"Closer.subclass"] == i))
    p <- abs(info[info[,"Closer.subclass"] == i, i])
    ss[i] <- mean(p[order(p, decreasing = TRUE)][1:(length(p)*0.05)])
  }
  if(all(!is.na(label)))
    infoSubclass <- data.frame(Samples=kk, FeaturesUP=gg1, FeaturesDOWN=gg2, SpeValue = ss,
                               row.names = paste(label,"Subclass",seq(from=1,to=k),sep=" "))
  else
    infoSubclass <- data.frame(Samples=kk, FeaturesUP=gg1, FeaturesDOWN=gg2, SpeValue = ss,
                               row.names = paste("Subclass",seq(from=1,to=k),sep=" "))
  if(k>1)
    info[,"h.Best"] <- apply(info[,c(paste("Scl",1:k,sep=""),"Closer.subclass")],1,function(x) x[x[k+1]])
  else
    info[,"h.Best"] <- info[,1]

  # Renaming some feature IDs and removing duplicated features.
  info[,c("ID")] <- as.character(id.names[rownames(info)])
  rownames(ca_ev)[rownames(ca_ev)%in%rownames(mx)][!(
    duplicated(id.names[rownames(ca_ev[rownames(ca_ev)%in%rownames(mx),])]) | duplicated(id.names[
      rownames(ca_ev[rownames(ca_ev)%in%rownames(mx),])], fromLast = TRUE))] <- id.names[rownames(ca_ev[rownames(ca_ev)%in%rownames(mx),])][!(
        duplicated(id.names[rownames(ca_ev[rownames(ca_ev)%in%rownames(mx),])]) | duplicated(id.names[rownames(ca_ev[rownames(ca_ev)%in%rownames(mx),])], fromLast = TRUE))]
  info <- info[!rownames(info)%in%rownames(info[info[,c("ID")]%in%info[duplicated(info[,c("ID")]),c("ID")],])[
    which(unlist(strsplit(rownames(info[info[,c("ID")]%in%info[duplicated(info[,c("ID")]),c("ID")],]),split = "deco",fixed = TRUE))=="DOWN")/2],]
  rownames(ca_res$inner.prod)[!(duplicated(
    id.names[rownames(ca_res$inner.prod)]) | duplicated(id.names[rownames(ca_res$inner.prod)], fromLast = TRUE))] <- id.names[
      rownames(ca_res$inner.prod)][!(duplicated(id.names[rownames(ca_res$inner.prod)]) | duplicated(id.names[
        rownames(ca_res$inner.prod)], fromLast = TRUE))]

  # Calculating h statistic per feature per subclass of samples.
  ord <- 1:length(info.dend.feat$cluster)
  names(ord) <- rownames(ca_res$inner.prod)[rev(info.dend.feat$dend$order)]
  patt <- info.dend.feat$cluster
  for(i in 1:max(info.dend.feat$cluster))
    patt[info.dend.feat$cluster == i] <- which(order(sapply(1:max(info.dend.feat$cluster), function(x)
      mean(which(info.dend.feat$cluster[rev(info.dend.feat$dend$order)] == x)))) == i)
  info.dend.feat$cluster <- patt
  info <- data.frame(info, h.Range = apply(info[,1:k],1,function(x) diff(range(x))),
                     Dendrogram.group = info.dend.feat$cluster[rownames(info)], Dendrogram.order = ord[rownames(info)])
  rownames(info) <- as.character(info[,c("ID")])
  rankingH <- cbind(t(interleave(t(apply(info[,1:k],2,function(x) rank(-abs(x), ties.method = "max"))),
                                 t(info[,1:k]))),info[,"h.Range"],
                    info[,"Dendrogram.group"], ord[rownames(info)])
  colnames(rankingH) <- c(paste(rep(c("Ranking","h"),k),rep(paste("Scl",1:k,sep = ""), each = 2),sep = "."),
                          "h.Range","Dendrogram.group","Dendrogram.order")

  # Find out sample membership to any subclass.
  samplesSubclass <- cbind(info.dend.samp$cluster)
  colnames(samplesSubclass) <- "Subclass"
  for(i in 1:dim(samplesSubclass)[1])
    samplesSubclass[i,1] <- rownames(infoSubclass)[as.numeric(samplesSubclass[i,1])]

  # Returning results.
  if(all(!is.null(data))){
    names(ca_res)[names(ca_res) == "inner.prod"] <- "h"
    results <- list(ca_res, info, rankingH, infoSubclass, mx, ca_ev, samplesSubclass, ca_res$Inertia[nd,3],
                    method.dend, info.dend.samp, info.dend.feat, inner, dispersion)
    names(results) <- c("NSCA", "infoFeature", "rankingFeature.h","infoSubclass", "incidenceMatrix",
                        "Biplot.coordinates", "samplesSubclass", "var",
                        "clus.method", "hclustSamp", "hclustFeat",
                        "inner","dispersion")
  }else{
    results <- list(ca_res, info, rankingH, infoSubclass, mx, ca_ev, samplesSubclass, ca_res$Inertia[nd,3],
                    method.dend, info.dend.samp, info.dend.feat)
    names(results) <- c("NSCA", "infoFeature", "rankingFeature.h","infoSubclass", "incidenceMatrix",
                        "Biplot.coordinates", "samplesSubclass", "var",
                        "clus.method", "hclustSamp", "hclustFeat")
  }

  names(results$var) <- "Variability explained by NSCA"
  message(paste(format(Sys.time(), "\r %H:%M:%S"),"-- Summary of clustering analysis:"))

  print(infoSubclass)
  print(results$var)

  return(results)
}


##############################
##############################
###
.cluster.stat <- function (d = NULL, clustering, alt.clustering = NULL, noisecluster = FALSE,
                          silhouette = TRUE, G2 = FALSE, G3 = FALSE, wgap = TRUE, sepindex = TRUE,
                          sepprob = 0.1, sepwithnoise = TRUE, compareonly = FALSE,
                          aggregateonly = FALSE)
{
  if (!is.null(d))
    d <- as.dist(d)
  cn <- max(clustering)
  clusteringf <- as.factor(clustering)
  clusteringl <- levels(clusteringf)
  cnn <- length(clusteringl)
  if (cn != cnn) {
    warning("clustering renumbered because maximum != number of clusters")
    for (i in 1:cnn) clustering[clusteringf == clusteringl[i]] <- i
    cn <- cnn
  }
  n <- length(clustering)
  noisen <- 0
  cwn <- cn
  if (noisecluster) {
    noisen <- sum(clustering == cn)
    cwn <- cn - 1
  }
  diameter <- average.distance <- median.distance <- separation <- average.toother <- cluster.size <- within.dist <- between.dist <- numeric(0)
  for (i in 1:cn) cluster.size[i] <- sum(clustering == i)
  pk1 <- cluster.size/n
  pk10 <- pk1[pk1 > 0]
  h1 <- -sum(pk10 * log(pk10))
  corrected.rand <- vi <- NULL
  if (!is.null(alt.clustering)) {
    choose2 <- function(v) {
      out <- numeric(0)
      for (i in 1:length(v)) out[i] <- ifelse(v[i] >= 2,
                                              choose(v[i], 2), 0)
      out
    }
    cn2 <- max(alt.clustering)
    clusteringf <- as.factor(alt.clustering)
    clusteringl <- levels(clusteringf)
    cnn2 <- length(clusteringl)
    if (cn2 != cnn2) {
      warning("alt.clustering renumbered because maximum != number of clusters")
      for (i in 1:cnn2) alt.clustering[clusteringf == clusteringl[i]] <- i
      cn2 <- cnn2
    }
    nij <- table(clustering, alt.clustering)
    dsum <- sum(choose2(nij))
    cs2 <- numeric(0)
    for (i in 1:cn2) cs2[i] <- sum(alt.clustering == i)
    sum1 <- sum(choose2(cluster.size))
    sum2 <- sum(choose2(cs2))
    pk2 <- cs2/n
    pk12 <- nij/n
    corrected.rand <- (dsum - sum1 * sum2/choose2(n))/((sum1 +
                                                          sum2)/2 - sum1 * sum2/choose2(n))
    pk20 <- pk2[pk2 > 0]
    h2 <- -sum(pk20 * log(pk20))
    icc <- 0
    for (i in 1:cn) for (j in 1:cn2) if (pk12[i, j] > 0)
      icc <- icc + pk12[i, j] * log(pk12[i, j]/(pk1[i] *
                                                  pk2[j]))
    vi <- h1 + h2 - 2 * icc
  }
  if (compareonly) {
    out <- list(corrected.rand = corrected.rand, vi = vi)
  }
  else {
    dmat <- as.matrix(d)
    within.cluster.ss <- 0
    overall.ss <- nonnoise.ss <- sum(d^2)/n
    if (noisecluster)
      nonnoise.ss <- sum(as.dist(dmat[clustering <= cwn,
                                      clustering <= cwn])^2)/sum(clustering <= cwn)
    ave.between.matrix <- separation.matrix <- matrix(0,
                                                      ncol = cn, nrow = cn)
    di <- list()
    for (i in 1:cn) {
      cluster.size[i] <- sum(clustering == i)
      di <- as.dist(dmat[clustering == i, clustering ==
                           i])
      if (i <= cwn) {
        within.cluster.ss <- within.cluster.ss + sum(di^2)/cluster.size[i]
        within.dist <- c(within.dist, di)
      }
      if (length(di) > 0)
        diameter[i] <- max(di)
      else diameter[i] <- NA
      average.distance[i] <- mean(di)
      median.distance[i] <- median(di)
      bv <- numeric(0)
      for (j in 1:cn) {
        if (j != i) {
          sij <- dmat[clustering == i, clustering ==
                        j]
          bv <- c(bv, sij)
          if (i < j) {
            separation.matrix[i, j] <- separation.matrix[j,
                                                         i] <- min(sij)
            ave.between.matrix[i, j] <- ave.between.matrix[j,
                                                           i] <- mean(sij)
            if (i <= cwn & j <= cwn)
              between.dist <- c(between.dist, sij)
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
    ch <- between.cluster.ss * (n - noisen - cwn)/(within.cluster.ss *
                                                     (cwn - 1))
    clus.avg.widths <- avg.width <- NULL
    if (silhouette) {
      sii <- silhouette(clustering, dmatrix = dmat)
      sc <- summary(sii)
      clus.avg.widths <- sc$clus.avg.widths
      if (noisecluster)
        avg.width <- mean(sii[clustering <= cwn, 3])
      else avg.width <- sc$avg.width
    }
    g2 <- g3 <- cn2 <- cwidegap <- widestgap <- sindex <- NULL
    if (G2) {
      splus <- sminus <- 0
      for (i in 1:nwithin) {
        splus <- splus + sum(within.dist[i] < between.dist)
        sminus <- sminus + sum(within.dist[i] > between.dist)
      }
      g2 <- (splus - sminus)/(splus + sminus)
    }
    if (G3) {
      sdist <- sort(c(within.dist, between.dist))
      sr <- nwithin + nbetween
      dmin <- sum(sdist[1:nwithin])
      dmax <- sum(sdist[(sr - nwithin + 1):sr])
      g3 <- (sum(within.dist) - dmin)/(dmax - dmin)
    }
    pearsongamma <- cor(c(within.dist, between.dist), c(rep(0,
                                                            nwithin), rep(1, nbetween)))
    dunn <- min(separation[1:cwn])/max(diameter[1:cwn], na.rm = TRUE)
    acwn <- ave.between.matrix[1:cwn, 1:cwn]
    dunn2 <- min(acwn[upper.tri(acwn)])/max(average.distance[1:cwn],
                                            na.rm = TRUE)
    if (wgap) {
      cwidegap <- rep(0, cwn)
      for (i in 1:cwn) if (sum(clustering == i) > 1)
        cwidegap[i] <- max(hclust(as.dist(dmat[clustering ==
                                                 i, clustering == i]), method = "single")$height)
      widestgap <- max(cwidegap)
    }
    if (sepindex) {
      psep <- rep(NA, n)
      if (sepwithnoise | !noisecluster) {
        for (i in 1:n) psep[i] <- min(dmat[i, clustering !=
                                             clustering[i]])
        minsep <- floor(n * sepprob)
      }
      else {
        dmatnn <- dmat[clustering <= cwn, clustering <=
                         cwn]
        clusteringnn <- clustering[clustering <= cwn]
        for (i in 1:(n - noisen)) psep[i] <- min(dmatnn[i,
                                                        clusteringnn != clusteringnn[i]])
        minsep <- floor((n - noisen) * sepprob)
      }
      sindex <- mean(sort(psep)[1:minsep])
    }
    if (!aggregateonly)
      out <- list(n = n, cluster.number = cn, cluster.size = cluster.size,
                  min.cluster.size = min(cluster.size[1:cwn]),
                  noisen = noisen, diameter = diameter, average.distance = average.distance,
                  median.distance = median.distance, separation = separation,
                  average.toother = average.toother, separation.matrix = separation.matrix,
                  ave.between.matrix = ave.between.matrix, average.between = average.between,
                  average.within = average.within, n.between = nbetween,
                  n.within = nwithin, max.diameter = max(diameter[1:cwn],
                                                         na.rm = TRUE), min.separation = sepwithnoise *
                    min(separation) + (!sepwithnoise) * min(separation[1:cwn]),
                  within.cluster.ss = within.cluster.ss, clus.avg.silwidths = clus.avg.widths,
                  avg.silwidth = avg.width, g2 = g2, g3 = g3, pearsongamma = pearsongamma,
                  dunn = dunn, dunn2 = dunn2, entropy = h1, wb.ratio = average.within/average.between,
                  ch = ch, cwidegap = cwidegap, widestgap = widestgap,
                  sindex = sindex, corrected.rand = corrected.rand,
                  vi = vi)
    else out <- list(n = n, cluster.number = cn, min.cluster.size = min(cluster.size[1:cwn]),
                     noisen = noisen, average.between = average.between,
                     average.within = average.within, max.diameter = max(diameter[1:cwn],
                                                                         na.rm = TRUE), min.separation = sepwithnoise *
                       min(separation) + (!sepwithnoise) * min(separation[1:cwn]),
                     ave.within.cluster.ss = within.cluster.ss/(n - noisen),
                     avg.silwidth = avg.width, g2 = g2, g3 = g3, pearsongamma = pearsongamma,
                     dunn = dunn, dunn2 = dunn2, entropy = h1, wb.ratio = average.within/average.between,
                     ch = ch, widestgap = widestgap, sindex = sindex,
                     corrected.rand = corrected.rand, vi = vi)
  }
  out
}

##############################
##############################
###

NSCA <- function(N, v = 80){
  I <- nrow(N) # Number of rows of table
  J <- ncol(N) # Number of columns of table
  Inames <- dimnames(N)[1] # Row category names
  Jnames <- dimnames(N)[2] # Column category names
  n <- sum(N) # Total number of classifications in the table
  p <- N * (1/n) # Matrix of joint relative proportions
  Imass <- as.matrix(apply(p, 1, sum))
  Jmass <- as.matrix(apply(p, 2, sum))
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
  dimnames(fbip) <- list(paste(Inames[[1]]), paste(1:min(I,J)))


  f <- sva$u %*% dmu # Row Principal Coordinates
  g <- dJh %*% sva$v %*% dmu # Column Principal Coordinates
  dimnames(f) <- list(paste(Inames[[1]]), paste(1:min(I,J)))
  dimnames(g) <- list(paste(Jnames[[1]]), paste(1:min(I,J)))


  Principal.Inertia <- diag(t(f[,1:min(I-1,J-1)]) %*% f[,1:min(I-1,J-1)])
  Total.Inertia <- sum(Principal.Inertia)
  tau.num <- Total.Inertia # Numerator of Goodman-Kruskal tau
  tau.denom <- 1 - sum(Imass^2) # Denominator of Goodman-Kruskal tau
  tau <- tau.num/tau.denom # Goodman-Kruskal tau index
  Ctau <- (n-1) * (I-1) * tau # Light & Margolin's C-statistic
  Percentage.Inertia <- 100 * (Principal.Inertia/tau.num)
  Cumm.Inertia <- cumsum(Percentage.Inertia)
  Inertia <- cbind(Principal.Inertia, Percentage.Inertia, Cumm.Inertia)
  dimnames(Inertia)[1] <- list(paste("Axis", 1:min(I - 1, J - 1), sep = " "))
  q.value <- 1 - pchisq(Ctau, df = (I - 1) * (J - 1))
  inner.prod <- fbip[,1:which(Cumm.Inertia >= v)[1]] %*% t(g[,1:which(Cumm.Inertia >= v)[1]])
  rownames(inner.prod) <- rownames(fbip)
  tau.num.j <- apply(apply(g^2,2,function(x) Jmass*(x)/tau.num),1,sum)
  names(tau.num.j) <- rownames(Jmass)

  list(N = N, f = f, fbip = fbip, g = g, tau = tau, tau.num.j = tau.num.j,
       di = apply(f,1, function(x) sum(x^2)), dj = apply(g,1, function(x) sum(x^2)),
       Cstat = Ctau, Total.Inertia = Total.Inertia, P.Value = q.value,
       Inertia = Inertia, inner.prod = inner.prod)
}

######################################################

range01 <- function(x) {(x-min(x))/(max(x)-min(x))}
