########################################################################################
### DECO ###
# algorithm to find differential significant markers in heterogeneous cohorts
# using jackknife-like and Non-Symmetrical Correspondence Analysis approach
#
# Authors: Francisco J. Campos-Laborie, Jose Manuel Sanchez-Santos & Javier De Las Rivas
# Bioinformatics & Functional Genomics Group
# Cancer Research Center (CiC-IBMCC, CSIC/USAL)
# Salamanca, Spain
# http://www.cicancer.org/
# http://bioinfow.dep.usal.es
#
# version beta: v1.0
# All Rights Reserved. License to be established.
# Last update: 03-11-2016
########################################################################################

# Run in R
# R version 3.2.3 (2015-12-10)
# Platform: x86_64-redhat-linux-gnu (64-bit)
# require(Biobase)
# require(AnnotationDbi)
# require(limma)
# require(snowfall)
# require(foreign)
# require(plsdof)
# require(cluster)
# require(RColorBrewer)
# require(gplots)
# require(made4)
# require(ade4)
# require(gdata)
# require(locfit)
# require(org.Hs.eg.db)
# require(sfsmisc)
# require(scatterplot3d)
#

############################################
slots <- c("data.frame","list","data.frame","factor","numeric","character","numeric","numeric","numeric","call","call")
names(slots) <- c("featureTable","NSCAcluster","incidenceMatrix","classes","pos.iter","control",
                  "q.val", "rep.thr", "samp.perc","subsampling.call","deco.call")
setClass(Class = "deco", slots = slots, sealed = FALSE)


############################################
options(repos = c(CRAN="http://cran.r-project.org"))


###################################
###         decoRDA           ###
###################################
## Recursive function to select differentially expressed features using LIMMA R package (eBayes method) along the data.

decoRDA <- function(data, classes = NA, control = NA, r = NULL,
                    q.val = 0.01, iterations = 10000, cpus = 2, parallel = FALSE,
                    temp.path = tempdir(), annot = FALSE, id.type = NA,
                    attributes = NA, rm.xy = FALSE, pack.db = "org.Hs.eg.db")
{
  call <- match.call()
  options(warn=-1)

  # Setting up temporary directory to write intermediate results.
  if(temp.path == tempdir() || !(temp.path %in% names(dir()))){
    temp.path <- tempdir(); message("'temp' folder will be created in current work directory to store internal loop data.")}

  # Assessing 'data' input.
  if(is(data, "ExpressionSet"))
    data <- exprs(data)
  if(!is.matrix(data))
    stop("Data must be a matrix or ExpressionSet")
  data <- as.matrix(data)

  # Annotating IDs to chromosome location (locus) in order to remove those placed on X or Y.
  # It is proposed to avoid significant results from unbalanced gender contrasts.
  if(rm.xy & annot){
    if(!exists(pack.db)) stop(paste("ERROR: Annotation package '",pack.db,"' has not been loaded."))
    else requireNamespace(pack.db, quietly = TRUE)
    if(id.type %in% columns(x = get(pack.db)))
    {
      infogenes <- AnnotateDECO(ids = rownames(data),
                                id.type = id.type, attributes = c("CHR","CHR"),
                                pack.db = pack.db, verbose = FALSE)
      if(length(which(infogenes[,"CHR"] %in% c("X","Y"))) > 1){
        data <- data[rownames(infogenes)[!infogenes[,"CHR"] %in% c("X","Y")],]
        message(paste(format(Sys.time(), " %H:%M:%S"),"-- Features located in X or Y chromosome have been filtered:\n",
                      length(which(infogenes[,"CHR"] %in% c("X","y"))),"features have been discarded."))
      }else
        message(paste(format(Sys.time(), " %H:%M:%S")," All features are placed on X or Y chromosome, filter will not be applied."))
    }else
      message(paste(format(Sys.time(), " %H:%M:%S")," It was not possible to assign chromosome location to input IDs."))
  }
  # Defining classes and sample sizes for 'supervised' design.
  multi <- FALSE
  if(all(!(is.na(classes))))
  {
    if(!(all(names(classes) %in% colnames(data))) | is.null(names(classes)))
      stop("ERROR: Names of samples in classes vector not match with data.")
    classes <- sort(as.factor(classes))
    # If multiclass vector has been proposed:
    if(length(levels(classes)) > 2){
      message(paste(format(Sys.time(), " %H:%M:%S")," 'Control' label will be NA."))
      control <- NA
      combMulti <- t(combn(levels(classes),2))
      # Removing auto-combination
      combMulti <- combMulti[!apply(combMulti, 1, function(x) x[1] == x[2]),]
      # Calculating maximal resampling size per class
      if(is.null(r) || r <= 0)
        r = round(sqrt(min(table(classes))),digits = 0)
      maxSub <- subsamplingProb(x = r, n1 = min(table(classes)), iter = 0)$n.all
      if(iterations == 0)
        iterations <- maxSub
      if(is.na(maxSub))
        stop("ERROR: Resampling size is higher than minimum class size.")
      multIter <- round(iterations/dim(combMulti)[1], 0)
      if(multIter > maxSub)
        multIter <- maxSub
      message(" NOTE: More than two classes of samples. Multiclass analysis must run ",dim(combMulti)[1]," (or multiple) iterations at least:\n ",
              multIter*dim(combMulti)[1]," random iterations (",multIter," rounds) will be calculated.\n")
      c1 <- apply(combMulti, 1, function(y) names(classes[classes == y[1]]))
      c2 <- apply(combMulti, 1, function(y) names(classes[classes == y[2]]))
      combi <- matrix(data = 0, ncol = 2*r, nrow = 1)
      # Calculating combinations of samples without mixing classes.
      for(i in seq_len(multIter)){
        if(is.list(c1)){
          combi1 <- matrix(unlist(lapply(c1, function(x) which(colnames(data) %in% sample(x, size = r, replace = FALSE)))), ncol = r, byrow = TRUE)
          combi2 <- matrix(unlist(lapply(c2, function(x) which(colnames(data) %in% sample(x, size = r, replace = FALSE)))), ncol = r, byrow = TRUE)
        }else{
          combi1 <- t(apply(c1, 2, function(x) which(colnames(data) %in% sample(x, size = r, replace = FALSE))))
          combi2 <- t(apply(c2, 2, function(x) which(colnames(data) %in% sample(x, size = r, replace = FALSE))))
        }
        combi <- rbind(combi, cbind(combi1, combi2))
      }
      results <- cbind(combi[2:dim(combi)[1],], nFeatures=c(rep(0,multIter*dim(combMulti)[1])))
      colnames(results) <- c(seq_len(2*r), "nFeatures")
      n1 <- dim(data)[2]
      n2 <- dim(data)[2]
      categories.control <- seq_len(dim(data)[2])
      categories.case <- seq_len(dim(data)[2])
      multi <- TRUE
    }
    # Binary contrast. 'Control' label defines '0' class for eBayes algorithm.
    else{
      if(is.na(control))
      {
        cl1 <- levels(classes)[1]
        cl2 <- levels(classes)[2]
        control <- cl1
      }else
      {
        if(!(control %in% levels(classes))){stop("Control label not found in vector classes provided.")}
        cl1 <- control
        cl2 <- levels(classes)[which(levels(classes) != control)]
      }
      n1 <- length(which(classes == cl1))
      n2 <- length(which(classes == cl2))
      classes <- classes[c(which(classes == cl1),which(classes == cl2))]
      data <- data[,names(classes)]
      categories.control <- seq_len(n1)
      categories.case <- (n1+1):(n1+n2)
      message("\n Classes vector defined as input. SUPERVISED analysis will be carry out.\n")}
  }else
  {
    n1 <- dim(data)[2]
    n2 <- dim(data)[2]
    categories.control <- seq_len(dim(data)[2])
    categories.case <- seq_len(dim(data)[2])
    message("\n Classes vector not defined as input. UNSUPERVISED analysis will be carry out.\n")
  }
  # Setting up subsampling size 'r'.
  if(is.null(r) || r <= 0)
    r = round(sqrt(min(c(n1,n2))),digits = 0)
  if(r < 3){r = 3; message(" Resampling size set to 3.\n")}
  if(r > n1 || r > n2){if(parallel){sfStop()};stop("Resampling size can not be higher than samples.")}
  # Message to user
  if(all(is.na(classes)))
    message(paste(format(Sys.time(), " %H:%M:%S"),"-- Resampling design:\n Unsupervised analysis for",
                  dim(data)[2],"samples.\n Resampling size:",r,"\n adj.p.value threshold:",q.val))
  else if(!multi)
    message(paste(format(Sys.time(), " %H:%M:%S"),"-- Resampling design:\n Binary analysis for",
                  dim(data)[2],"samples.\n Control or 0 --> '",cl1,"' with",n1,"samples.\n Case or 1    --> '",cl2,"' with",
                  n2,"samples.\n","Resampling size for both classes:",r,"\n adj.p.value threshold:",q.val))
  else
    message(paste(format(Sys.time(), " %H:%M:%S"),"-- Resampling design:\n Multiple supervised analysis for",
                  dim(data)[2],"samples and",length(levels(classes)),"classes:\n",paste(levels(classes),collapse = "/"),
                  ".\n Resampling size for all classes:",r,"\n adj.p.value threshold:",q.val))
  cont <- 0
  user <- "n"
  # Asking for number of iterations if it was not defined.
  if(iterations == 0 || is.null(iterations) & !multi)
  {
    n.all <- subsamplingProb(x = r, iter = 1000, n1 = n1, n2 = n2)$n.all
    user <- readline(prompt=message(paste("\n All possible combinations:",n.all,"\n Run ALL combinations? (y/n):")))
    if(user == "n")
    {
      iterations <- as.integer(readline(prompt=message(paste("\n Define number of random combinations to calculate:"))))
      if(!(is.integer(iterations)))
        stop("Non-integer argument passed as number of combinations.")
    }
  }
  # Calculate random or all possible combinations between samples if 'binary' or 'unsupervised' has been proposed.
  if(user == "y" & !multi)
  {
    allCombControl <- combn(n1,r)
    allCombCase <- combn(n2,r)
    nSamples <- nrow(allCombControl)+nrow(allCombCase)
    ncomb <- ncol(allCombControl)*ncol(allCombCase)
    results <- matrix(nrow=ncomb,ncol=nSamples)
    colnames(results) <- seq_len(nSamples)
    cont <- 1
    # Building combination matrix.
    for(iCombCONTROL in seq_len(ncol(allCombControl)))
    {
      for(iCombCASE in seq_len(ncol(allCombCase)))
      {
        samples <- c(categories.control[allCombControl[,iCombCONTROL]],categories.case[allCombCase[,iCombCASE]])
        if(length(which(duplicated(samples))) > 0){samples[duplicated(samples)] <-
          sample(which(!(seq_len(n1+n2) %in% samples)),1)}
        results[cont,seq_len(nSamples)] <- samples
        cont <- cont + 1
      }
    }
    results <- cbind(results,nFeatures=c(rep(0,ncomb)))
    message(paste(format(Sys.time(), "\n %H:%M:%S"),"-- Number of total iterations:",ncomb))
  }else if(!multi)
  {
    results <- matrix(nrow = iterations,ncol = 2*r+1)
    colnames(results) <- c(seq(from=1,to=2*r),"nFeatures")
    # Building combination matrix.
    for(i in seq_len(iterations))
    {
      if(all(!(is.na(classes))))
      {
        results[i, seq_len(r)] <- sample(categories.control,r,replace = FALSE)
        results[i,(r+1):(2*r)] <- sample(categories.case,r,replace = FALSE)
      }
      else
        results[i, seq_len(r*2)] <- sample(categories.control,2*r,replace = FALSE)
    }
    results[,c("nFeatures")] <- 0
    message(paste(format(Sys.time(), "\n %H:%M:%S"),"-- Randomly selected",iterations,"iterations."))
  }else
    message(paste(format(Sys.time(), "\n %H:%M:%S"),"-- Number of total iterations:", dim(results)[1]))

  # Parallel computing and temporary directory.
  if(parallel)
    sfInit(parallel=parallel,cpus=cpus)
  dir.create(file.path(temp.path,"temp"), showWarnings = TRUE)
  # Creating final variable containing all subsampling results.
  res <- list(data=as.matrix(data), results=results, subStatFeature=NULL, incidenceMatrix=NULL,
              classes=classes, resampleSize=r, control=control,
              pos.iter = 0, q.val = q.val, call = call)
  # Some variables definition to parallel subsampling.
  top_eje <- matrix(ncol=7)
  colnames(top_eje) <- c("logFC","AveExpr","t","P.Value","adj.P.Val","B","ID")
  top_eje[,2:7] <- 0
  UP <- as.data.frame(data)
  if(all(is.na(classes)) | multi){UP[,seq_len(n1)] <- 0}else{UP[,seq_len(n1+n2)] <- 0}
  DOWN <- UP
  mx <- UP
  results <- cbind(results,counter = seq(seq_len(dim(results)[1])))

  # Subsampling loop with 'sfApply' or 'apply'.
  if(parallel==TRUE & dim(results)[1] > 1)
  {
    sfExport('results','data','r','q.val','mx','top_eje')
    sfExport(list = list('lmFit','makeContrasts','contrasts.fit',
                         'eBayes','topTable'),
             namespace = "limma")
    sfExport(list = list('write.dta'), namespace = "foreign")
    limma1 <- sfApply(results,1,function(results)
    {
      CASE = CONTROL = NULL
      # Taking different combinations from previous matrix and running LIMMA.
      samples <- results[seq_len(r*2)]
      counter <- results[(r*2)+2]
      design <- cbind(CONTROL=rep(c(1,0),each=r),CASE=rep(c(0,1),each=r))
      fit <- lmFit(data[,samples], design)
      cont.mat <- makeContrasts(CASEvsCONTROL=CASE-CONTROL,levels=design)
      fit2 <- contrasts.fit(fit,cont.mat)
      fit3 <- eBayes(fit2)
      top <- topTable(fit3,adjust.method="fdr",number=nrow(data),p.value=q.val)
      pos <- length(rownames(top))
      # Generating intermediate files in temporary dir. If any DE feature is found for
      # any combination, no file would be created.
      if(pos > 0)
      {
        top <- data.frame(top, ID = rownames(top))
        # colnames(top) <- c("logFC","AveExpr","t","P.Value","adj.P.Val","B","ID")
        top <- top[seq_len(pos),]
        top_eje <- rbind(top_eje,top[,colnames(top_eje)])
        # Differential features and incidence matrix files are separated.
        resultfile1 <- paste(temp.path,"/temp/diff",counter,".dta",sep="")
        resultfile2 <- paste(temp.path,"/temp/incid",counter,".dta",sep="")
        write.dta(top_eje[2:length(rownames(top_eje)),],file=resultfile1)
        write.dta(as.data.frame(samples),file=resultfile2)
      }
      else{resultfile1 <- NA;resultfile2 <- NA}
      return(c(resultfile1,resultfile2))
    })
    limma1 <- limma1[!(is.na(limma1))]

  }else
  {
    limma1 <- apply(results,1,function(results)
    {
      CASE = CONTROL = NULL
      # Taking different combinations from previous matrix and running LIMMA.
      samples <- results[seq_len(r*2)]
      counter <- results[(r*2)+2]
      design <- cbind(CONTROL=rep(c(1,0),each=r),CASE=rep(c(0,1),each=r))
      fit <- lmFit(data[,samples], design)
      cont.mat <- makeContrasts(CASEvsCONTROL=CASE-CONTROL,levels=design)
      fit2 <- contrasts.fit(fit,cont.mat)
      fit3 <- eBayes(fit2)
      top <- topTable(fit3,adjust.method="fdr",number=nrow(data),p.value=q.val)
      pos <- length(rownames(top))

      if(pos > 0)
      {
        top <- data.frame(top, ID = rownames(top))
        # colnames(top) <- c("logFC","AveExpr","t","P.Value","adj.P.Val","B","ID")
        rownames(top) <- as.character(top[,"ID"])
        top <- top[seq_len(pos),]
        top_eje <- rbind(top_eje,top[,colnames(top_eje)])
        # Differential features and incidence matrix files are separated.
        resultfile1 <- paste(temp.path,"/temp/diff",counter,".dta",sep="")
        resultfile2 <- paste(temp.path,"/temp/incid",counter,".dta",sep="")
        write.dta(top_eje[2:length(rownames(top_eje)),],file=resultfile1)
        write.dta(as.data.frame(samples),file=resultfile2)
      }
      else{resultfile1 <- NA;resultfile2 <- NA}
      return(c(resultfile1,resultfile2))
    })
    limma1 <- limma1[!(is.na(limma1))]
  }
  # Stop subsampling for non-differential signal.
  if(length(limma1) == 0){
    if(parallel)
      sfStop()
    stop("WARNING: Any combination shows DE features with these resampling size and p-value threshold.
         Try again with another parameters.")}
  message(paste(format(Sys.time(), "\r %H:%M:%S"),"-- Summarizing results..."))
  mx0 <- mx
  mx20 <- mx
  pb <- txtProgressBar(style = 3, min = 1, max = length(limma1)-1)
  # Reading and summarizing intermediate results.
  for(i in seq(from = 1,to = length(limma1)-1,by = 2))
  {
    mx <- mx0
    mx2 <- mx20
    resultfile1 <- limma1[i];resultfile2 <- limma1[i+1]
    if(is.matrix(limma1)){resultfile1 <- limma1[1,i];resultfile2 <- limma1[2,i]}
    # Joining all differential event statistics.
    top_eje1 <- read.dta(file = resultfile1)
    colnames(top_eje1) <- colnames(top_eje)
    rownames(top_eje1) <- make.names(rep(LETTERS,dim(data)[1])[seq_len(dim(top_eje1)[1])],unique=TRUE)
    samples <- as.vector(t(read.dta(file = resultfile2)))
    counter <- as.numeric(unlist(strsplit(unlist(
      strsplit(resultfile1, split = "diff"))[2], split = ".", fixed = TRUE))[1])
    # results[counter,"nFeatures"] <- dim(top_eje1)[1]
    # Generating global incidence matrix with UP and DOWN events.
    mx[rownames(mx) %in% top_eje1[top_eje1[,c("logFC")]<0,c("ID")],samples[seq_len(r)]] <- mx[
      rownames(mx) %in% top_eje1[top_eje1[,c("logFC")]<0,c("ID")],samples[seq_len(r)]]+1
    mx[rownames(mx) %in% top_eje1[top_eje1[,c("logFC")]>0,c("ID")],samples[(r+1):(r*2)]] <- mx[
      rownames(mx) %in% top_eje1[top_eje1[,c("logFC")]>0,c("ID")],samples[(r+1):(r*2)]]+1
    mx2[rownames(mx2) %in% top_eje1[top_eje1[,c("logFC")]<0,c("ID")],samples[(r+1):(r*2)]] <- mx2[
      rownames(mx2) %in% top_eje1[top_eje1[,c("logFC")]<0,c("ID")],samples[(r+1):(r*2)]]+1
    mx2[rownames(mx2) %in% top_eje1[top_eje1[,c("logFC")]>0,c("ID")],samples[seq_len(r)]] <- mx2[
      rownames(mx2) %in% top_eje1[top_eje1[,c("logFC")]>0,c("ID")],samples[seq_len(r)]]+1

    UP <- UP + mx
    DOWN <- DOWN + mx2
    top_eje <- rbind(top_eje,top_eje1)
    rownames(top_eje1) <- NULL
    setTxtProgressBar(pb, i)
  }
  close(pb)


  # Matrix of combinations and number of events.
  ord <- order(results[,"nFeatures"], decreasing = TRUE)
  res$results <- results[ord,-(r*2+2)]
  res$pos.iter <- length(limma1)/2
  # Making up incidence matrix with separated UP, DOWN and MIXED events.
  both.gene <- names(which(apply(UP,1,function(x) sum(x[seq_len(n1)]) > 0 &&
                                   sum(x[(n1+1):(n1+n2)]) > 0) == TRUE))
  UP <- UP[apply(UP,1,sum) > 0,]
  DOWN <- DOWN[apply(DOWN,1,sum) > 0,]
  MULTI <- UP + DOWN
  if(length(both.gene) > 0)
  {
    both.gene.mx <- UP[apply(UP,1,sum) > 0,][rownames(UP)%in%both.gene,]
    rownames(both.gene.mx) <- paste(rownames(both.gene.mx),sep = "deco","UP")
    UP <- UP[!(rownames(UP)%in%both.gene),]
    DOWN <- DOWN[!(rownames(DOWN)%in%both.gene),]
  }
  # Vector containing direction of differential event if 'binary'
  UpDw <- c(rep("UP",length(rownames(UP)[apply(UP[,seq_len(n1)],1,sum) == 0])),
            rep("DOWN",length(rownames(UP)[apply(UP[,seq_len(n1)],1,sum) > 0])),
            rep("MIXED",length(both.gene)))
  names(UpDw) <- c(rownames(UP)[apply(UP[,seq_len(n1)],1,sum) == 0],
                   rownames(UP)[apply(UP[,seq_len(n1)],1,sum) > 0],both.gene)
  rownames(UP) <- paste(rownames(UP),sep = "deco","UP")
  rownames(MULTI) <- paste(rownames(MULTI),sep = "deco","UP")
  rownames(DOWN) <- paste(rownames(DOWN),sep = "deco","DOWN")
  if(all(!is.na(classes)) & !multi)
    res$incidenceMatrix <- interleave(UP, DOWN, append.source = TRUE,
                                      sep = "deco", drop = FALSE)
  else
    res$incidenceMatrix <- MULTI
  if(length(both.gene) > 0)
    res$incidenceMatrix <- rbind(res$incidenceMatrix, both.gene.mx)
  res$incidenceMatrix <- res$incidenceMatrix[,colnames(res$incidenceMatrix) %in% colnames(data)]

  # Creating table containing all statistical information about features.
  j <- length(limma1)/2
  top_eje <- top_eje[2:length(rownames(top_eje)),]
  res$subStatFeature <- as.data.frame(matrix(ncol=12,nrow=length(unique(top_eje[,c("ID")]))))
  colnames(res$subStatFeature) <- c("ID","UpDw","Avrg.logFC","Best.adj.P.Val",
                                    "Repeats","FR.Repeats","RF.Positive.Repeats","Chi.Square","P.Val.ChiSq",
                                    "ChiSq.adj.P.Val.FDR","ChiSq.adj.P.Val.Holm")
  res$subStatFeature[,c("ID")] <- unique(top_eje[,c("ID")])
  res$subStatFeature[,3:11] <- 0
  ncomb <- dim(results)[1]

  if(parallel)
    sfExport('top_eje','ncomb','res')
  # Filling 'diffgenes' table. It also could be executed in parallel.
  if(dim(res$subStatFeature)[1]>1 && parallel==TRUE)
  {
    limma2 <- sfApply(res$subStatFeature,1,function(name)
    {
      res_ <- vector(length = 12)
      names(res_) <- names(name)
      lines <- top_eje[as.character(top_eje[,c("ID")]) %in% as.character(name[c("ID")]),]
      # Statistics from this RDA step.
      if(length(lines[,c("ID")])>0)
      {
        res_[c("Avrg.logFC")] <- .colMeans(lines[,c("logFC")],
                                           m = length(lines[,c("ID")]), n = 1)
        res_[c("Best.adj.P.Val")] <- lines[1,c("adj.P.Val")]
        res_[c("Repeats")] <- length(lines[,1])
        res_[c("FR.Repeats")] <- res_[c("Repeats")]/ncomb
        res_[c("RF.Positive.Repeats")] <- res_[c("Repeats")]/j
        res_[c("Chi.Square")] <- (-2)*sum(log(lines[,c("adj.P.Val")]))
        res_[c("P.Val.ChiSq")] <-  1-pchisq(as.numeric(res_[c("Chi.Square")]),
                                            df = 2*length(lines[,1]))
        res_[c("ID")] <- name[c("ID")]
      }

      return(res_)
    })
    message(paste(format(Sys.time(), " %H:%M:%S"),"-- Computed",
                  length(unique(top_eje[,c("ID")])),"features.\n"))
    res$subStatFeature <- t.data.frame(as.data.frame(limma2))
    res$subStatFeature[,3:11] <- apply(res$subStatFeature[,3:11],2,as.numeric)
  }else
  {
    for(i in seq_len(length(unique(top_eje[,c("ID")]))))
    {
      lines <- top_eje[top_eje[,c("ID")] %in% as.character(res$subStatFeature[i,c("ID")]),]
      # Statistics from this RDA step.
      if(length(lines[1,c("ID")])>0)
      {
        res$subStatFeature[i,c("Avrg.logFC")] <- .colMeans(lines[,c("logFC")],
                                                           m = length(lines[,c("ID")]), n = 1)
        res$subStatFeature[i,c("Best.adj.P.Val")] <- lines[1,c("adj.P.Val")]
        res$subStatFeature[i,c("Repeats")] <- length(lines[,1])
        res$subStatFeature[i,c("FR.Repeats")] <- res$subStatFeature[i,c("Repeats")]/ncomb
        res$subStatFeature[i,c("RF.Positive.Repeats")] <- res$subStatFeature[i,c("Repeats")]/j
        res$subStatFeature[i,c("Chi.Square")] <- (-2)*sum(log(lines[,c("adj.P.Val")]))
        res$subStatFeature[i,c("P.Val.ChiSq")] <-  1-pchisq(as.numeric(res$subStatFeature[i,c("Chi.Square")]),
                                                            df = 2*length(lines[,1]))
      }
    }
    message(paste(format(Sys.time(), " %H:%M:%S"),"-- Computed",
                  length(unique(top_eje[,c("ID")])),"features.\n"))
  }
  # Adjusting p.values from Chi.Square.
  res$subStatFeature[,c("ChiSq.adj.P.Val.FDR")] <-  p.adjust(
    res$subStatFeature[,c("P.Val.ChiSq")],method="fdr",n=length(res$subStatFeature[,c("P.Val.ChiSq")]))
  res$subStatFeature[,c("ChiSq.adj.P.Val.Holm")] <-  p.adjust(
    res$subStatFeature[,c("P.Val.ChiSq")],method="holm",n=length(res$subStatFeature[,c("P.Val.ChiSq")]))
  res$subStatFeature <- res$subStatFeature[,c("ID","UpDw","Avrg.logFC","Best.adj.P.Val","Repeats",
                                              "FR.Repeats","RF.Positive.Repeats","Chi.Square","P.Val.ChiSq",
                                              "ChiSq.adj.P.Val.FDR","ChiSq.adj.P.Val.Holm")]
  # Annotating features using 'AnnotateDECO' function.
  # Further information in vignette.
  if(annot==TRUE)
  {
    infogenes <- AnnotateDECO(ids = as.character(res$subStatFeature[,c("ID")]),
                              id.type = id.type, attributes = attributes,
                              pack.db = pack.db)
    res$subStatFeature <- cbind(res$subStatFeature,infogenes[as.vector(res$subStatFeature[,c("ID")]),])
  }
  # Making up all statistics.
  res$subStatFeature <- droplevels.data.frame(as.data.frame(res$subStatFeature))
  UpDw[as.character(res$subStatFeature[,c("ID")])]
  res$subStatFeature[,3:11] <- apply(res$subStatFeature[,3:11],2,as.numeric)
  if(all(!is.na(classes)) & length(levels(classes)) == 2)
    res$subStatFeature[,c("UpDw")] <- UpDw[match(as.character(res$subStatFeature[,c("ID")]),
                                                 names(UpDw))]
  else
    res$subStatFeature <- res$subStatFeature[,colnames(res$subStatFeature)!="UpDw"]

  # Standard.Chi.Square
  res$subStatFeature[,"Standard.Chi.Square"] <- ((as.numeric(res$subStatFeature[,c("Chi.Square")])-
                                                    (2*as.numeric(res$subStatFeature[,c("Repeats")])))/
                                                   sqrt(4*as.numeric(res$subStatFeature[,c("Repeats")])))

  message(paste(format(Sys.time(), " %H:%M:%S"),"-- Done.\n"))
  # Removing temporary dir and all intermediate files generated.
  unlink(paste(temp.path,"/temp",sep=""),recursive = TRUE)
  # Parallel execution.
  if(parallel)
    sfStop()

  options(warn=0)
  return(res)
}


#################################
###        decoNSCA           ###
#################################
## Function to carry out factorial analysis with NSCA (Non-Symmetrical Correspondence Analysis) of
## the incidenceMatrix generated by 'DECO' subsampling function.

decoNSCA <- function(sub, rm.complete = FALSE, v = 80, k.control = NULL,
                     k.case = NULL, rep.thr = 3,
                     samp.perc = 0.05,
                     method = "ward.D", parallel = FALSE, cpus = 2,
                     overlap = TRUE)
{
  call <- match.call()
  options(warn=-1)
  # Agglomeration method for hierarchical clustering of samples.
  if(all(!(is.null(method))) && !(method %in% c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")))
    stop("ERROR: Input 'method.heatmap' have to be one of 'hclust'{stats} function methods:
         ward.D, ward.D2, single, complete, average, mcquitty, median or centroid")
  rownames(sub$subStatFeature) <- as.character(sub$subStatFeature[,c("ID")])
  data <- sub$data
  sfInit(cpus = cpus, parallel = parallel)

  # Variable with original IDs used to substitute modified IDs of incidenceMatrix.
  message(paste(format(Sys.time(), "\r %H:%M:%S"),"Applying repeats threshold...\n"))
  g.names <- sfSapply(rownames(sub$incidenceMatrix), function(x) unlist(strsplit(x,split = "deco",fixed = TRUE))[1])
  if(all(is.na(sub$classes)) || length(levels(sub$classes)) > 2){
    f <- sfApply(sub$incidenceMatrix,1,function(x) length(which(x >= rep.thr)))
    names(f) <- g.names[names(g.names)%in%names(f)]
  }
  else{
    f <- sfSapply(unique(g.names), function(x)
      sum(apply(sub$incidenceMatrix[grepl(rownames(sub$incidenceMatrix), pattern = x, fixed = TRUE),],1,function(x)
        length(which(x >= rep.thr)))))
  }
  # 'Repeats' threshold.
  # All features under those thresholds will be removed.
  if(length(which(f > samp.perc*dim(sub$incidenceMatrix)[2])) > 0)
    message(paste(" NOTE: Repeats index threshold have been applied:\n",length(which(f <= samp.perc*dim(sub$incidenceMatrix)[2])),
                  "features removed from",dim(sub$subStatFeature)[1],"total."))
  sub$subStatFeature <- sub$subStatFeature[names(f[which(f > samp.perc*dim(sub$incidenceMatrix)[2])]),]
  sub$subStatFeature <- cbind(sub$subStatFeature, Repeats.index = f[names(f)%in%rownames(sub$subStatFeature)]/dim(sub$incidenceMatrix)[2]*100)
  sub$subStatFeature <- sub$subStatFeature[order(as.character(sub$subStatFeature[,c("ID")]),decreasing = TRUE),]
  # Removing features from incidenceMatrix with lower 'Repeats' than threshold.
  sub$incidenceMatrix <- sub$incidenceMatrix[rownames(sub$incidenceMatrix) %in%
                                               names(g.names[g.names %in% as.character(sub$subStatFeature$ID)]),]
  if(length(which(f > samp.perc*dim(sub$incidenceMatrix)[2])) < 10)
    stop("After applying 'Repeats' filter, there are not enough features (10 features minimum) to input NSCA.")
  # Setting up both classes and control for 'supervised' design.
  if(all(!(is.na(sub$classes))) & !length(levels(sub$classes)) > 2)
  {
    # Instructions just for 'supervised' RDA design.
    if(is.na(sub$control))
    {
      cl1 <- levels(sub$classes)[1]
      cl2 <- levels(sub$classes)[2]
    }else
    {
      cl1 <- sub$control
      cl2 <- levels(sub$classes)[which(levels(sub$classes) != sub$control)]
      sub$classes <- sub$classes[c(which(sub$classes == cl1),which(sub$classes == cl2))]
    }
    n1 <- length(which(sub$classes == cl1))
    n2 <- length(which(sub$classes == cl2))
    sub$incidenceMatrix <- sub$incidenceMatrix[,names(sub$classes)]
    # Calculating raw 'delta.signal' differences between classes for all differential features.
    delta.signal <- round(apply(
      sub$data[rownames(sub$data)%in%as.character(sub$subStatFeature[,c("ID")]),],1, function(x)
        mean(x[(n1+1):(n1+n2)])-mean(x[seq_len(n1)])),3)
    sub$subStatFeature <- data.frame(sub$subStatFeature[order(sub$subStatFeature[,c("ID")]),],
                                     delta.signal = delta.signal[order(names(delta.signal))])
    # Standard deviation per class
    sd.Ctrl <- apply(sub$data[rownames(sub$subStatFeature),seq_len(n1)],1,sd)
    sd.Case <- apply(sub$data[rownames(sub$subStatFeature),(n1+1):(n1+n2)],1,sd)
    # Classification of feature profiles: 'ideal', 'generic', 'specific, and 'both'.
    if(overlap){
      message(paste(format(Sys.time(), "\r %H:%M:%S"),"-- Calculating overlapping signal per feature..."))
      if(parallel){
        sfExport('sub', 'cl1', 'data')
        sfExport('overlapFeature', namespace = "deco")
        overlap <- sfSapply(rownames(sub$subStatFeature), function(x) {
          overlapFeature(id = x, data = sub$data, classes = sub$classes,
                         control = cl1, analysis = "Binary")
        })
      }else
        overlap <- sapply(rownames(sub$subStatFeature), function(x)
          overlapFeature(id = x, data = sub$data, classes = sub$classes,
                         control = cl1, analysis = "Binary"))
      prof <- sapply(overlap, function(x) max(which(c(0,0.2,0.4,0.75) <= x)))
      profile <- prof
      profile[prof %in% c(3,4)] <- "Minority"
      profile[prof == 2] <- "Majority"
      profile[prof == 1] <- "Complete"
      profile[sub$subStatFeature[names(profile),c("UpDw")] == "MIXED"] <- "Mixed"
      names(profile) <- names(overlap)
      profile <- profile[order(names(profile))]
    }
    else{
      overlap <- rep(NA,dim(sub$subStatFeature)[1])
      profile <- rep(NA,dim(sub$subStatFeature)[1])
      names(overlap) <- rownames(sub$subStatFeature)
      names(profile) <- rownames(sub$subStatFeature)
    }

    # Filter to remove 'Ideal' features.
    if(rm.complete){
      message(paste("NOTE: Features that mark differences among classes
                    (Generic profile) will be discarded for NSCA.\n",
                    length(profile[profile == "Complete"]),"features removed from",
                    dim(sub$subStatFeature)[1],"total."))
      sub$incidenceMatrix <- sub$incidenceMatrix[rownames(
        sub$incidenceMatrix) %in% names(g.names[g.names %in% as.vector(
          sub$subStatFeature[profile != "Complete",c("ID")])]),]
    }
    UpDw <- sub$subStatFeature$UpDw
    names(UpDw) <- as.character(sub$subStatFeature$ID)
  }else
  {
    n1 <- dim(data)[2]
    n2 <- dim(data)[2]
    SD <- apply(sub$data,1,function(x) sd(na.omit(x)))
  }
  # Applying NSCA to all samples or both classes.
  z <- 1
  while(z %in% seq_len(2))
  {
    # NSCA for 'unsupervised' and 'multiclass' design.
    if(all(is.na(sub$classes)) || length(levels(sub$classes)) > 2)
    {
      # Removing features or samples with any differential event.
      mx <- sub$incidenceMatrix[which(apply(sub$incidenceMatrix,1,sum)>0),]
      mx <- mx[,which(apply(sub$incidenceMatrix,2,sum)>0)]
      message(paste(format(Sys.time(), "\r %H:%M:%S"),"-- Calculating all samples together..."))
      # NSCA
      nsc.res <- NSCAcluster(mx = mx, data = sub$data, k = k.control, method.dend = method, id.names = g.names, v = v)
      k <- k.control
      diffMX <- cbind(sub$subStatFeature, nsc.res$infoFeature[rownames(sub$subStatFeature),])
      diffMX <- cbind(diffMX, sd = SD[rownames(sub$subStatFeature)])
      # Order for columns of final table with statistical information.
      coln <- c("ID","SYMBOL","hgnc_symbol","Repeats","Repeats.index","FR.Repeats","Avrg.logFC",
                "Standard.Chi.Square","P.Val.ChiSq","ChiSq.adj.P.Val.FDR",
                "ChiSq.adj.P.Val.Holm","sd","Tau.feature","Dendrogram.group",
                "h.Best","h.Range","GENENAME")
      ord <- na.omit(match(c(coln, colnames(diffMX)[!colnames(diffMX)%in%coln]), colnames(diffMX)))
      z <- 3
    }
    # NSCA for 'control' samples
    if(z==1)
    {
      # Removing features or samples with any differential event.
      mx <- sub$incidenceMatrix[which(apply(sub$incidenceMatrix[,seq_len(n1)],1,sum)>0), seq_len(n1)]
      mx <- mx[,which(apply(sub$incidenceMatrix[, seq_len(n1)],2,sum)>0)]
      message(paste(format(Sys.time(), "\r %H:%M:%S"),"-- Calculating control group..."))
      # NSCA
      nsc.res1 <- NSCAcluster(mx = mx, data = sub$data, k = k.control, method.dend = method, UpDw = UpDw, dir = "UP",
                              id.names = g.names, v = v, label = cl1, raw = apply(sub$data[,colnames(mx)], 1, mean))
      colnames(nsc.res1$infoFeature) <- paste(colnames(nsc.res1$infoFeature),"Ctrl",sep=".")
      diffMX <- cbind(sub$subStatFeature, sd.Ctrl = sd.Ctrl[rownames(sub$subStatFeature)],
                      nsc.res1$infoFeature[rownames(sub$subStatFeature),c("Tau.feature.Ctrl","Dendrogram.group.Ctrl"
                                                                          ,"h.Best.Ctrl","h.Range.Ctrl")])
    }
    # NSCA for 'case' samples
    if(z==2)
    {
      # Removing features or samples with any differential event.
      mx <- sub$incidenceMatrix[which(apply(sub$incidenceMatrix[,(n1+1):(n1+n2)],1,sum)>0),(n1+1):(n1+n2)]
      mx <- mx[,which(apply(sub$incidenceMatrix[,(n1+1):(n1+n2)],2,sum)>0)]
      message(paste(format(Sys.time(), "\r %H:%M:%S"),"-- Calculating case group..."))
      # NSCA
      nsc.res2 <- NSCAcluster(mx = mx, data = sub$data, k = k.case, method.dend = method, UpDw = UpDw, dir = "DOWN",
                              id.names = g.names, v = v, label = cl2, raw = apply(sub$data[,colnames(mx)], 1, mean))
      colnames(nsc.res2$infoFeature)[seq_len(length(colnames(nsc.res2$infoFeature))-1)] <-
        paste(colnames(nsc.res2$infoFeature)[seq_len(length(colnames(nsc.res2$infoFeature))-1)],"Case",sep=".")
      diffMX <- cbind(diffMX, sd.Case = sd.Case[rownames(sub$subStatFeature)], overlap.Ctrl.Case = overlap[rownames(sub$subStatFeature)],
                      nsc.res2$infoFeature[rownames(sub$subStatFeature),c("Tau.feature.Case","Dendrogram.group.Case",
                                                                          "h.Best.Case","h.Range.Case")])
      # Standard.Chi.Square calculation.
      diffMX <- cbind(diffMX, Profile = profile[rownames(sub$subStatFeature)])
      # Order for columns of final table with statistical information.
      coln <- c("ID","SYMBOL","hgnc_symbol","UpDw","Profile","overlap.Ctrl.Case","Repeats","Repeats.index","FR.Repeats",
        "delta.signal","Avrg.logFC","Standard.Chi.Square","P.Val.ChiSq","ChiSq.adj.P.Val.FDR",
        "ChiSq.adj.P.Val.Holm","sd.Ctrl","Tau.feature.Ctrl",
        "Dendrogram.group.Ctrl","h.Best.Ctrl","h.Range.Ctrl","sd.Case",
        "Tau.feature.Case","Dendrogram.group.Case","h.Best.Case",
        "h.Range.Case","GENENAME")
      ord <- na.omit(match(c(coln, colnames(diffMX)[!colnames(diffMX)%in%coln]), colnames(diffMX)))
    }
    z <- z + 1
  }
  sfStop()
  # Making up of final statistical table.
  diffg <- data.frame(diffMX[,c(ord)],diffMX[,na.action(ord)])
  colnames(diffg) <- c(colnames(diffMX)[c(ord)],colnames(diffMX)[na.action(ord)])
  # Creating final result object.
  if(!(all(is.na(sub$classes))) & !length(levels(sub$classes)) > 2)
  {
    ord <- with(diffg, order(Profile, -Standard.Chi.Square))
    diffg <- diffg[ord,]
    results <- new("deco", featureTable = as.data.frame(diffg[,!duplicated(colnames(diffg))]),
                   NSCAcluster = list(Control = nsc.res1, Case = nsc.res2),
                   incidenceMatrix = sub$incidenceMatrix,
                   classes = as.factor(sub$classes), pos.iter = sub$pos.iter,
                   control = as.character(sub$control), q.val = sub$q.val,
                   rep.thr = rep.thr, samp.perc = samp.perc,
                   subsampling.call = sub$call, deco.call = call)
  }else{
    ord <- order(apply(apply(cbind(diffg$Standard.Chi.Square,diffg$Standard.Chi.Square,diffg$Repeats,
                                   diffg$sd,diffg$h.Range),2,function(x)rank(-x)), 1, mean))
    diffg <- diffg[ord,]
    results <- new("deco", featureTable = as.data.frame(diffg[,!duplicated(colnames(diffg))]),
                   NSCAcluster = list(All = nsc.res), incidenceMatrix = sub$incidenceMatrix,
                   classes = as.factor(sub$classes), pos.iter = sub$pos.iter,
                   control = as.character(sub$control), q.val = sub$q.val,
                   rep.thr = rep.thr, samp.perc = samp.perc,
                   subsampling.call = sub$call, deco.call = call)
  }
  options(warn=0)
  return(results)
}
