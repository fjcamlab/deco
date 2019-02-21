# deco
R package
---------------------------------------------------
DECO:
Decomposing heterogeneous population cohorts for patient stratification and discovery of biomarkers using omic data profiling.
---------------------------------------------------
We present DECO, a new bioinformatic method to explore and find differences in heterogeneous large datasets usually produced in biological or biomedical omic-wide studies. This method makes a comprehensive analysis of multidimensional datasets, consisting on a collection of samples where hundreds or thousands of features have been measured with a large-scale high-throughput technology, for example, a genomic or proteomic technique. The method finds the differences in the profiles of the features along the samples and identifies the associations between them, showing the features that best mark a given class or category as well as possible sample outliers that do not follow the same pattern of the majority of the corresponding cohort. Interestingly, it can be used for the comparison of two or more classes of samples or for unsupervised comparisons. **DECO** allows the discovery of multiple classes or categories and is quite adequate for patients stratification. 

The statistical procedure followed in both parts of the method are detailed in the original [publication](#references) [1]. A detailed **vignette** is included to explain how to use **DECO** for the analysis of multidimensional datasets, which may include heterogeneous samples or categories. The aim is to improve characterization and stratification of complex sample series, mostly focusing on large patient cohorts, where the existence of outlier or mislabeled samples is quite possible.

**DECO** performs a recursive exploration of differential signal changes between samples, finding variables assigned to: 
**(i)** the main classes or groups of samples that are in the studied cohorts 
**(ii)** significant variation or alteration among certain individuals (related or not to an a-priori known class) 
**(iii)** outlier patterns within feature profiles 
**(iv)** sample outliers (i.e. individuals that behave in a different way to the main groups and have specific markers). 

![Workflow](https://github.com/fjcamlab/deco/blob/master/vignettes/figures/flow.jpg)

---------------------------------------------------
**INSTALLATION**

The **deco** R source package can be directly downloaded from **[Bioconductor repository](https://bioconductor.org/packages/deco/)** or **[GitHub repository](https://github.com/fjcamlab/deco)**. This R package contains a experimental dataset as example, two pre-run R objects and all functions needed to run a DECO analysis.

```r
## Bioconductor repository
if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
  }
BiocManager::install("deco")

## GitHub repository using devtools
BiocManager::install("devtools")
devtools::install_github("fjcamlab/deco")
```

---------------------------------------------------
**EXAMPLE OF DIFFERENTIAL ANALYSIS USING DECO**

**DECO** R package has been designed to perform a deep analysis of heterogeneous populations through two main steps: (i) RDA, including a subsampling procedure (`decoRDA()`) based on LIMMA and (ii) a posterior integration with NSCA (`decoNSCA()`) to find out new subclasses of samples. As a consequence of these two steps, **DECO** will calculate a new *h-statistic* which replaces the original omic data to improve the sample stratification. The *h-statistic* integrates both omic dispersion and predictor-response information given by NSCA.

Pipeline using two categories of samples to compare:
```r
### Loading the R packages
library(deco)
library(BiocParallel) # for parallel computation

# Computing in shared memory
bpparam <- MulticoreParam()

########################
# Loading example data #
########################
## Data from two subtypes (ALK+ and ALK-) of Anaplastic Large Cell Leukemia (ALCL).
data(ALCLdata)

# It includes a subset from an ALCL transcriptomic dataset (GEO id: GSE65823).
# to see the example SummarizedExperiment object
ALCL

## Classes vector to run a binary analysis to compare both classes.
classes.ALCL <- colData(ALCL)[,"Alk.positivity"]
names(classes.ALCL) <- colnames(ALCL)

#######################################################################
# RUNNING SUBSAMPLING OF DATA: BINARY design (two classes of samples) #
#######################################################################
# if annotation and rm.xy == TRUE, then
library(Homo.sapiens)

sub.ma.3r.1K <- decoRDA(data = assay(ALCL), classes = classes.ALCL, q.val = 0.01,
                rm.xy = TRUE, r = NULL, control = "pos", annot = FALSE, bpparam = bpparam,
                id.type = "ENSEMBL", iterations = 10000, pack.db = "Homo.sapiens")

#########################################################################################
# RUNNING NSCA STEP: Looking for subclasses within a category/class of samples compared #
#########################################################################################

deco.results.ma <- decoNSCA(sub = sub.ma.3r.1K, v = 80, method = "ward.D", bpparam = bpparam,
                       k.control = 3, k.case = 3, samp.perc = 0.05, rep.thr = 10)

# Phenotypical data from TCGA RNAseq samples.
colData(ALCL)

# h-statistic matrix used to stratify samples
hMatrix <- NSCAcluster(deco.results.ma)$Control$NSCA$h   #for control samples
hMatrix <- NSCAcluster(deco.results.ma)$Case$NSCA$h      #for case samples

```

Finally, the third main function `decoReport()` will generate a PDF file containing a summary of the analysis including: top-discriminant features, new subclasses of samples found, and several plots showing any relevant result of the analysis.

```r
########################################################
# PDF report with feature-sample patterns or subgroups #
########################################################
## Generate PDF report with relevant information and several plots.

## Binary example (ALK+ vs ALK-)
decoReport(deco.results.ma, sub.ma.3r.1K,
          pdf.file = "report_example_microarray_binary.pdf",
          info.sample = as.data.frame(colData(ALCL)[,8:10]),
          cex.names = 0.3, print.annot = TRUE)
```

---------------------------------------------------
**REFERENCES**

1: Campos-Laborie FJ, Risueño A, Roson-Burgo B, Droste C, Fontanillo C, Ortiz-Estevez M, Trotter MW, Sánchez-Santos JM and De Las Rivas J (2019). **DECO: decompose heterogeneous population cohorts for patient stratification and discovery of sample biomarkers using omic data profiling.** Article in revision.

2: Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W and Smyth GK (2015). **limma powers differential expression analyses for RNA-sequencing and microarray studies.** *Nucleic Acids Res.*, 43:e47. \doi{doi:10.1093/nar/gkv007}.



---------------------------------------------------
**NEWS**

Changes in version 0.99.46 (21-02-2019)
+ Bug fixes.
+ Annotation adapted to new "OrganismDbi" related packages. 
+ Three new diagnostic (plot) functions.
+ Enlarged vignette.

Changes in version 0.99.42 (17-12-2018)
+ Bug fixes.
+ Vignette converted into HTML format.
+ Accepted in Bioconductor. 

Changes in version 0.99.0 (06-11-2018)
+ Submitted to Bioconductor.

