\name{selfea-package}
\alias{selfea-package}
\alias{selfea}
\docType{package}
\title{
Selfea: R package for reliable feature selection using Cohen's effect sizes
}
\description{
Functions using Cohen's effect sizes (Cohen, Jacob. Statistical power analysis for the behavioral sciences. Academic press, 2013) are provided for reliable feature selection in biology data analysis.  In addition to Cohen's effect sizes, p-values are calculated and adjusted from quasi-Poisson GLM, negative binomial GLM and Normal distribution ANOVA.  Significant features (genes, RNAs or proteins) are selected by adjusted p-value and minimum Cohen's effect sizes, calculated to keep certain level of statistical power of biology data analysis given p-value threshold and sample size.
}
\details{
\tabular{ll}{
Package: \tab selfea\cr
Type: \tab Package\cr
Version: \tab 1.0.1\cr
Date: \tab 2015-06-20\cr
License: \tab GPL-2\cr
}
}
\author{
Lang Ho Lee, Arnold Saxton, Nathan Verberkmoes

Maintainer: Lang Ho Lee <llee27@utk.edu>
}
\references{
Lang Ho Lee, Arnold Saxton, Nathan Verberkmoes, Selfea: A R package for reliable feature selection in process
}
\keyword{ Quasi-Poisson, Negative Binomial, Cohen's f, Cohen's w, Cohen's f2, Cohen's d, ANOVA }
\seealso{

\code{\link{get_statistics_from_dataFrame}} ,
\code{\link{get_statistics_from_file}} ,
\code{\link{top_table}} ,
\code{\link{ttest_cohens_d}}

}
\examples{
library(selfea)

## Test to calculate p-value of Student's t-test and Cohen's d
values <- c(8,10,8,8,11,29,26,22,27,26)
groups <- c("U200","U200","U200","U200","U200","U600","U600","U600","U600","U600")
list_result <- ttest_cohens_d (values, groups, 0.05, 0.90)

## Test selfea for single protein expression
values <- c(6,8,10,29,26,22)
groups <- c("U200","U200","U200","U600","U600","U600")
experiments <- c("exp1","exp2","exp3","exp4","exp5","exp6")

df_expr <- data.frame(ID="Protein_1",exp1=6,exp2=8,exp3=10,exp4=29,exp5=26,exp6=22)
df_group <- data.frame(Col_Name=experiments,Group=groups)
list_result <- get_statistics_from_dataFrame(df_expr,df_group)
top_table(list_result)

## Load Gregori's data and test Selfea

## Josep Gregori, Laura Villareal, Alex Sanchez, Jose Baselga, Josep Villanueva (2013). 
## An Effect Size Filter Improves the Reproducibility 
## in Spectral Counting-based Comparative Proteomics. 
## Journal of Proteomics, DOI http://dx.doi.org/10.1016/j.jprot.2013.05.030')

## Description:
## Each sample consists in 500ng of standard yeast lisate spiked with 
## 100, 200, 400 and 600fm of a mix of 48 equimolar human proteins (UPS1, Sigma-Aldrich).
## The dataset contains a different number of technical replimessagees of each sample

## Import Gregori data
## data(example_data2)  ## if you want to test whole Gregori dataset
data(example_data1)  ## example_data1 has only 50 proteins for fast run

df_contrast <- example_data
df_group <- example_group

## calculate statistics including Cohen's effect sizes and p-values
## To see detail of method option, read R document about get_statistics_from_dataFrame.
list_result <- get_statistics_from_dataFrame(df_contrast,df_group,padj = 'fdr')

## get significant features by desired statistical power and alpha
## For this example, we set p-value threshold = 0.05, power = 0.84
## To see detail of method option, read R document about top_table.
significant_qpf <- top_table(list_result,pvalue=0.05,power_desired=0.84,method='QPF')
}
