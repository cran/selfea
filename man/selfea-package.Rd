\name{selfea-package}
\alias{selfea-package}
\alias{selfea}
\docType{package}
\title{
Selfea: A R package for reliable feature selection
}
\description{
description
}
\details{
\tabular{ll}{
Package: \tab selfea\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-03-31\cr
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
\keyword{ Quasi-Poisson, Negative Binomial, Cohen's f, Cohen's w, ANOVA }
\seealso{

\code{\link{get_statistics_from_file}}
\code{\link{top_table}}
\code{\link{get_statistics_from_dataFrame}}

}
\examples{
library(selfea)

## Load Gregori's data and test Selfea

## Josep Gregori, Laura Villareal, Alex Sanchez, Jose Baselga, Josep Villanueva (2013). 
## An Effect Size Filter Improves the Reproducibility 
## in Spectral Counting-based Comparative Proteomics. 
## Journal of Proteomics, DOI http://dx.doi.org/10.1016/j.jprot.2013.05.030')

## Description:
## Each sample consists in 500ng of standard yeast lisate spiked with 
## 100, 200, 400 and 600fm of a mix of 48 equimolar human proteins (UPS1, Sigma-Aldrich).
## The dataset contains a different number of technical replimessagees of each sample

## import Gregori data
data(example_data1)
df_contrast <- example_data
df_group <- example_group

## calculate statistics including Cohen's effect sizes and p-values
## To see detail of method option, check \code{\link{get_statistics_from_dataFrame}}.
list_result <- get_statistics_from_dataFrame(df_contrast,df_group,padj = 'fdr')

## get significant features by desired statistical power and alpha
## For this example, we set p-value threshold = 0.05, power = 0.84
## To see detail of method option, check \code{\link{top_table}}.
significant_qpf <- top_table(list_result,pvalue=0.05,power_desired=0.84,method='QPF')
}
