\name{example_data2}
\alias{example_data2}
\docType{data}
\title{
Yeast lisate samples spiked with human proteins
}
\description{
The spectral counts matrix has samples in the columns, and proteins in the rows.
Each sample consists in 500ng of standard yeast lisate spiked with 100, 200, 400 and 600fm of a mix of 48 equimolar human proteins (UPS1, Sigma-Aldrich). The dataset contains a different number of technical replicates of each sample.
}
\usage{
data(example_data2)
}
\format{
  Two data frames, df_contrast (protein expression profile) and df_group (experiment group information).
}
\references{
Josep Gregori, Laura Villareal, Alex Sanchez, Jose Baselga, Josep Villanueva (2013). An Effect Size Filter Improves the Reproducibility in Spectral Counting-based Comparative Proteomics. Journal of Proteomics, 
DOI \url{http://dx.doi.org/10.1016/j.jprot.2013.05.030}

Laurent Gatto and Kathryn S. Lilley, MSnbase - an R/Bioconductor package for isobaric tagged mass spectrometry data visualization, processing and quantitation, Bioinformatics 28(2), 288-289 (2012).
}

\examples{
data(example_data1)
}
\keyword{datasets}
