\name{concordance-package}
\alias{concordance-package}
\alias{concordance}
\docType{package}
\title{
Calculate and test concordance between datasets
}
\description{
This package calculates the concordance (similarity) between
data sets. Concordance measures are based on the correlation structures
of data sets and can be used to test if a subset of data represents a
larger set or if two data sets can be fused (merged) to increase
power in regression challenges.
}
\details{
\tabular{ll}{
Package: \tab cocordance\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2015-01-02\cr
License: \tab Apache License 2.0\cr
}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
Michael J. Kane
Maintainer: Michael J. Kane <kaneplusplus@gmail.com>
}
%\references{
%~~ Literature or other references for background information ~~
%}
\keyword{data-fusion}
\keyword{data-merging}
\keyword{regression-diagnostics}
\examples{
# Test the concordance of the Sepal.Length and Sepal.Width 
# variables for setosa and virginica in iris.
concordance.test(iris[iris$Species=="setosa",c("Sepal.Length", "Sepal.Width")],
                 iris[iris$Species=="virginica",c("Sepal.Length","Sepal.Width")]
}
