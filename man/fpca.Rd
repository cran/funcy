%
%  Copyright (C) 2011-2015 Christina Yassouridis
% 
%

\name{fpca}
\alias{fpca}

\title{
  Functional principal component analysis
}
\description{
  Performs a functional principal component analysis.
}
\usage{
fpca(data=data, dimBase=4, fpcCtrl=NULL)
}

\arguments{
  \item{data}{ Dataset in \code{"Format1"} for irregular or
  \code{"Format2"} for regular time measurements (see \code{\link{formatFuncy}}).}

\item{dimBase}{ Dimension of the basis.}

\item{fpcCtrl}{ Control object, see \code{\link{fpcCtrl}}.  } }
  \details{ A functional principal component analysis is
  executed. Therefore, a smoothed mean curve and a smoothed covariance
  matrix are built and upon these, eigenvalues and eigenvectors are
  calculated. If \code{average}=TRUE, data is stored in a sparse matrix
  of dimension (nr_curves x nr_union_time_points) and matrix operations
  speed up the calculation. However, if the curves do not have many
  common time points, average should be set to FALSE. If the smoothing
  parameter (either for the mean or for the covariance) specified in
  \code{\link{fpcCtrl}} led to NA-values in the smoothing process, it is
  automatically increased by 10\% of its actual value and a warning is
  printed. This process is repeated until smoothing is possible. }

\value{
 \item{yreg}{Smoothed dataset evaluated at union time points.}
 \item{time}{Vector of union time points.}
 \item{meanfcn}{Smoothed mean function.}
 \item{covfcn}{Smoothed covariance matrix.}
 \item{base}{Basis functions.}
 \item{eigval}{Eigenvalues.}
 \item{coeffs}{Coefficients for the basis functions.}
 \item{varprop}{Proportion of the variance.}

}

\author{
  Christina Yassouridis
}


\seealso{
  \code{\link{fpcCtrl}}
  }
  \examples{
##Generate a regular dataset
set.seed(2804)
ds <- sampleFuncy(obsNr=40, timeNr=15, reg=TRUE)

##Execute functional principal component analysis
res <- fpca(Data(ds))
matplot(res$base,type='l')

##Generate an irregular dataset
set.seed(2804)
ds <- sampleFuncy(obsNr=30, k=5, timeNrMin=5, timeNrMax=7, reg=FALSE)

##Execute functional principal component analysis
res <- fpca(Data(ds))
matplot(res$base, type='l', main="First 4 basis functions.")
}

