%
%  Copyright (C) 2011-2015 Christina Yassouridis
% 
%

\name{fpcCtrl-class}
\Rdversion{0.8}
\alias{fpcCtrl-class}
\alias{fpcCtrl}
\alias{fpcCtrlMbc}

\docType{class}
\title{Class "fpcCtrl"}
\description{
  Hyperparameters for functional principal component analysis. See \code{\link{fpca}}.
}
\section{Objects from the class}{
  Objects can be created by calls of the form
  \code{new("fpcCtrl")}. In addition, named lists can be coerced to
  \code{fpcCtrl} objects, names are completed if unique.
}
\section{Slots}{
  Objects of class \code{fpcCtrl} have the following slots:
  \describe{
    \item{\code{select}:}{Character string, one of \code{"automatic"} or
      \code{"manual"} specifying whether the bandwidth for smoothing is
      given or should be calculated based on the data.}
    \item{\code{h1Dim}:}{If \code{select}=\code{"manual"}, bandwidth for
      smoothing the mean.}
    \item{\code{h2Dim}:}{If \code{select}=\code{"manual"}, bandwidth for
      smoothing the covariance.}
    \item{\code{sm1Dim}:}{Character string, name of the mean smoothing function. One of \code{"sm.regression"} or \code{"sm1"}.}
    \item{\code{sm2Dim}:}{Character string, name of the covariance smoothing function. One of \code{"sm.regression"} or \code{"sm2"}.}
    \item{\code{coeffsCalc}:}{Character string, specifying how to
      calculate the coefficients. One of \code{"estimate"}
      or \code{"integrate"}.}
    \item{\code{nrMaxTime}:}{Maximum number of evaluation time points for the
      covariance matrix.}
    \item{\code{average}:}{If TRUE, matrix calculation is used and
      speeds up the calculation if curves have many evaluation time points in common.}
    }
}
\section{Details}{
    \code{"sm.regression"} is a nonparametric regression estimate from
    the R-package \pkg{sm}. \code{"sm1"} and \code{"sm2"} are kernel
    smoothers defined by Chiou2007.
    The coefficients for the basis functions can be computed by
    numerical integration (\code{coeffCalc}=\code{"integrate"}) or by a
    sparse estimation defined by Mueller2005 (\code{coeffCalc}=\code{"estimate"}). 
}

\references{
Chiou J-M and Li P-L. "Functional clustering and identifying substructures
of longitudinal data".  Journal of the Royal Statistical Society: Series
B. 69 (4). pp. 679--699. 2007

F. Yao and H-G. Mueller and J-L. Wang. Functional Data Analysis for
Sparse Longitudinal Data. Journal of the American Statistical
Association. 100 (470). 577--590. 2005
}
\author{Christina Yassouridis}
\keyword{classes}

\examples{
##Have a look at the defaults
new("fpcCtrl")

##Coerce a list
mycont = list(select="automatic", sm1Dim="sm1", sm2Dim="sm2", nrMaxTime=20)
as(mycont, "fpcCtrl")

##Default values for coefficients calculation procedure
new("fpcCtrl")@coeffsCalc 
}
