%
%  Copyright (C) 2011-2015 Christina Yassouridis
% 
%

\name{funcyCtrl}
\docType{class}
\alias{funcyCtrl}
\alias{funcyCtrl-class}
\alias{coerce,list,funcyCtrl-method}
\title{Class \code{"funcyCtrl"}}
\description{
  Hyperparameters for functional cluster algorithms.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("funcyCtrl")}. In addition, named lists can be
  coerced to \code{funcyCtrl} objects, names are completed if unique (see examples).
}
\section{Slots}{
  Objects of class \code{funcyCtrl} have the following slots:
  \describe{
    \item{\code{baseType}:}{Type of basis functions, one of \code{"eigenbasis"}, \code{"splines"},
      \code{"exponential"},
      \code{"fourier"},
      \code{"power"},
      \code{"polynomial"}
      .}
    \item{\code{dimBase}:}{Dimension of the basis functions.}
    \item{\code{flexDim}:}{If \code{TRUE}, dimension can vary between
      clusters (if supported by the algorithm). \code{dimBase} is therefore the
      maximum dimension.}
    \item{\code{init}:}{Algorithm for initial clustering, one of
      \code{"kmeans"}, \code{"random"} or \code{"hclust"}.}
    \item{\code{nrep}:}{Number of replications for initial clustering.}
    \item{\code{seed}:}{Seed number.}
    \item{\code{thd}:}{Threshold if fpca was integrated into method.}
    \item{\code{redDim}:}{Reduced dimension if coefficients are additionally projected onto lower subspace.}
  }
  
  Objects of class \code{funcyCtrlMbc} inherit from
  \code{funcyCtrl} and have the following additional slots:
  \describe{
    \item{\code{eps}:}{Convergence threshold for EM-algorithm.} 
    \item{\code{maxit}:}{Maximum number of iterations.}
    \item{\code{hard}:}{Hard classification?}
  }
}

\author{
  Christina Yassouridis
}

\seealso{
  \code{\link{funcit}}
}

\examples{
##Show slots
showClass("funcyCtrl")

##Define new parameters
mycont = new("funcyCtrl", baseType="fourier", dimBase=4, flexDim=TRUE,
init="hclust")
mycont
}

\keyword{classes}
