%
%  Copyright (C) 2011-2015 Christina Yassouridis
% 
%

\name{funcyOut-class}
\Rdversion{0.8}
\docType{class}
\alias{funcyOut}
\alias{funcyOut-class}
\title{Class \code{"funcyOut"}}
\description{
List entry of \code{\link{funcyOutList}} for an object created by calls to the function \code{\link{funcit}}.

}
\section{Slots}{
  \describe{
    \item{\code{control}:}{Object of class \code{\link{funcyCtrl}}.}
    \item{\code{methodName}:}{Name of the method. }
    \item{\code{kOut}:}{Number of output clusters (can be smaller than
      k, if method did not find k clusters). }
    \item{\code{dimBaseOut}:}{Output dimensions of the basis functions
      (only relevant if flexDim in \code{\link{funcyCtrl}} is \code{TRUE}).}
    \item{\code{time}:}{Vector of measurement time points.}
    \item{\code{cluster}:}{Vector of cluster outcomes.}
    \item{\code{centers}:}{Matrix of center functions with centers
      stored in columns.}
    \item{\code{props}:}{Vector of cluster proportions.}
    \item{\code{dist2centers}:}{Distance to the centers. Matrix of
      dimension nr_curves x nr_clusters. }
    \item{\code{cldist}:}{Matrix of dimension nr_curves x 2. Distances
      to closest and second closest cluster center.}
    \item{\code{calcTime}:}{Calculation time, object of class \code{"proc_time"}. }
    \item{\code{plotParams}:}{Plot parameters, only relevant for call to
    \code{\link{plot}}.}
    \item{\code{correctCl}:}{Numeric, Rand index if correct clusters
      were given as input.}
  }
  Objects of class \code{funcyOut-iterSubspace} inherit from
  \code{funcyOut} and have the following additional slots:
  \describe{
    \item{\code{groupDimBase}:}{Cluster specific dimensions.}
	\item{\code{prms}:}{List of model specific parameters.}
	\item{\code{nrIter}:}{Number of iterations.}
      }
Object of class \code{funcyOut-Mbc} inherit from
\code{funcyOut-iterSubspace} and can have the following additional
slots:
\describe{
  \item{\code{probs}:}{Class probabilities for each curve, matrix of
    dimension nr_curves x nr_clusters.}
  \item{\code{AIC}:}{AIC.}
  \item{\code{BIC}:}{BIC.}
  \item{\code{loglik}:}{Log-likelihood.}
}
Objects of class \code{funcyOutMbc-fclust} inherit from
\code{funcyOut-Mbc} and can have the following additional slots:
\describe{
\item{\code{fit}:}{Output needed for the method specific plots.}
  }

Objects of class \code{funcyOutMbc-fscm} inherit from
\code{funcyOut-Mbc} and can have the following additional slots:
\describe{
  \item{\code{trends}:}{Cluster trends.}
  \item{\code{location}:}{Location matrix of the curves.}
}
    }
    

\author{
Christina Yassouridis
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{funcyOutList}}
}


\keyword{classes}
