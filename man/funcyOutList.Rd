%
%  Copyright (C) 2011-2015 Christina Yassouridis
% 
%

\name{funcyOutList-class}
\Rdversion{0.8}
\docType{class}
\alias{funcyOutList-class}
\alias{funcyOutList}

\title{Class \code{"funcyOutList"}}
\description{
  Return object, created by calls of the function \code{\link{funcit}}.
}

\section{Slots}{
  \describe{
    \item{\code{call}:}{Method call of \code{\link{funcit}}.}
    \item{\code{models}:}{List of all \code{\link{funcyOut}}-objects.}
    \item{\code{data}:}{Input data.}
    \item{\code{timeNr}:}{Number of time points.}
    \item{\code{reg}:}{Regular or irregular data (see \code{\link{formatFuncy}}). }
    \item{\code{k}:}{Number of clusters.}
    \item{\code{methodName}:}{Method names.}
    \item{\code{allClusters}:}{Matrix of all cluster results. Result for
      each method in column.}
    \item{\code{randIndex}:}{Matrix of Rand indices showing the
      similarity between the methods. If true cluster
      membership was given, correct classification on diagonal.}
    \item{\code{votedCluster}:}{Cluster, which majority of methods voted for.}
    \item{\code{accordance}:}{Percentage of methods voting for the voted
    Cluster.}
  }
}

\author{
Christina Yassouridis. 
}

\seealso{
  \code{\link{funcyOut}}
}

\examples{
set.seed(2808)
ds <- sampleFuncy(obsNr=30, k=4, timeNr=7)
data <- Data(ds)
clusters <- Cluster(ds)

res <- funcit(data=data, clusters=clusters, seed=2808,
              methods=1:5, k=4, parallel=TRUE)
class(res)
summary(res)
}

\keyword{classes}
