%
%  Copyright (C) 2011-2015 Christina Yassouridis
% 
%

\name{dist2centers}
\alias{dist2centers}

\title{
Calculates curve distances to the centers.
}
\description{
Distance to all cluster centers is calculated for each curve. Dataset
can be in general format "Format1" or matrix format "Format2" (see \code{\link{formatFuncy}}).
}
\usage{
dist2centers(data, centers)
}

\arguments{
  \item{data}{
    Dataset in \code{"Format1"} or \code{"Format2"}.
  }
  \item{centers}{
    Matrix with cluster center curves in
    columns. 
  }
}
\details{
Calculates the sum of the squared distances between the curves and the cluster
centers on the available time points and divides it by the total number
of time points for the corresponding curve. If the dataset is in irregular format \code{"Format1"}, center curves must
correspond to evaluations on the unique union of all time points.
}
\value{
Matrix of dimension nr_curves x nr_clusters.
}

\author{
Christina Yassouridis
}


\examples{
##Generate dataset
set.seed(2804)
ds <- sampleFuncy(obsNr=100, k=4, timeNrMin=5, timeNrMax=10, reg=FALSE)
data <- Data(ds)
##Get unique union of all time points
time <- sort(unique(data[,3]))
##Generate center curves
c1 <- sin(time)+rnorm(30)
c2 <- cos(time)+rnorm(30)
c3 <- time^2+rnorm(3)
centers <- cbind(c1,c2,c3)
##Calculate distance to the centers
dist2centers(data, centers)
}
