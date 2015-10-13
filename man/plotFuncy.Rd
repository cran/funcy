%
%  Copyright (C) 2011-2015 Christina Yassouridis
% 
%

\name{plotFuncy}
\alias{plotFuncy}
\alias{plotFuncy,matrix-method}
\alias{plotFuncy,sampleFuncy-method}

\title{
Plot functional data. 
}
\description{
Plot functional data, color curves according to cluster labels and add
center curves.
}
\usage{
\S4method{plotFuncy}{matrix}(data, regTime=NULL, col=NULL,
          ctr=NULL, ctrOnly=FALSE, ctrCols=NULL, showLegend=FALSE,
          legendPlace="bottomleft", lty=3, lwd=NULL, \dots)

\S4method{plotFuncy}{sampleFuncy}(data, regTime=NULL, col=NULL,
          ctr=NULL, ctrOnly=FALSE, ctrCols=NULL, showLegend=FALSE,
          legendPlace="bottomleft", lty=3, lwd=NULL, \dots)
}


\arguments{
  \item{data}{
    For irregular time measurements, a matrix consisting of curveID in first
    column, curve evaluations in second and time evaluation points in third
    column (Format1).
    For regular time measurements, a matrix consisting of curves in rows and time evaluations in
    columns (Format2).
    Alternatively an object created by \code{\link{sampleFuncy}}.
  }
  \item{regTime}{
    If Format2, optional vector representing the
    evaluation time points. If regTime=NULL equidistant time points from 1 to
    number of curves are used. 
  }
  \item{col}{
    Numeric vector of cluster labels. 
  }
  \item{ctr}{
     Numeric matrix of cluster centers
   }
   \item{ctrOnly}{
    Plot only centers, no curves?
  }
  \item{ctrCols}{
    Numeric vector specifying the color of the center curves. 
   
  }
  \item{showLegend}{
    Whether to show cluster legend. 
      }
  \item{legendPlace}{
    Where to place cluster legend.
  }
  \item{lty}{
    Line type of the functional data. 
  }
  \item{lwd}{
    Line width of the functional data. 
  }
  \item{\dots}{
    Further plotting parameters. 
  }
}


\author{
Christina Yassouridis
}


\examples{
##cluster the data with methods for regular sets
##sample a regular dataset
set.seed(2001)
ds <- sampleFuncy(obsNr=40, k=4, timeNr=20, reg=TRUE)
data <- Data(ds)
clusters <- Cluster(ds)

##plot sampled functions. 
plotFuncy(data, col=clusters, lty=1, showLegend=TRUE, legendPlace="topleft")

##sample an irregular dataset
set.seed(2001)
ds <- sampleFuncy(obsNr=40, k=4, timeNrMin=2, timeNrMax=10, reg=FALSE)
data <- Data(ds)
clusters <- Cluster(ds)

##plot sampled functions. 
plotFuncy(data, col=clusters, lty=1, showLegend=TRUE, legendPlace="topleft")


##cluster functions with method fitfclust
res <- funcit(data=data, clusters=clusters,
              methods="fitfclust",
              k=4)

##plot including the centers
plotFuncy(data, col=clusters, ctr=Center(res), showLegend=TRUE, legendPlace="topleft")

}

\keyword{plot}
