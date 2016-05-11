\name{plot-methods}
\docType{methods}
\alias{plot,funcyOutList,missing-method}
\alias{plot,funcyOut,missing-method}
\alias{plot,funcyOut,ANY-method}
\alias{plot,funcyOutMbc-fitfclust,missing-method}
\alias{plot,funcyOutMbc-fscm,missing-method}


\title{Methods for Function \code{plot} in Package \pkg{funcy}}
\description{
  Plots clustered curves and/or cluster centers and other results. 
}
\usage{
\S4method{plot}{funcyOut,missing}(x, y, data, type="all",
     showLegend=TRUE, legendPlace="bottomleft", main, \dots)

\S4method{plot}{funcyOutList,missing}(x, y, data=NULL,
     select=NULL, type="all", showLegend=TRUE,
     legendPlace="bottomleft", main, \dots)
}
\arguments{
  \item{x}{An object of class \code{"funcyOut"} or \code{"funcyOutList"}.}
  \item{y}{Not used.}
  \item{data}{Data to include in plot. If the cluster object \code{x} is
    of type \code{"funcyOutList"} that was created with
 \code{save.data=TRUE}, then these are used by default.}
  \item{select}{Select the methods, you want to generate the plot
    for.}
  \item{type}{Plot type, see details.}
  \item{showLegend}{If \code{TRUE}, cluster legend is shown.}
  \item{legendPlace}{Legend placement.}
  \item{main}{Plot title, can be missing.}
  \item{\dots}{Further plotting parameters}
}
\details{
  If data was clustered by \code{\link{funcit}} with
  \code{save.data}=\code{TRUE}, different plots can be used. Some plots
  are available for all methods, others depend on method which was
  used. The plot types are listed below. 
  If method specific plots are used, method must be extracted by
  \code{select}=\code{method name}, see examples. 
  \describe{
    \item{"all methods":}{
      \describe{
	\item{}{}
	\item{all:}{
	  Plots data and cluster centers.}
	\item{centers:}{
	  Plots only cluster centers.}
	\item{shadow:}{
	  Creates a shadow plot (see function \code{shadow} in package
	  \pkg{flexclust} - Leisch 2010).}
	\item{dist2centers:}{
	  Multiple plots for each cluster. Thickness of lines corresponds
	  to the proximity to the cluster centers. Thicker lines means
	  curve is closer to its center.}
	\item{fpc:}{
	  Only if \code{baseType}=\code{"eigenbasis"} in \code{funcyCtrl}.
	  Plots the smoothed mean function, covariance matrix and
	  eigenbasis.} 
      }
    }
  \item{"fitfclust":}{
    \describe{
      \item{}{}
      \item{discrim}{
	Plots discriminant functions to show the time points of
	maximum discrimination between clusters (see James2003).}
      \item{conf}{
	Plots confidence intervals for the curves.}
    }
  }
      
  \item{"fscm":}{
    \describe{
      \item{}{}
      \item{overview:}{
	Plots curve locations, temporal trends and overall trends (see Serban2012).
	For the spatial coefficients, dots are colored according to
	spatial dependency from yellow to blue. Darker dots mean
	stronger depencency.}
    }
  }
  
}
}

\references{
    Friedrich Leisch. Neighborhood graphs, stripes and shadow plots for
    cluster visualization. Statistics and
    Computing. 20(4). 457--469. 2010

     Gareth James and Catherine A. Sugar. Clustering for Sparsely
     Sampled Functional Data. Journal of the American Statistical
     Association. 98 (462). 297--408. 2003

     Nicoleta Serban and Huijing Jiang.Clustering Random Curves Under Spatial Interdependence With Application to Service Accessibility. Technometrics. 54 (2). 108--119. 2012
   }


   \examples{
set.seed(2804)
ds <- sampleFuncy(obsNr=60, k=4, timeNrMin=5, timeNrMax=10, reg=FALSE)
data <- Data(ds)
clusters <- Cluster(ds)
res <- funcit(data=data, clusters=clusters,
              methods=c("fitfclust","distclust", "iterSubspace") ,
              k=4, parallel=TRUE)
plot(res)
plot(res, select="fitfclust", type="conf")
plot(res, select="fitfclust", type="discrim")
plot(res, select="distclust", type="shadow")
}
   
\keyword{plot}
