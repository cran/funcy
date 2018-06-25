%
%  Copyright (C) 2011-2015 Christina Yassouridis
% 
%

\name{funcit}
\alias{funcit}

\title{Functional Cluster Analysis}
\description{
  Main function for clustering functional data according to one or several of seven algorithms.
}
\usage{
funcit(data, k, methods=c("fitfclust","distclust", "iterSubspace",
       "funclust", "funHDDC", "fscm", "waveclust"), seed=NULL, regTime=NULL,
       clusters=NULL, funcyCtrl=NULL, fpcCtrl=NULL, parallel=FALSE,
       save.data=TRUE, ...)
}
\arguments{
  \item{data}{
    Data in format "Format1" or format "Format2" (see \code{\link{formatFuncy}}).
  }
  \item{k}{
    Number of clusters.
  }
  \item{methods}{
    \describe{
      \item{"fitfclust":}{Model based cluster algorithm - based on a functional
	  mixed mixture model. Allows irregular measurements, eigenbasis possible.}
	\item{"distFPCA":}{
	Cluster algorithm - based on a distance measure. Allows irregular
	measurements, eigenbasis possible.}
      \item{"iterSubspace":}{Model based cluster algorithm - based on a subspace
	projection. Allows irregular measurements, eigenbasis possible,
	dimension between clusters can vary.}
      \item{"funclust":}{Model based cluster algorithm - based on a functional mixed mixture model.}
      \item{"funHDDC":}{Model based cluster algorithm - based on a
	functional mixed mixture model. Dimension between clusters can
vary.}
     \item{"fscm":}{Model based cluster algorithm - based on a functional mixed
	mixture model. Curves can dependent on location. A matrix
	\code{location} is then an optional input parameter (see
	Details).}
      \item{"waveclust":}{Model based cluster algorithm - based on a
	functional mixed mixture model. Wavelet basis is the only possible.}
    }
For a detailed description of the methods please see the references. 
  }
  \item{seed}{
    Seed for initial clustering. See \code{\link{funcyCtrl}}.
  }  
  \item{regTime}{
    If data is in "Format2", optional vector representing the
    time points (see \code{\link{formatFuncy}}). If regTime=NULL and \code{format}="Format2", equidistant time points from 1 to number of curves are used. 
  }
  \item{clusters}{
    Optional vector of true cluster labels.
  }
  \item{funcyCtrl}{
    A control object of class \code{\link{funcyCtrl}}. If a model based
    clustering algorithm is used, further parameters can be specified by
    using the extended class \code{\link{fpcCtrlMbc}}.
  }
  \item{fpcCtrl}{
    A control object of class \code{\link{fpcCtrl}}. Only used for
    eigenbasis calculation (\code{baseType}=\code{"eigenbasis"} in \code{\link{funcyCtrl}}).
  }
  \item{parallel}{
    If \code{TRUE}, package \pkg{parallel} is used for parallel
    computing.
  }
  \item{save.data}{
    Save a copy of the \code{data} in the return object? Must be set to TRUE
    in order to use plot function \code{\link{plot}}. 
  }
  \item{\dots}{
    Additional optional model specific parameters. Works only if exactly one method
    is called in \code{methods}. The parameters are the
    following: 
    \describe{
      \item{"fitfclust"}{
	\describe{
	  \item{}{}
 	  \item{p:}{
 	    Rank of the covariance matrix \eqn{\Gamma}, must be at
 	    least \code{dimBase}.}
 	  \item{pert:}{
 	    Adds a ridge term to the least squares fit, helps if
 	    only few observations per curve were registered.}
	}
      }

      \item{"distclust"}{
	\describe{
	  \item{}{}
	  \item{method:}{
	    One of \code{"hclust"} or \code{"pam"} specifying how
	    distance matrix is processed.}
	}
      }
      
      \item{"iterSubspace"}{
	\describe{
	  \item{}{}
	  \item{simplif:}{\code{FALSE}, if curve affiliation is tested
	    again by projecting the curve onto the current subspace
	    created without the current curve (leave-one-out-curve-estimation).}
	}
      }
      
      \item{"funclust"}{
	\describe{
	  \item{}{}
	  \item{nbInit:}{The number of small-EM used to determine the initialization of
            the main EM-like algorithm.}
	  \item{nbIterInit:}{The maximum number of iterations for each
	    small-EM.}
	}
      }
      
      \item{"funHDDCWrapper"}{
	\describe{
	  \item{}{}
	  \item{model:}{
	    The chosen model among "AkjBkQkDk", "AkjBQkDk",\cr
	    "AkBkQkDk","AkBQkDk","ABkQkDk","ABQkDk". \cr
	    See (Bouveyron & Jacques, 2011) for details.}
      }}
      
      \item{"fscm"}{
	\describe{
	  \item{}{}
	  \item{location:}{
	    A two-dimensional matrix of the curve
	    locations (coordinates).}
	  \item{knn:}{
	    Number of neighbors each curve depends on.}
	  \item{useCode:}{
	    \code{"R"} or {"C"}. If C is installed, a lot faster than
	    R.}
	  \item{verbose:}{\code{TRUE}, if number of iterations and
	sigma, theta and f are to be printed.}
      }}
      
      \item{"waveclust"}{
	\describe{
	  \item{}{}
	  \item{gamma:}{
	    One of \code{"group"}, \code{"scale.location"},
	    \code{"group.scale.location"} or \code{"constant"}.}
	  \item{init:}{
	    One of \code{"rEM"} or \code{"SEM"} for random or
	    stochastic EM.}
	  \item{plotLoglik:}{
	    \code{TRUE}, if log-likelihood is to be plotted. 
	    }
      }}
    }
  }
}


\details{
  \code{funcit} is the core function to execute one or more methods to cluster functional
  data. Functional data can be measured on a regular or on an irregular
  grid. While for regular datasets, all curves are measured on the same
  time points, for irregular datasets, number or/and location of time
  points can differ (see \code{\link{formatFuncy}} for different formats). Only algorithms \code{"fitfclust"},\code{"distclust"} and
  \code{"iterSubspace"} are applicable to irregular datasets. 
  All methods are based on the projection of the curves onto a
  basis defined in \code{\link{funcyCtrl}} and building mixed effects
  models of the basis coefficients.
}


\value{
  Returns an object of class \code{\link{funcyOutList}}.
}
	 

\references{
  \describe{
    \item{"fitfclust":}{Gareth James and Catherine A. Sugar. Clustering for Sparsely
 Sampled Functional Data. Journal of the American Statistical
 Association. 98 (462). 297--408. 2003}
    \item{"distFPCA":}{Jie Peng and Hans-Georg Mueller. Distance-based clustering of sparsely observed stochastic processes, with applications to online auctions. The Annals
 of Applied Statistics. 2 (3). 1056--1077, 2008}
 \item{"iterSubspace":}{Chiou Jeng-Min and Pai-Ling Li. Functional clustering and
 identifying substructures of longitudinal data. Journal of the
 Royal Statistical Society: Series B. 69 (4). 679--699. 2007}
 \item{"waveclust":}{Madison Giacofci and Sophie Lambert-Lacroix and Guillemette Marot
 and Franck Picard. Wavelet-based clustering for mixed-effects
  functional models in high dimension. Biometrics. 69. 31--40. 2011}
 \item{"fscm":}{Nicoleta Serban and Huijing Jiang.Clustering Random Curves Under Spatial Interdependence With Application to Service Accessibility. Technometrics. 54 (2). 108--119. 2012 }
 \item{"funclust":}{Julien Jacques and Cristian Preda. Funclust: a curves clustering
 method using functional random variables density
 approximation. Neurocomputing. 112. 164–171. 2013}
 \item{"funHDDC":}{Charles Bouveyron and Julien and Jacques. Model-based clustering
 of time series in group-specific functional subspaces. Advances
 in Data Analysis and Classification. 5 (4). 281--300. 2011}
    }

}
\author{Christina Yassouridis}


\examples{
##Cluster the data with methods for regular sets
##Sample a regular dataset
set.seed(2804)
ds <- sampleFuncy(obsNr=50, k=4, timeNr=8, reg=TRUE)

##Cluster the functions with all available methods. 
res <- funcit(data=Data(ds), clusters=Cluster(ds),
              methods=c(1,2,3,4), seed=2404,
              k=4)
summary(res)
Cluster(res)

##Additional method specific parameters for method fitfclust
res <- funcit(data=Data(ds), clusters=Cluster(ds), methods="fitfclust", seed=2405,
              k=4, p=5, pert=0)


##Cluster the data with methods for irregular sets
##Sample an irregular dataset
set.seed(2804)
ds <- sampleFuncy(obsNr=50, k=4, timeNrMin=4, timeNrMax=8,
                  reg=FALSE)
data <- Data(ds)
clusters <- Cluster(ds)

res <- funcit(data=data, clusters=clusters,
              methods=c("fitfclust","distclust", "iterSubspace"), seed=2406,
              k=4, parallel=TRUE)

summary(res)
Cluster(res)
plot(res)

##Two reallife examples
\dontrun{
data("genes")
data <- genes$data
clusters <- genes$clusters

##Cluster the functions with all available methods. 
res <- funcit(data=data, clusters=clusters,
              methods=1:7, seed=2404,
              k=4)
summary(res)
Cluster(res)
}

\dontrun{
data("electricity")
res <- funcit(data=electricity, methods=c("fitfclust","distclust",
"waveclust"), seed=2406, k=5, parallel=TRUE)
plot(res, legendPlace="topleft")
}
}

\keyword{cluster}
