%
%  Copyright (C) 2011-2015 Christina Yassouridis
% 
%

\name{relabel}
\alias{relabel}
\alias{relabelMethods}

\title{
Relabel cluster IDs. 
}
\description{
Relabels cluster IDs of two or more cluster configurations according to
the minimal distance between their centers. 
}
\usage{
relabel(cl1, cl2, ctr1, ctr2)
relabelMethods(methodNames=NULL, cls, ctrs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cl1}{
    Cluster IDs of the first configuration. 
  }
  \item{cl2}{
    Cluster IDs of the second configuration. 
  }
  \item{ctr1}{
    Numeric matrix of cluster centers of the first configuration. 
  }
  \item{ctr2}{
    Numeric matrix of cluster centers of the second configuration. 
  }
  \item{methodNames}{
    Character vector of names for the different cluster methods. 
  }
  \item{cls}{
   Numeric or character matrix of cluster outcomes for different
    methods.
  }
  \item{ctrs}{
    List of centers for the different methods.
  }
}
\details{
The two configurations for the method \code{relabel} do not have to have the same number of
observations neither the same number of clusters. The configuration with
less observations has to be put on place \code{cl1}.\cr
\code{relabelMethods} can be used for more than two
configurations. Cluster outputs must therefore be saved in a matrix \code{cls}
and thus have the same number of observations. The relabeling works only
correctly if the number of classes is the same. The relabeling of the
methods follows the following scheme:
\enumerate{
  \item
  Methods are sorted in that way that one of the two most similars is on first
  place.
  \item
  Methods are successively added in the order of the highest similarity to one of
  the already added methods.
  \item
  Once the order is fixed, the methods are relabeled after the ones they
  are most similar to.
  }
}
\value{
  For method \code{relabel}:
  \item{recode}{Recoding scheme of the second cluster labels \code{cl2}.}
  \item{cluster}{New cluster labels for \code{cl2}.}
  \item{centers}{Cluster centers \code{ctrs2} in the new order. }
  
  For method \code{relabelMethods}:
  \item{allClusters}{Matrix of new cluster labels.}
  \item{allCenters}{List of cluster centers in the new order.}
  \item{fromTo}{Recoding scheme of the methods.}
}

\author{
Christina Yassouridis
}

\examples{
##Generate dataset
k <- 6
set.seed(2004)
ds <- sampleFuncy(obsNr=50, timeNrMin=3, timeNrMax=10, reg=FALSE, k=k, sd=.5)

##Cluster with different methods
res1 <- funcit(methods="fitfclust", data=Data(ds), k=k, reg=FALSE)
res2 <- funcit(methods="iterSubspace", data=Data(ds), k=k, reg=FALSE)
res3 <- funcit(methods="distclust", data=Data(ds), k=k, reg=FALSE)

##Relabel two configurations
relabel(Cluster(res3),Cluster(res1),Center(res3),Center(res1))

##Make matrix of clutser configurations
cls <- cbind(Cluster(res1),Cluster(res2),Cluster(res3))
##Make list of Centers
ctrs <- list(Center(res1), Center(res2), Center(res3))

##Relabel cluster configurations
rel <- relabelMethods(cls=cls, ctrs=ctrs)

##Compare
cls
rel$allClusters

}
