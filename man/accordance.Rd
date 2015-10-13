%
%  Copyright (C) 2011-2015 Christina Yassouridis
% 
%

\name{accordance}
\alias{accordance}
\title{
Accordance for cluster outcomes from different methods. 
}
\description{
Votes for the cluster with maximum agreement and shows relative number of
methods that voted for this cluster.  
}
\usage{
accordance(cls, ctrs=NULL, relabel=FALSE)
}

\arguments{
  \item{cls}{
    Numeric or character matrix of cluster outcomes for different
    methods. Cluster outcomes for each method are in columns. 
}
\item{ctrs}{
  If \code{relabel}=\code{TRUE}, a list of centers for the different
  methods, if \code{NULL}, clusters must have been labeled in an
  appropriate way before.
}
\item{relabel}{
  If \code{TRUE}, clusters are relabeled according to maximum agreement
  between their center curves. \code{ctrs} must not be \code{NULL} in
  this case. 
}
}
\details{
  If \code{relabel}=\code{FALSE}, clusters must have been relabeled in an appropriate way
  before, so that cluster labels representing the same or the most
  similar clusters are identical.
  If \code{relabel}=\code{TRUE}, centers must be given as list to
  \code{ctrs}. Each list entry consists of a matrix with cluster centers stored in columns.
}
\value{
 \item{votedCluster}{Cluster with maximum agreement between the
  methods.}
\item{accordance}{Percentage of methods that voted for the cluster.}

}

\author{
Christina Yassouridis
}


\examples{
##Generate dataset
k <- 3
set.seed(2511)
ds <- sampleFuncy(obsNr=30, timeNr=10, reg=TRUE, k=k, sd=.4)
 
##Cluster dataset
res1 <- funcit(methods=1:3, data=Data(ds), k=4, clusters=Cluster(ds))

##Retrieve clusters and centers
cl <- Cluster(res1)
ctr <- Center(res1)

accordance(cl, ctr, relabel=TRUE)
}
