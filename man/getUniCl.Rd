\name{getUniCl}
\alias{getUniCl}
\title{
  Unique cluster labels
}
\description{
  For repeated measurements on the same ID, build unique cluster labels according to the vector of IDs.
}
\usage{
getUniCl(id, clusters, reduce=TRUE)
}

\arguments{
  \item{id}{
    Vector of IDs.
}
  \item{clusters}{
    Vector of cluster labels. 
  }
  \item{reduce}{
    \code{TRUE} if cluster labels shall be reduced according to unique IDs.
    \code{FALSE} if cluster labels shall be duplicated according to IDs. 

  }
}
\details{
  Data might have been stored in 4 columns: curveID, curve evaluations,
  time points and cluster labels. Cluster labels were therefore repeated
  for each curve evaluation point. Method \code{\link{funcit}} accepts
  dataset only in formats "Format1" and "Format2" and an optional vector
  of cluster true labels \code{clusters} of length nr_curves. \code{getUniCl} can be
  applied to columns curveID and repeated labels to reduce them
  to the number of curves. 
}
\value{
A vector of either reduced (\code{reduce}=\code{TRUE}) or duplicated (\code{reduce}=\code{FALSE}) cluster labels. 
}

\author{
Christina Yassouridis
}

\examples{
##Generate dataset
nr_time <- sample(1:5, 100,  replace=TRUE)
clusters <- sample(1:4, 100, replace=TRUE)
IDs <- rep(1:100, nr_time)
cls <- rep(clusters,nr_time)
IDs
cls

##Get reduced cluster vector according to IDs
unicl <- getUniCl(IDs,cls)
unicl

##Get original cluster vector from reduced one
dupcl <- getUniCl(IDs,unicl, reduce=FALSE)
dupcl
}

