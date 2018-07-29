%
%  Copyright (C) 2011-2015 Christina Yassouridis
% 
%

\name{makeCommonTime}
\alias{makeCommonTime}

\title{
Create appropriate evaluation time points.  
}
\description{
 Create an appropriate grid of common time points based on
 the available time points from the curves. 
}
\usage{
makeCommonTime(data, timeNr, plot = TRUE)
}

\arguments{
  \item{data}{
    Dataset in general format Format1 (see \code{\link{funcit}}).
}
  \item{timeNr}{
Number of time points. 
}
  \item{plot}{
If \code{TRUE}, original and new time points are plotted.
}
}
\details{
The evaluation time points of all curves together are clustered and new
  time points are represented by the cluster centers. 

}
\value{
A vector of new time points. 
}

\author{
Christina Yassouridis
}

\examples{
##sample an irregular dataset
set.seed(2804)
ds <- sampleFuncy(obsNr=100, k=4, timeNrMin=5, timeNrMax=10, reg=FALSE)
makeCommonTime(Data(ds), 10)
}

\references{
    Christina Yassouridis and Dominik Ernst and Friedrich Leisch. 
    Generalization, Combination and Extension of Functional Clustering Algorithms: 
    The R Package funcy. Journal of Statistical Software. 85 (9). 1--25. 2018
}
