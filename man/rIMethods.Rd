%
%  Copyright (C) 2011-2015 Christina Yassouridis
% 
%


\name{rIMethods}
\alias{rIMethods}
\title{
Rand index for cluster configurations of different methods. 
}
\description{
Calculates a matrix of Rand indices for the configurations of all
cluster method combinations.
}
\usage{
rIMethods(methodNames = NULL, cls, trueCluster = NULL)
}

\arguments{
  \item{methodNames}{
    Character Vector of names for the methods.
  }
  \item{cls}{
    Numeric or character matrix of cluster outcomes for different
    methods.
  }
  \item{trueCluster}{
    Vector of true clusters if knows. 
  }
}
\details{
If \code{trueCluster} is given, the Rand index between the true cluster
outcome and the clusters calculated by the methods is on diagonal. 
}
\value{
A matrix of Rand indices showing the similarity of the methods. 
}

\author{
Christina Yassouridis
}

\examples{
set.seed(1234)
ds <- sampleFuncy(obsNr=80, timeNr=20, reg=TRUE, k=4, sd=.4)
res1 <- funcit(methods=1:3, data=Data(ds), k=4, clusters=Cluster(ds))

cls <- Cluster(res1)
rIMethods(methodNames=c("method1","method2","method3"),cls=cls, trueCluster=Cluster(ds))
}

\keyword{Rand index}
