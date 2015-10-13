%
%  Copyright (C) 2011-2015 Christina Yassouridis
% 
%

\name{label2lowerk}
\alias{label2lowerk}

\title{
Relabeling of clusters. 
}
\description{
Relabels arbitrary cluster labels for k classes to 1 to k.
}
\usage{
label2lowerk(cluster)
}

\arguments{
  \item{cluster}{
    Original cluster labels.
  }
}

\value{
New cluster labels.
}

\author{
  Christina Yassouridis
}

\examples{
##Generating cluster labels for 4 clusters
cl <- rep(sample(2:10, 4), sample(3:5,4, replace=TRUE))
cl

##Labels them to 1:4
label2lowerk(cl)
}
\keyword{relabel}
