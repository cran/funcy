\name{regFuncy}
\alias{regFuncy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Converts irregular data into regular.
}
\description{
Converts irregular data into regular by projecting it to a basis or
interpolating it. 
}
\usage{
regFuncy(data, timeNr = 10, method = "project", baseType = NULL, nbasis = 4, plot = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    Dataset in "Format1" (see \code{\link{funcit}}).
  }
  \item{timeNr}{
    Number of time points the regular dataset shall be evaluated at.
  }
  \item{method}{
    Method to tranform regularization with one of "project",
    "interpolate", "pace". See Details. 
  }
  \item{baseType}{
    Base type. Only used if method="project".
  }
  \item{nbasis}{
    Number of basis functions.
  }
  \item{plot}{
    Plot the result?
  }
}
\details{
  Data is either interpolated or projected to a basis. If method="pace",
  regularization is done via functional principal component analysis by
  conditional estimation (Müller2005). This is especially useful if
  the dataset is sparse. "pace" is implemented in the supplementary
  package "funcyOctave" which is loaded automatically but works only on
  Unix-platforms. 
  For all methods, curve evaluation takes place on time points calculated by
  \code{\link{makeCommonTime}}.
}
\value{
  \item{data}{Numeric matrix of "Format1" (see \code{\link{formatFuncy}}).}
  \item{time}{Vector of evualuation time points.}

}
\references{
  F. Yao and H.G. Müller and J.L. Wang. Functional data analysis for
  sparse longitudinal data. J. American Statistical
  Association. 100. 577--590. 2005
  URL: http://www.stat.ucdavis.edu/PACE/

}
\author{
Christina Yassouridis
}


\examples{
##Generate irregular dataset
set.seed(2705)
ds <- sampleFuncy(reg=FALSE, obsNr=100, timeNrMin=5, timeNrMax=10)


reg <- regFuncy(Data(ds), timeNr=10, baseType="splines",nbasis=5, method="project")

reg <- regFuncy(Data(ds),timeNr=10, method="interpolate")

\dontrun{
reg <- regFuncy(Data(ds), timeNr=10, baseType="eigenbasis", nbasis=5,
method="project")
}

}


\keyword{format}
