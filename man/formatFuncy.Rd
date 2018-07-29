%
%  Copyright (C) 2011-2015 Christina Yassouridis
% 
%

\name{formatFuncy}
\alias{formatFuncy}
\alias{formatFuncy,list,character-method}
\alias{formatFuncy,matrix,character-method}

\title{Switch between different data formats.}
\description{Datasets can be stored in different formats.
  \describe{
    \item{"Format1":}{General format for regular and irregular data. One long matrix with curve-ID in first column, curve evaluations in
      second and time points in third column.}
    \item{"Format2":}{Format for regular data only. Matrix of dimension nr_time x nr_curves.}
    \item{"Format3":}{Format for regular and irregular datasets. List of three
      matrices, all with same dimension.
      One evaluation matrix \code{Yin} where curves are stored
      in rows. One matrix of time points \code{Tin}. One incidence matrix \code{isobs} with entry 1, if curve has evaluation at this time point.}
  }
}
  \usage{
\S4method{formatFuncy}{list,character}(data, format="Format1", regTime=NULL)
\S4method{formatFuncy}{matrix,character}(data, format="Format1", regTime=NULL)
}

\arguments{
  \item{data}{
    Data in format \code{"Format1"}, \code{"Format2"} or \code{"Format3"}.
  }

  \item{format}{
    Format to transform data to. One of \code{"Format1"} or \code{"Format3"}.\cr
    \code{"Format2"} is possible only if dataset in \code{data} is regular.
  }
  \item{regTime}{
    Optional vector of time points if original data is in \code{"Format2"}.
  }
}
\details{
  Data, especially for irregular time points is often stored in different
  ways. To switch back and forth between data formats, formatFuncy can
    be used. For regular datasets in format \code{"Format2"}, an optional
    vector of evaluation time points can be given as
    \code{regTime}.
}
\value{
  Numeric matrix if \code{format}=\code{"Format1"} or
  \code{format}=\code{"Format2"}. If \code{format}=\code{"Format3"}, a
  list of the three matrices \code{Yin}, \code{Tin}, \code{isobs} and
  two vectors \code{N} and \code{t_all}. \code{N} stores the number of
  time points for each curve and \code{t_all} is a vector of unique
  time points (see \code{Arguments}).
}

\author{
Christina Yassouridis
}


\examples{
##Generate regular dataset
set.seed(2005)
ds <- sampleFuncy(obsNr=100, k=6, timeNr=20, reg=TRUE)
Data(ds)

##Format dataset to Format1
newdat <-formatFuncy(data=Data(ds), format="Format1")
newdat

##Back to matrix out of Format1
formatFuncy(newdat, format="Format2")

##To Format3
formatFuncy(newdat, format="Format3")

##Generate irregular dataset
set.seed(2005)
ds <- sampleFuncy(obsNr=100, k=5, timeNrMin=5, timeNrMax=10, reg=FALSE)

res <- formatFuncy(Data(ds), format="Format3", reg=FALSE)
res

##Back to Format1
formatFuncy(data=res, format="Format1", reg=FALSE)
}

\references{
    Christina Yassouridis and Dominik Ernst and Friedrich Leisch. 
    Generalization, Combination and Extension of Functional Clustering Algorithms: 
    The R Package funcy. Journal of Statistical Software. 85 (9). 1--25. 2018
    }