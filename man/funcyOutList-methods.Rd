%
%  Copyright (C) 2011-2015 Christina Yassouridis
% 
%

\name{funcyOutList-methods}
\Rdversion{0.8}
\docType{methods}
\title{Methods for the output of \code{\link{funcit}}, an object of
  class \code{\link{funcyOutList}}. }
\alias{[[,funcyOutList,ANY,missing-method}
\alias{calcTime}
\alias{calcTime,funcyOutList-method}
\alias{Center}
\alias{Center,funcyOutList-method}
\alias{Cluster}
\alias{Cluster,funcyOutList-method}
\alias{props}
\alias{props,funcyOutList-method}
\alias{randIndex,funcyOutList,ANY-method}
\alias{summary,funcyOutList-method}
\description{Apply a function to \code{\link{funcyOutList}} which is the
  result from calling the function \code{\link{funcit}}.}
\usage{
\S4method{calcTime}{funcyOutList}(object)
\S4method{Center}{funcyOutList}(object)
\S4method{Cluster}{funcyOutList}(object)
\S4method{props}{funcyOutList}(object)
\S4method{randIndex}{funcyOutList,ANY}(x)
\S4method{summary}{funcyOutList}(object)
}

\arguments{
  \item{x, object}{object \code{\link{funcyOutList}} as a result of
    function \code{\link{funcit}}.}
}

\value{
  \describe{
    \item{calcTime}{Numeric matrix of dimension nr_methods x 5. Calculation times for the different methods.}
    \item{Cluster}{Numeric matrix of dimension nr_curves x
  methods. Cluster results for all methods.}
    \item{Center}{List of matrices of cluster centers for the different
  methods. Centers are stored in columns.}
    \item{props}{data.frame, proportion of the clusters for all methods. Can include NAs if a method reduced cluster number.}
    \item{randIndex}{\code{signature(x="funcyOutList", y="missing")}: Quadratic matrix of Rand indices showing the similarity between the methods. If true cluster membership was given, correct
      classification on diagonal.}
    \item{summary}{ Summary showing method call, cluster proportions, Rand indices and
      calculation time.}
  }
}

\author{Christina Yassouridis}
\keyword{calcTime, Center, Cluster, props, randIndex, summary}
