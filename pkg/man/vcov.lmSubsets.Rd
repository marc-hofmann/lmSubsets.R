\name{vcov.lmSubsets}


\alias{vcov.lmSubsets}
\alias{vcov.lmSelect}


\title{Extract the variance-covariance matrix from a subset regression}


\description{

  Return the variance-covariance matrix for the specified submodel.

}


\usage{
\method{vcov}{lmSubsets}(object, size, best = 1, ...)

\method{vcov}{lmSelect}(object, best = 1, ...)
}


\arguments{

  \item{object}{\code{"lmSubsets"}, \code{"lmSelect"}---a subset
    regression}

  \item{size}{\code{integer}---the submodel size}

  \item{best}{\code{integer}---the submodel position}

  \item{...}{ignored}

}


\value{

  \code{double[,]}---the variance-covariance matrix

}


\seealso{

  \itemize{

    \item{\code{\link[=lmSubsets]{lmSubsets()}} for all-subsets
      regression}

    \item{\code{\link[=lmSelect]{lmSelect()}} for best-subset
      regression}

    \item{\code{\link[=vcov]{vcov()}} for the S3 generic}

  }

}
