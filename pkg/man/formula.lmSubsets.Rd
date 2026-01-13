\name{formula.lmSubsets}


\alias{formula.lmSubsets}
\alias{formula.lmSelect}


\title{Extract a formula from a subset regression}


\description{

  Return the formula for the specified submodel.

}


\usage{
\method{formula}{lmSubsets}(x, size, best = 1, ...)

\method{formula}{lmSelect}(x, best, ...)
}


\arguments{

  \item{x}{\code{"lmSubsets"}, \code{"lmSelect"}---a subset regression}

  \item{size}{\code{integer}---the submodel size}

  \item{best}{\code{integer}---the submodel position}

  \item{...}{ignored}

}


\value{

  \code{"formula"}---the submodel formula

}


\seealso{

  \itemize{

    \item{\code{\link[=lmSubsets]{lmSubsets()}} for all-subsets
      regression}

    \item{\code{\link[=lmSelect]{lmSelect()}} for best-subset
      regression}

    \item{\code{\link[=formula]{formula()}} for the S3 generic}

  }

}
