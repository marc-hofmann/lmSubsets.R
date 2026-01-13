\name{fitted.lmSubsets}


\alias{fitted.lmSubsets}
\alias{fitted.lmSelect}


\title{Extract the fitted values from a subset regression}


\description{

  Return the fitted values for the specified submodel.

}


\usage{
\method{fitted}{lmSubsets}(object, size, best = 1, ...)

\method{fitted}{lmSelect}(object, best = 1, ...)
}


\arguments{

  \item{object}{\code{"lmSubsets"}, \code{"lmSelect"}---a subset
    regression}

  \item{size}{\code{integer}---the submodel size}

  \item{best}{\code{integer}---the submodel position}

  \item{...}{ignored}

}


\value{

  \code{double[]}---the fitted values

}


\seealso{

  \itemize{

    \item{\code{\link[=lmSubsets]{lmSubsets()}} for all-subsets
      regression}

    \item{\code{\link[=lmSelect]{lmSelect()}} for best-subset
      regression}

    \item{\code{\link[=fitted]{fitted()}} for the S3 generic}

  }

}
