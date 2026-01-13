\name{residuals.lmSubsets}


\alias{residuals.lmSubsets}
\alias{residuals.lmSelect}


\title{Extract the residuals from all-subsets regression}


\description{

  Return the residuals for the specified submodel.

}


\usage{
\method{residuals}{lmSubsets}(object, size, best = 1, ...)

\method{residuals}{lmSelect}(object, best = 1, ...)
}


\arguments{

  \item{object}{\code{"lmSubsets"}, \code{"lmSelect"}---a subset
    regression}

  \item{size}{\code{integer}---the submodel size}

  \item{best}{\code{integer}---the submodel position}

  \item{...}{ignored}

}


\value{

  \code{double[]}---the residuals

}


\seealso{

  \itemize{

    \item{\code{\link[=lmSubsets]{lmSubsets()}} for all-subsets
      regression}

    \item{\code{\link[=lmSelect]{lmSelect()}} for best-subset
      regression}

    \item{\code{\link[=residuals]{residuals()}} for the S3 generic}

  }

}
