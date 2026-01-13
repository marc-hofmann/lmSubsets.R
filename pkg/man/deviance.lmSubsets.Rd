\name{deviance.lmSubsets}


\alias{deviance.lmSubsets}
\alias{deviance.lmSelect}


\title{Extract the deviance from a subset regression}


\description{

  Return the deviance for the specified submodels.

}


\usage{
\method{deviance}{lmSubsets}(object, size, best = 1, ..., na.rm = TRUE, drop = TRUE)

\method{deviance}{lmSelect}(object, best = 1, ..., na.rm = TRUE, drop = TRUE)
}


\arguments{

  \item{object}{\code{"lmSubsets"}, \code{"lmSelect"}---a subset
    regression}

  \item{size}{\code{integer[]}---the submodel sizes}

  \item{best}{\code{integer[]}---the submodel positions}

  \item{...}{ignored}

  \item{na.rm}{\code{logical}---if \code{TRUE}, remove \code{NA}
    entries}

  \item{drop}{\code{logical}---if \code{TRUE}, simplify structure}

}


\value{

  \code{double[]}, \code{"data.frame"}---the submodel deviances

}


\seealso{

  \itemize{

    \item{\code{\link[=lmSubsets]{lmSubsets()}} for all-subsets
      regression}

    \item{\code{\link[=lmSelect]{lmSelect()}} for best-subset
      regression}

    \item{\code{\link[=deviance]{deviance()}} for the S3 generic}

  }

}
