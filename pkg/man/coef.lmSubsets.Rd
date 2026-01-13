\name{coef.lmSubsets}


\alias{coef.lmSubsets}
\alias{coef.lmSelect}


\title{Extract the ceofficients from a subset regression}


\description{

  Return the coefficients for the specified submodels.

}


\usage{
\method{coef}{lmSubsets}(object, size, best = 1, ..., na.rm = TRUE, drop = TRUE)

\method{coef}{lmSelect}(object, best = 1, ..., na.rm = TRUE, drop = TRUE)
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

  \code{double[,]}, \code{"data.frame"}---the submodel coefficients

}


\seealso{

  \itemize{

    \item{\code{\link[=lmSubsets]{lmSubsets()}} for all-subsets
      regression}

    \item{\code{\link[=lmSelect]{lmSelect()}} for best-subset
      regression}

    \item{\code{\link[=coef]{coef()}} for the S3 generic}

  }

}
