\name{sigma.lmSubsets}


\alias{sigma.lmSubsets}
\alias{sigma.lmSelect}


\title{Extract the residual standard deviation from a subset regression}


\description{

  Return the residual standard deviation for the specified submodels.

}


\usage{
\method{sigma}{lmSubsets}(object, size, best = 1, ..., na.rm = TRUE, drop = TRUE)

\method{sigma}{lmSelect}(object, best = 1, ..., na.rm = TRUE, drop = TRUE)
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

  \code{double[]}---the residual standard deviations

}


\seealso{

  \itemize{

    \item{\code{\link[=lmSubsets]{lmSubsets()}} for all-subsets
      regression}

    \item{\code{\link[=lmSelect]{lmSelect()}} for best-subset
      regression}

    \item{\code{\link[=sigma]{sigma()}} for the S3 generic}

  }

}
