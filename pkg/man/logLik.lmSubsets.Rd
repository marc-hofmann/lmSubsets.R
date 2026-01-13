\name{logLik.lmSubsets}


\alias{logLik.lmSubsets}
\alias{logLik.lmSelect}


\title{Extract the log-likelihood from a subset regression}


\description{

  Return the log-likelihood of the the specified submodels.

}


\usage{
\method{logLik}{lmSubsets}(object, size, best = 1, ..., na.rm = TRUE, drop = TRUE)

\method{logLik}{lmSelect}(object, best = 1, ..., na.rm = TRUE, drop = TRUE)
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

  \code{double[]}---the log-likelihoods

}


\seealso{

  \itemize{

    \item{\code{\link[=lmSubsets]{lmSubsets()}} for all-subsets
      regression}

    \item{\code{\link[=lmSelect]{lmSelect()}} for best-subset
      regression}

    \item{\code{\link[=logLik]{logLik()}} for the S3 generic}

  }

}
