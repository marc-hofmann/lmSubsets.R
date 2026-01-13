\name{BIC.lmSubsets}


\alias{BIC.lmSubsets}
\alias{BIC.lmSelect}


\title{Extract BIC values from a subset regression}


\description{

  Evaluate the Bayesian information criterion (BIC) for the specified
  submodels.

}


\usage{
\method{BIC}{lmSubsets}(object, size, best = 1, ..., na.rm = TRUE, drop = TRUE)

\method{BIC}{lmSelect}(object, best = 1, ..., na.rm = TRUE, drop = TRUE)
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

  \code{double[]}---the BIC values

}



\seealso{

  \itemize{

    \item{\code{\link[=lmSubsets]{lmSubsets()}} for all-subsets
      regression}

    \item{\code{\link[=lmSelect]{lmSelect()}} for best-subset
      regression}

    \item{\code{\link[=BIC]{BIC()}} for the S3 generic}

  }

}
