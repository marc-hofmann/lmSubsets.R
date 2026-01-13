\name{AIC.lmSubsets}


\alias{AIC.lmSubsets}
\alias{AIC.lmSelect}


\title{Extract AIC values from a subset regression}


\description{

  Evaluate Akaike's information criterion (AIC) for the specified
  submodels.

}


\usage{
\method{AIC}{lmSubsets}(object, size, best = 1, ..., k = 2, na.rm = TRUE, drop = TRUE)

\method{AIC}{lmSelect}(object, best = 1, ..., k = 2, na.rm = TRUE, drop = TRUE)
}


\arguments{

  \item{object}{\code{"lmSubsets"}, \code{"lmSelect"}---a subset
    regression}

  \item{size}{\code{integer[]}---the submodel sizes}

  \item{best}{\code{integer[]}---the submodel positions}

  \item{...}{ignored}

  \item{k}{\code{double}---the penalty per model parameter}

  \item{na.rm}{\code{logical}---if \code{TRUE}, remove \code{NA}
    entries}

  \item{drop}{\code{logical}---if \code{TRUE}, simplify structure}

}


\value{

  \code{double[]}---the AIC values

}


\seealso{

  \itemize{

    \item{\code{\link[=lmSubsets]{lmSubsets()}} for all-subsets
      regression}

    \item{\code{\link[=lmSelect]{lmSelect()}} for best-subset
      regression}

    \item{\code{\link[=AIC]{AIC()}} for the S3 generic}

  }

}
