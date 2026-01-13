\name{variable.names.lmSubsets}


\alias{variable.names.lmSubsets}
\alias{variable.names.lmSelect}


\title{Extract variable names from a subset regression}


\description{

  Return the variable names for the specified submodels.

}


\usage{
\method{variable.names}{lmSubsets}(object, size, best = 1, ..., na.rm = TRUE, drop = TRUE)

\method{variable.names}{lmSelect}(object, best = 1, ..., na.rm = TRUE, drop = TRUE)
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

  \code{logical[,]}, \code{"data.frame"}---the variable names

}


\seealso{

  \itemize{

    \item{\code{\link[=lmSubsets]{lmSubsets()}} for all-subsets
      regression}

    \item{\code{\link[=lmSelect]{lmSelect()}} for best-subset
      regression}

    \item{\code{\link[=variable.names]{variable.names()}} for the S3
      generic}

  }

}
