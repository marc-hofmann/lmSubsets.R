\name{model.matrix.lmSubsets}


\alias{model.matrix.lmSubsets}
\alias{model.matrix.lmSelect}


\title{Extract a model matrix from a subset regression}


\description{

  Returns the model matrix for the specified submodel.

}


\usage{
\method{model.matrix}{lmSubsets}(object, size, best = 1, ...)

\method{model.matrix}{lmSelect}(object, best, ...)
}


\arguments{

  \item{object}{\code{"lmSubsets"}, \code{"lmSelect"}---a subset
    regression}

  \item{size}{\code{integer}---the submodel size}

  \item{best}{\code{integer}---the submodel position}

  \item{...}{forwarded to \code{\link[=model.frame]{model.frame()}}}

}


\value{

  \code{double[,]}---the model matrix

}


\seealso{

  \itemize{

    \item{\code{\link[=lmSubsets]{lmSubsets()}} for all-subsets
      regression}

    \item{\code{\link[=lmSelect]{lmSelect()}} for best-subset
      regression}

    \item{\code{\link[=model.matrix]{model.matrix()}} for the S3
      generic}

  }

}
