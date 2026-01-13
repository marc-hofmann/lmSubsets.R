\name{model.frame.lmSubsets}

\alias{model.frame.lmSubsets}
\alias{model.frame.lmSelect}


\title{Extract the model frame from a subset regression}


\description{
  Return the model frame.
}


\usage{
\method{model.frame}{lmSubsets}(formula, ...)

\method{model.frame}{lmSelect}(formula, ...)
}


\arguments{

  \item{formula}{\code{"lmSubsets"}, \code{"lmSelect"}---a subset
    regression}

  \item{...}{forwarded to \code{model.frame()}}

}


\value{

  \code{"data.frame"}---the model frame

}


\seealso{

  \itemize{

    \item{\code{\link[=lmSubsets]{lmSubsets()}} for all-subsets
      regression}

    \item{\code{\link[=lmSelect]{lmSelect()}} for best-subset
      regression}

    \item{\code{\link[=model.frame]{model.frame()}} for the S3 generic}

  }

}
