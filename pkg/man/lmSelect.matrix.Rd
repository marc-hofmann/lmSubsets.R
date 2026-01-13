\name{lmSelect.matrix}


\alias{lmSelect.matrix}


\title{Best-subset regression}


\description{

  Matrix interface to best-variable-subset selection in ordinary linear
  regression.

}


\usage{
\method{lmSelect}{matrix}(formula, y, intercept = TRUE, ...)
}


\arguments{

  \item{formula}{\code{"matrix"}---the model matrix}

  \item{y}{\code{double[]}---the model response}

  \item{intercept}{\code{logical[]}---if \code{FALSE}, remove intercept
    term}

  \item{...}{forwarded to
    \code{\link[=lmSelect.default]{lmSelect.default()}}}

}


\details{

  This is a utility interface.  Use the standard formula interface
  wherever possible.

}


\value{

  \code{\link[=lmSelect]{"lmSelect"}}---a best-subset regression

}


\seealso{

  \itemize{

    \item{\code{\link[=lmSelect]{lmSelect()}} for the S3 generic}

    \item{\code{\link[=lmSelect.default]{lmSelect.default()}} for the
      standard formula interface}

  }

}


\keyword{regression}
