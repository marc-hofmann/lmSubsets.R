\name{lmSubsets.matrix}


\alias{lmSubsets.matrix}


\title{All-subsets regression}


\description{

  Matrix interface to all-variable-subsets selection in ordinary linear
  regression.

}


\usage{
\method{lmSubsets}{matrix}(formula, y, intercept = TRUE, ...)
}


\arguments{

  \item{formula}{\code{"matrix"}---the model matrix}

  \item{y}{\code{double[]}---the model response}

  \item{intercept}{\code{logical}---if \code{FALSE}, remove intercept
    term}

  \item{...}{forwarded to
    \code{\link[=lmSubsets.default]{lmSubsets.default()}}}

}


\details{

  This is a utility interface.  Use the standard formula interface
  wherever possible.

}


\value{

  \code{\link[=lmSubsets]{"lmSubsets"}}---an all-subsets regression

}


\seealso{

  \itemize{

    \item{\code{\link[=lmSubsets]{lmSubsets()}} for the S3 generic}

    \item{\code{\link[=lmSubsets.default]{lmSubsets.default()}} for the
      standard formula interface}

    }

}


\examples{
data("AirPollution", package = "lmSubsets")

x <- as.matrix(AirPollution)

lm_mat <- lmSubsets(x, y = "mortality")
lm_mat
}


\keyword{regression}
