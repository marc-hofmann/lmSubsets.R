\name{lmSelect.lmSubsets}


\alias{lmSelect.lmSubsets}


\title{Best-subset regression}


\description{

  Coerce an all-subsets regression.

}


\usage{
\method{lmSelect}{lmSubsets}(formula, penalty = "BIC", ...)
}


\arguments{

  \item{formula}{\code{"lmSubsets"}---an all-subsets regression}

  \item{penalty}{\code{double}, \code{character},
    \code{"function"}---penalty per model parameter}

  \item{...}{ignored}

}


\details{

  Computes a best-subset regression from an all-subsets regression.

}


\value{

  \code{\link[=lmSelect]{"lmSelect"}}---a best-subset regression

}


\seealso{

  \itemize{

    \item{\code{\link[=lmSelect]{lmSelect()}} for the S3 generic}

    \item{\code{\link[=lmSubsets]{lmSubsets()}} for all-subsets
      regression}

  }

}


\examples{
data("AirPollution", package = "lmSubsets")

lm_all <- lmSubsets(mortality ~ ., data = AirPollution, nbest = 20)

lm_best <- lmSelect(lm_all)
lm_best
}


\keyword{regression}
