\name{refit.lmSubsets}


\alias{refit.lmSubsets}
\alias{refit.lmSelect}


\title{Refit a subset regression}


\description{

  Fit the specified submodel and return the obtained \code{"lm"} object.

}


\usage{
\method{refit}{lmSubsets}(object, size, best = 1, ...)

\method{refit}{lmSelect}(object, best = 1, ...)
}


\arguments{

  \item{object}{\code{"lmSubsets"}, \code{"lmSelect"}---a subset
    regression}

  \item{size}{\code{integer}---the submodel size}

  \item{best}{\code{integer}---the submodel position}

  \item{...}{ignored}

}


\value{

  \code{"lm"}---the fitted model

}


\seealso{

  \itemize{

    \item{\code{\link[=lmSubsets]{lmSubsets()}} for all-subsets
      regression}

    \item{\code{\link[=lmSelect]{lmSelect()}} for best-subset
      regression}

    \item{\code{\link[=refit]{refit()}} for the S3 generic}

  }

}


\examples{
## load data
data("AirPollution", package = "lmSubsets")

## fit subsets
lm_all <- lmSubsets(mortality ~ ., data = AirPollution)

## refit best model
lm5 <- refit(lm_all, size = 5)
summary(lm5)
}
