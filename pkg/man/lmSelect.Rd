\name{lmSelect}


\alias{lmSelect}
\alias{lmSelect.default}
\alias{print.lmSelect}


\title{Best-subset regression}


\description{

  Best-variable-subset selection in ordinary linear regression.

}


\usage{
lmSelect(formula, ...)

\method{lmSelect}{default}(formula, data, subset, weights, na.action,
         model = TRUE, x = FALSE, y = FALSE, contrasts = NULL,
         offset, ...)
}


\arguments{

  \item{formula, data, subset, weights, na.action, model, x, y,
    contrasts, offset}{standard formula interface}

  \item{...}{forwarded to \code{lmSelect_fit()}}

}


\details{

  The \code{lmSelect()} generic provides various methods to conveniently
  specify the regressor and response variables.  The standard formula
  interface (see \code{\link[=lm]{lm()}}) can be used, or the model
  information can be extracted from an already fitted \code{"lm"}
  object.  The model matrix and response can also be passed in directly.

  After processing the arguments, the call is forwarded to
  \code{\link[=lmSelect_fit]{lmSelect_fit()}}.

}


\value{

  \code{"lmSelect"}---a \code{list} containing the components returned
  by \code{lmSelect_fit()}

  Further components include \code{call}, \code{na.action},
  \code{weights}, \code{offset}, \code{contrasts}, \code{xlevels},
  \code{terms}, \code{mf}, \code{x}, and \code{y}.  See
  \code{\link[=lm]{lm()}} for more information.

}


\seealso{

  \itemize{

    \item{\code{\link[=lmSelect.matrix]{lmSelect.matrix()}} for the
      matrix interface}

    \item{\code{\link[=lmSelect.lmSubsets]{lmSelect.lmSubsets()}} for
      coercing an all-subsets regression}

    \item{\code{\link[=lmSelect_fit]{lmSelect_fit()}} for the low-level
      interface}

    \item{\code{\link[=lmSubsets]{lmSubsets()}} for all-subsets
      regression}

  }

}


\examples{
## load data
data("AirPollution", package = "lmSubsets")


###################
##  basic usage  ##
###################

## fit 20 best subsets (BIC)
lm_best <- lmSelect(mortality ~ ., data = AirPollution, nbest = 20)
lm_best

## summary statistics
summary(lm_best)

## visualize
plot(lm_best)


########################
##  custom criterion  ##
########################

## the same as above, but with a custom criterion:
M <- nrow(AirPollution)

ll <- function (rss) {
  -M/2 * (log(2 * pi) - log(M) + log(rss) + 1)
}

aic <- function (size, rss, k = 2) {
  -2 * ll(rss) + k * (size + 1)
}

bic <- function (size, rss) {
  aic(size, rss, k = log(M))
}

lm_cust <- lmSelect(mortality ~ ., data = AirPollution,
                    penalty = bic, nbest = 20)
lm_cust
}


\keyword{regression}
