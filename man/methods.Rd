\name{methods}


\alias{methods}

\alias{model_response}

\alias{variable.names.lmSubsets}
\alias{formula.lmSubsets}
\alias{model.frame.lmSubsets}
\alias{model.matrix.lmSubsets}
\alias{model_response.lmSubsets}
\alias{refit.lmSubsets}
\alias{deviance.lmSubsets}
\alias{sigma.lmSubsets}
\alias{logLik.lmSubsets}
\alias{AIC.lmSubsets}
\alias{BIC.lmSubsets}
\alias{coef.lmSubsets}
\alias{vcov.lmSubsets}
\alias{fitted.lmSubsets}
\alias{residuals.lmSubsets}

\alias{variable.names.lmSelect}
\alias{formula.lmSelect}
\alias{model.frame.lmSelect}
\alias{model.matrix.lmSelect}
\alias{model_response.lmSelect}
\alias{refit.lmSelect}
\alias{deviance.lmSelect}
\alias{sigma.lmSelect}
\alias{logLik.lmSelect}
\alias{AIC.lmSelect}
\alias{BIC.lmSelect}
\alias{coef.lmSelect}
\alias{vcov.lmSelect}
\alias{fitted.lmSelect}
\alias{residuals.lmSelect}


\title{Methods for 'lmSubsets' and 'lmSelect' Objects}


\description{
  Extractor methods for \code{"lmSubsets"} and \code{"lmSelect"}
  objects.
}


\usage{
\method{variable.names}{lmSubsets}(object, size, best = 1, \ldots, na.rm = TRUE, drop = TRUE)
\method{variable.names}{lmSelect}(object, best = 1, \ldots, na.rm = TRUE, drop = TRUE)

\method{formula}{lmSubsets}(x, size, best = 1, \ldots)
\method{formula}{lmSelect}(x, best, \ldots)

\method{model.frame}{lmSubsets}(formula, \ldots)
\method{model.frame}{lmSelect}(formula, \ldots)

\method{model.matrix}{lmSubsets}(object, size, best = 1, \ldots)
\method{model.matrix}{lmSelect}(object, best, \ldots)

\method{model_response}{lmSubsets}(data, \ldots)
\method{model_response}{lmSelect}(data, \ldots)

\method{refit}{lmSubsets}(object, size, best = 1, \ldots)
\method{refit}{lmSelect}(object, best = 1, \ldots)

\method{deviance}{lmSubsets}(object, size, best = 1, \ldots, na.rm = TRUE, drop = TRUE)
\method{deviance}{lmSelect}(object, best = 1, \ldots, na.rm = TRUE, drop = TRUE)

\method{sigma}{lmSubsets}(object, size, best = 1, \ldots, na.rm = TRUE, drop = TRUE)
\method{sigma}{lmSelect}(object, best = 1, \ldots, na.rm = TRUE, drop = TRUE)

\method{logLik}{lmSubsets}(object, size, best = 1, \ldots, na.rm = TRUE, drop = TRUE)
\method{logLik}{lmSelect}(object, best = 1, \ldots, na.rm = TRUE, drop = TRUE)

\method{AIC}{lmSubsets}(object, size, best = 1, \ldots, k = 2, na.rm = TRUE, drop = TRUE)
\method{AIC}{lmSelect}(object, best = 1, \ldots, k = 2, na.rm = TRUE, drop = TRUE)

\method{BIC}{lmSubsets}(object, size, best = 1, \ldots, na.rm = TRUE, drop = TRUE)
\method{BIC}{lmSelect}(object, best = 1, \ldots, na.rm = TRUE, drop = TRUE)

\method{coef}{lmSubsets}(object, size, best = 1, \ldots, na.rm = TRUE, drop = TRUE)
\method{coef}{lmSelect}(object, best = 1, \ldots, na.rm = TRUE, drop = TRUE)

\method{vcov}{lmSubsets}(object, size, best = 1, \ldots)
\method{vcov}{lmSelect}(object, best = 1, \ldots)

\method{fitted}{lmSubsets}(object, size, best = 1, \ldots)
\method{fitted}{lmSelect}(object, best = 1, \ldots)

\method{residuals}{lmSubsets}(object, size, best = 1, \ldots)
\method{residuals}{lmSelect}(object, best = 1, \ldots)
}


\arguments{
  \item{object, formula, data, x}{An object of class \code{"lmSubsets"} or
    \code{"lmSelect"}.}
  \item{size}{The submodel size.}
  \item{best}{The submodel ranking.}
  \item{\dots}{Forwarded arguments.}
  \item{k}{AIC penalty.}
  \item{drop}{Control shape of return value.}
  \item{na.rm}{Remove missing submodels.}
}


\details{
  The extractor methods work for \code{"lmSubsets"} and
  \code{"lmSelect"} objects, i.e., objects that have been generated
  using the formula interface.

  If \code{drop == FALSE}, the extractor methods \code{variable.names},
  \code{deviance}, \code{sigma}, \code{logLik}, \code{AIC}, \code{BIC}
  and \code{coef} return a \code{data.frame} object.  If \code{drop ==
  TRUE}, the return value is a \code{logical} matrix
  (\code{variable.names}), a \code{numeric} matrix (\code{coef}), or a
  \code{numeric} vector.  If the \code{drop} parameter is not set
  explicitly when calling \code{variable.names} or \code{coef},
  one-dimensional values are represented in a compact form.

  If a desired extractor function is not available, \code{refit} can be
  called explicitly to obtain the corresponding \code{"lm"} object.
}


\seealso{\code{\link{lmSubsets}}, \code{\link{lmSelect}},
  \code{\link{refit}}.}


\examples{
## load data
data("AirPollution", package = "lmSubsets")

## fit subsets (5 best subsets per size)
lm_all <- lmSubsets(mortality ~ ., data = AirPollution, nbest = 5)

## extract information (for best submodel of size 3)
coef(lm_all, size = 3)
vcov(lm_all, size = 3)
residuals(lm_all, size = 3)
fitted(lm_all, size = 3)
model.matrix(lm_all, size = 3)

## select best (BIC) submodels
lm_best <- lmSelect(lm_all)

## extract information
deviance(lm_best)
logLik(lm_best)
AIC(lm_best)
BIC(lm_best, best = 1:5)

## refit model
lm5 <- refit(lm_all, size = 5)
summary(lm5)
## (Note that the p-values are not valid due to model selection.)
}


\keyword{regression}
