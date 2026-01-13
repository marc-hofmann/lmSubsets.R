\name{lmSelect_fit}


\alias{lmSelect_fit}


\title{Best-subset regression}


\description{

  Low-level interface to best-variable-subset selection in ordinary
  linear regression.

}


\usage{
lmSelect_fit(x, y, weights = NULL, offset = NULL, include = NULL,
             exclude = NULL, penalty = "BIC", tolerance = 0,
             nbest = 1, ..., pradius = NULL)
}


\arguments{

  \item{x}{\code{double[,]}---the model matrix}

  \item{y}{\code{double[]}---the model response}

  \item{weights}{\code{double[]}---the model weights}

  \item{offset}{\code{double[]}---the model offset}

  \item{include}{\code{logical[]}, \code{integer[]},
    \code{character[]}---the regressors to force in}

  \item{exclude}{\code{logical[]}, \code{integer[]},
    \code{character[]}---the regressors to force out}

  \item{penalty}{\code{double}, \code{character},
    \code{"function"}---the penalty per model parameter}

  \item{tolerance}{\code{double}---the approximation tolerance}

  \item{nbest}{\code{integer}---the number of best subsets}

  \item{...}{ignored}

  \item{pradius}{\code{integer}---the preordering radius}

}


\details{

  The best variable-subset model is determined, where the "best" model
  is the one with the lowest information criterion value.  The
  information criterion belongs to the \link[=AIC]{AIC} family.

  The regression data is specified with the \code{x}, \code{y},
  \code{weights}, and \code{offset} parameters.  See
  \code{\link[=lm.fit]{lm.fit()}} for further details.

  To force regressors into or out of the regression, a list of
  regressors can be passed as an argument to the \code{include} or
  \code{exclude} parameters, respectively.

  The information criterion is specified with the \code{penalty}
  parameter.  Accepted values are \code{"AIC"}, \code{"BIC"}, or a
  \code{"numeric"} value representing the penalty-per-model-parameter.
  A custom selection criterion may be specified by passing an \code{R}
  function as an argument.  The expected signature is \code{function
  (size, rss)}, where \code{size} is the number of predictors (including
  the intercept, if any), and \code{rss} is the residual sum of squares.
  The function must be non-decreasing in both parameters.

  An approximation \code{tolerance} can be specified to speed up the
  search.

  The number of returned submodels is determined by the \code{nbest}
  parameter.

  The preordering radius is given with the \code{pradius} parameter.

}


\value{

  A \code{list} with the following components:

  \item{NOBS}{\code{integer}---number of observations in model (before
    \code{weights} processing)}

  \item{nobs}{\code{integer}---number of observations in model (after
    \code{weights} processing)}

  \item{nvar}{\code{integer}---number of regressors in model}

  \item{weights}{\code{double[]}---model weights}

  \item{intercept}{\code{logical}---is \code{TRUE} if model contains an
    intercept term, \code{FALSE} otherwise}

  \item{include}{\code{logical[]}---regressors forced into the
    regression}

  \item{exclude}{\code{logical[]}---regressors forced out of the
    regression}

  \item{size}{\code{integer[]}---subset sizes}

  \item{ic}{information criterion}

  \item{tolerance}{\code{double}---approximation tolerance}

  \item{nbest}{\code{integer}---number of best subsets}

  \item{submodel}{\code{"data.frame"}---submodel information}

  \item{subset}{\code{"data.frame"}---selected subsets}

}


\references{

  Hofmann M, Gatu C, Kontoghiorghes EJ, Colubi A, Zeileis A (2020).
  lmSubsets: Exact variable-subset selection in linear regression for R.
  \emph{Journal of Statistical Software}, \strong{93}, 1--21.
  \doi{10.18637/jss.v093.i03}.

}


\seealso{

  \itemize{

    \item{\code{\link[=lmSelect]{lmSelect()}} for the high-level
      interface}

    \item{\code{\link[=lmSubsets_fit]{lmSubsets_fit()}} for all-subsets
      regression}

  }

}


\examples{
data("AirPollution", package = "lmSubsets")

x <- as.matrix(AirPollution[, names(AirPollution) != "mortality"])
y <-           AirPollution[, names(AirPollution) == "mortality"]

f <- lmSelect_fit(x, y)
f
}


\keyword{regression}
