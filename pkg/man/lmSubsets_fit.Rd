\name{lmSubsets_fit}


\alias{lmSubsets_fit}


\title{All-subsets regression}


\description{

  Low-level interface to all-variable-subsets selection in ordinary
  linear regression.

}


\usage{
lmSubsets_fit(x, y, weights = NULL, offset = NULL, include = NULL,
              exclude = NULL, nmin = NULL, nmax = NULL,
              tolerance = 0, nbest = 1, ..., pradius = NULL)
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

  \item{nmin}{\code{integer}---the minimum number of regressors}

  \item{nmax}{\code{integer}---the maximum number of regressors}

  \item{tolerance}{\code{double[]}---the approximation tolerances}

  \item{nbest}{\code{integer}---the number of best subsets}

  \item{...}{ignored}

  \item{pradius}{\code{integer}---the preordering radius}

}


\details{

  The best variable-subset model for every subset size is determined,
  where the "best" model is the one with the lowest residual sum of
  squares (RSS).

  The regression data is specified with the \code{x}, \code{y},
  \code{weights}, and \code{offset} parameters.  See
  \code{\link[=lm.fit]{lm.fit()}} for further details.

  To force regressors into or out of the regression, a list of
  regressors can be passed as an argument to the \code{include} or
  \code{exclude} parameters, respectively.

  The scope of the search can be limited to a range of subset sizes by
  setting \code{nmin} and \code{nmax}, the minimum and maximum number of
  regressors allowed in the regression, respectively.

  A \code{tolerance} vector can be specified to speed up the search,
  where \code{tolerance[j]} is the approximation tolerance applied to
  subset models of size \code{j}.

  The number of submodels returned for each subset size is determined by
  the \code{nbest} parameter.

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

  \item{tolerance}{\code{double[]}---approximation tolerances}

  \item{nbest}{\code{integer}---number of best subsets}

  \item{submodel}{\code{"data.frame"}---submodel information}

  \item{subset}{\code{"data.frame"}---variable subsets}

}


\references{

  Hofmann M, Gatu C, Kontoghiorghes EJ, Colubi A, Zeileis A (2020).
  lmSubsets: Exact variable-subset selection in linear regression for R.
  \emph{Journal of Statistical Software}, \strong{93}, 1--21.
  \doi{10.18637/jss.v093.i03}.

}


\seealso{

  \itemize{

    \item{\code{\link[=lmSubsets]{lmSubsets()}} for the high-level
      interface}

    \item{\code{\link[=lmSelect_fit]{lmSelect_fit()}} for best-subset
      regression}

  }

}


\examples{
data("AirPollution", package = "lmSubsets")

x <- as.matrix(AirPollution[, names(AirPollution) != "mortality"])
y <-           AirPollution[, names(AirPollution) == "mortality"]

f <- lmSubsets_fit(x, y)
f
}


\keyword{regression}
