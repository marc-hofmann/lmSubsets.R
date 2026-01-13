\name{lmSelect}


\alias{lmSelect}
\alias{lmSelect.default}
\alias{lmSelect.matrix}
\alias{lmSelect.lmSubsets}
\alias{lmSelect_fit}
\alias{lmSubsets_select}
\alias{print.lmSelect}
\alias{plot.lmSelect}
\alias{summary.lmSelect}
\alias{print.summary.lmSelect}


\title{Best-Subset Regression}


\description{Best-subset regression for ordinary linear models.}

\usage{
lmSelect(formula, \ldots)

\method{lmSelect}{default}(formula, data, subset, weights, na.action, model = TRUE,
         x = FALSE, y = FALSE, contrasts = NULL, offset, \ldots)

\method{lmSelect}{matrix}(formula, y, intercept = TRUE, \ldots)

\method{lmSelect}{lmSubsets}(formula, penalty = "BIC", \ldots)

lmSelect_fit(x, y, weights = NULL, offset = NULL, include = NULL,
             exclude = NULL, penalty = "BIC", tolerance = 0,
             nbest = 1, \ldots, pradius = NULL)
}


\arguments{
  \item{formula, data, subset, weights, na.action, model, x, y, contrasts,
    offset}{Standard formula interface.}
  \item{intercept}{Include intercept.}
  \item{include, exclude}{Force regressors in or out.}
  \item{penalty}{Penalty per parameter.}
  \item{tolerance}{Approximation tolerance.}
  \item{nbest}{Number of best subsets.}
  \item{\dots}{Forwarded to \code{lmSelect_fit}.}
  \item{pradius}{Preordering radius.}
}


\details{
  The \code{lmSelect} generic provides a convenient interface for best
  variable-subset selection in linear regression: The \code{nbest} best
  -- according to an information criterion of the AIC family -- subset
  models are returned.

  The information criterion is specified with the \code{penalty}
  parameter.  Accepted values are \code{"AIC"}, \code{"BIC"}, or a
  \code{numeric} value representing the penalty per model parameter (see
  \code{\link[stats]{AIC}}).

  A custom selection criterion may be specified by passing an R function
  as the \code{penalty} argument.  The expected signature is
  \code{function(size, rss)}, where \code{size} is the number of
  predictors (including intercept, if any), and \code{rss} the residual
  sum of squares.  The function must be non-decreasing in both
  parameters.

  A low-level matrix interface is provided by \code{lmSelect_fit}.

  See \code{\link{lmSubsets}} for further information.
}


\value{
  An object of class \code{"lmSelect"}, i.e., a list with the following
  components:
  \item{nobs, nvar}{Number of observations, of variables.}
  \item{intercept}{\code{TRUE} if model has intercept term;
    \code{FALSE} otherwise.}
  \item{include, exclude}{Included, excluded variables.}
  \item{size}{Subset sizes.}
  \item{tolerance}{Approximation tolerance.}
  \item{nbest}{Number of best subsets.}
  \item{submodel}{Submodel information.}
  \item{subset}{Selected variables.}

  Further components include \code{call}, \code{na.action},
  \code{weights}, \code{offset}, \code{contrasts}, \code{xlevels},
  \code{terms}, \code{mf}, \code{x}, and \code{y}.  See
  \code{\link[stats]{lm}} for more information.
}


\references{
  Hofmann M, Gatu C, Kontoghiorghes EJ, Colubi A, Zeileis A (2018).
  lmSubsets: Exact Variable-Subset Selection in Linear Regression for R.
  \emph{Journal of Statistical Software}.  Submitted.
}


\seealso{\code{\link{lmSubsets}}, \code{\link{summary}},
  \link{methods}.}


\examples{
## load data (with logs for relative potentials)
data("AirPollution", package = "lmSubsets")


###################
##  basic usage  ##
###################

## fit 20 best subsets (BIC)
lm_best <- lmSelect(mortality ~ ., data = AirPollution, nbest = 20)
lm_best

## equivalent to:
\dontrun{
lm_all <- lmSubsets(mortality ~ ., data = AirPollution, nbest = 20)
lm_best <- lmSelect(lm_all)
}

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
