\name{lmSubsets}


\alias{lmSubsets}
\alias{lmSubsets.default}
\alias{lmSubsets.matrix}
\alias{lmSubsets_fit}
\alias{print.lmSubsets}
\alias{plot.lmSubsets}
\alias{summary.lmSubsets}
\alias{print.summary.lmSubsets}


\title{All-Subsets Regression}


\description{All-subsets regression for linear models estimated by
  ordinary least squares (OLS).}


\usage{
lmSubsets(formula, \ldots)

\method{lmSubsets}{default}(formula, data, subset, weights, na.action, model = TRUE,
          x = FALSE, y = FALSE, contrasts = NULL, offset, \ldots)

\method{lmSubsets}{matrix}(formula, y, intercept = TRUE, \ldots)

lmSubsets_fit(x, y, weights = NULL, offset = NULL, include = NULL,
              exclude = NULL, nmin = NULL, nmax = NULL,
              tolerance = 0, nbest = 1, \ldots, pradius = NULL)
}


\arguments{
  \item{formula, data, subset, weights, na.action, model, x, y, contrasts,
    offset}{Standard formula interface.}
  \item{intercept}{Include intercept.}
  \item{include, exclude}{Force regressors in or out.}
  \item{nmin, nmax}{Minimum and maximum number of regressors.}
  \item{tolerance}{Approximation tolerance (vector).}
  \item{nbest}{Number of best subsets.}
  \item{\dots}{Forwarded to \code{lmSubsets.default} and \code{lmSubsets_fit}.}
  \item{pradius}{Preordering radius.}
}


\details{
  The \code{lmSubsets} generic provides various methods to conveniently
  specify the regressor and response variables.  The standard
  \code{formula} interface (see \code{\link[stats]{lm}}) can be used, or
  the information can be extracted from an already fitted \code{"lm"}
  object.  The regressor matrix and response variable can also be passed
  in directly (see Examples).

  The call is forwarded to \code{lmSubsets_fit}, which provides a
  low-level matrix interface.

  The \code{nbest} best subset models for every subset size are
  computed, where the "best" models are the models with the lowest
  residual sum of squares (RSS).  The scope of the search can be limited
  to a range of subset sizes by setting \code{nmin} and \code{nmax}.  A
  tolerance vector (expanded if necessary) may be specified to speed up
  the search, where \code{tolerance[j]} is the tolerance applied to
  subset models of size \code{j}.

  By way of \code{include} and \code{exclude}, variables may be forced
  in to or out of the regression, respectively.

  The extent to which variables are preordered is controlled with the
  \code{pradius} parameter.

  A set of standard extractor functions for fitted model objects is
  available for objects of class \code{"lmSubsets"}.  See
  \code{\link{methods}} for more details.

  The \code{summary} method can be called to obtain summary statistics.
}


\value{
  An object of class \code{"lmSubsets"}, i.e., a list with the
  following components:
  \item{nobs, nvar}{Number of observations, of variables.}
  \item{intercept}{\code{TRUE} if model has intercept term;
    \code{FALSE} otherwise.}
  \item{include, exclude}{Included, excluded regressors.}
  \item{size}{Subset sizes.}
  \item{tolerance}{Approximation tolerance (vector).}
  \item{nbest}{Number of best subsets.}
  \item{submodel}{Submodel information.}
  \item{subset}{Selected variables.}

  Further components include \code{call}, \code{na.action},
  \code{weights}, \code{offset}, \code{contrasts}, \code{xlevels},
  \code{terms}, \code{mf}, \code{x}, and \code{y}.  See
  \code{\link[stats]{lm}} for more information.
}


\references{
  Hofmann M, Gatu C, Kontoghiorghes EJ, Colubi A, Zeileis A (2019).
  lmSubsets: Exact Variable-Subset Selection in Linear Regression for R.
  \emph{Journal of Statistical Software}.  In press.

  Hofmann M, Gatu C, Kontoghiorghes EJ (2007).  Efficient Algorithms for
  Computing the Best Subset Regression Models for Large-Scale Problems.
  \emph{Computational Statistics \& Data Analysis}, \bold{52}, 16--29.

  Gatu C, Kontoghiorghes EJ (2006).  Branch-and-Bound Algorithms for
  Computing the Best Subset Regression Models.  \emph{Journal of
    Computational and Graphical Statistics}, \bold{15}, 139--156.
}


\seealso{\code{\link{lmSelect}}, \code{\link{summary}},
  \code{\link{methods}}.}


\examples{
## load data (with logs for relative potentials)
data("AirPollution", package = "lmSubsets")


###################
##  basic usage  ##
###################

## canonical example: fit all subsets
lm_all <- lmSubsets(mortality ~ ., data = AirPollution, nbest = 5)
lm_all

## plot RSS and BIC
plot(lm_all)

## summary statistics
summary(lm_all)


############################
##  forced in-/exclusion  ##
############################

lm_force <- lmSubsets(lm_all, include = c("nox", "so2"),
                      exclude = "whitecollar")
lm_force


########################
##  matrix interface  ##
########################

## same as above
x <- as.matrix(AirPollution)

lm_mat <- lmSubsets(x, y = "mortality")
lm_mat
}


\keyword{regression}
