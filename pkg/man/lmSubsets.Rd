\name{lmSubsets}


\alias{lmSubsets}
\alias{lmSubsets.default}
\alias{print.lmSubsets}


\title{All-subsets regression}


\description{

  All-variable-subsets selection in ordinary linear regression.

}


\usage{
lmSubsets(formula, ...)

\method{lmSubsets}{default}(formula, data, subset, weights, na.action,
          model = TRUE, x = FALSE, y = FALSE, contrasts = NULL,
          offset, ...)
}


\arguments{

  \item{formula, data, subset, weights, na.action, model, x, y,
    contrasts, offset}{standard formula interface}

  \item{...}{fowarded to \code{lmSubsets_fit()}}

}


\details{

  The \code{lmSubsets()} generic provides various methods to
  conveniently specify the regressor and response variables.  The
  standard formula interface (see \code{\link[=lm]{lm()}}) can be used,
  or the model information can be extracted from an already fitted
  \code{"lm"} object.  The model matrix and response can also be passed
  in directly.

  After processing of the arguments, the call is forwarded to
  \code{\link[=lmSubsets_fit]{lmSubsets_fit()}}.

}


\value{

  \code{"lmSubsets"}---a \code{list} containing the components returned
  by \code{lmSubsets_fit()}

  Further components include \code{call}, \code{na.action},
  \code{weights}, \code{offset}, \code{contrasts}, \code{xlevels},
  \code{terms}, \code{mf}, \code{x}, and \code{y}.  See
  \code{\link[=lm]{lm()}} for more information.

}


\seealso{

  \itemize{

    \item{\code{\link[=lmSubsets.matrix]{lmSubsets.matrix()}} for the
      \code{"matrix"} interface}

    \item{\code{\link[=lmSubsets_fit]{lmSubsets_fit()}} for the
      low-level interface}

    \item{\code{\link[=lmSelect]{lmSelect()}} for best-subset
      regression}

  }

}


\examples{
## load data
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
}


\keyword{regression}
