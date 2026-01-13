\name{plot.lmSubsets}


\alias{plot.lmSubsets}
\alias{plot.lmSelect}


\title{Plot a subset regression}


\description{

  Plot the deviance of the selected submodels, as well as a specified
  information criterion.

}


\usage{
\method{plot}{lmSubsets}(x, penalty = "BIC", xlim, ylim_rss, ylim_ic, type_rss = "o",
     type_ic = "o", main, sub, xlab, ylab_rss, ylab_ic, legend_rss,
     legend_ic, ann = par("ann"), axes = TRUE, lty_rss = c(1, 3),
     pch_rss = c(16, 21), col_rss = "black", bg_rss = "white",
     lty_ic = c(1, 3), pch_ic = c(16, 21), col_ic = "red",
     bg_ic = "white", ...)

\method{plot}{lmSelect}(x, xlim, ylim, type = "o", main, sub, xlab, ylab, legend,
     ann = par("ann"), axes = TRUE, lty = 1, pch = 16, col = "red",
     bg = "white", ...)
}


\arguments{

  \item{x}{\code{"lmSubsets"}, \code{"lmSelect"}---a subset regression}

  \item{penalty}{the information criterion}

  \item{xlim, ylim, ylim_rss, ylim_ic}{x and y limits}

  \item{type, type_rss, type_ic}{type of plot}

  \item{main, sub}{main and sub-title}

  \item{xlab, ylab, ylab_rss, ylab_ic}{axis titles}

  \item{legend, legend_rss, legend_ic}{plot legend}

  \item{ann}{annotate plot}

  \item{axes}{plot axes}

  \item{lty, lty_rss, lty_ic}{line type}

  \item{pch, pch_rss, pch_ic}{plotting character}

  \item{col, col_rss, col_ic}{color}

  \item{bg, bg_rss, bg_ic}{background color}

  \item{...}{further graphical parameters}

}


\value{

  \code{invisible(x)}

}



\seealso{

  \itemize{

    \item{\code{\link[=lmSubsets]{lmSubsets()}} for all-subsets
      regression}

    \item{\code{\link[=lmSelect]{lmSelect()}} for best-subset
      regression}

    \item{\code{\link[=plot]{plot()}} for the S3 generic}

  }

}


\examples{
## load data
data("AirPollution", package = "lmSubsets")


#################
##  lmSubsets  ##
#################

lm_all <- lmSubsets(mortality ~ ., data = AirPollution, nbest = 5)
plot(lm_all)


################
##  lmSelect  ##
################

lm_best <- lmSelect(mortality ~ ., data = AirPollution, nbest = 20)
plot(lm_best)
}


\keyword{regression}
