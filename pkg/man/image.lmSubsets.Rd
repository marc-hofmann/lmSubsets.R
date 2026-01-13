\name{image.lmSubsets}


\alias{image.lmSubsets}
\alias{image.lmSelect}


\title{Heatmap of a subset regression}


\description{

  Plot a heatmap of the specified submodels.

}


\usage{
\method{image}{lmSubsets}(x, size = NULL, best = 1, which = NULL, hilite, hilite_penalty,
      main, sub, xlab = NULL, ylab, ann = par("ann"), axes = TRUE,
      col = c("gray40", "gray90"), lab = "lab",
      col_hilite = cbind("red", "pink"), lab_hilite = "lab",
      pad_size = 3, pad_best = 1, pad_which = 3, axis_pos = -4,
      axis_tck = -4, axis_lab = -10, ...)

\method{image}{lmSelect}(x, best = NULL, which = NULL, hilite, hilite_penalty,
      main, sub = NULL, xlab = NULL, ylab, ann = par("ann"),
      axes = TRUE, col = c("gray40", "gray90"), lab = "lab",
      col_hilite = cbind("red", "pink"), lab_hilite = "lab",
      pad_best = 2, pad_which = 2, axis_pos = -4, axis_tck = -4,
      axis_lab = -10, ...)
}


\arguments{

  \item{x}{\code{"lmSubsets"}, \code{"lmSelect"}---a subset regression}

  \item{size, best}{submodels to be plotted}

  \item{which}{regressors to be plotted}

  \item{hilite, hilite_penalty}{submodels to be highlighted}

  \item{main, sub, xlab, ylab}{main, sub-, and axis titles}

  \item{ann}{annotate plot}

  \item{axes}{plot axes}

  \item{col, lab}{color and label style}

  \item{col_hilite, lab_hilite}{highlighting style}

  \item{pad_size, pad_best, pad_which}{padding}

  \item{axis_pos, axis_tck, axis_lab}{position of axes, tick length,
    and position of labels}

  \item{...}{ignored}

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

  }

}


\examples{
## data
data("AirPollution", package = "lmSubsets")


#################
##  lmSubsets  ##
#################

lm_all <- lmSubsets(mortality ~ ., data = AirPollution, nbest = 20)

## heatmap
image(lm_all, best = 1:3)

## highlight 5 best (BIC)
image(lm_all, best = 1:3, hilite = 1:5, hilite_penalty = "BIC")


################
##  lmSelect  ##
################

## default criterion: BIC
lm_best <- lmSelect(lm_all)

## highlight 5 best (AIC)
image(lm_best, hilite = 1:5, hilite_penalty = "AIC")

## axis labels
image(lm_best, lab = c("bold(lab)", "lab"), hilite = 1,
      lab_hilite = "underline(lab)")
}


\keyword{regression}