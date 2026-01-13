\name{image}


\alias{image}
\alias{image.lmSubsets}
\alias{image.lmSelect}


\title{Variable Selection Heatmaps}


\description{Visualization of variable subsets.}



\usage{

\method{image}{lmSubsets}(x, size = NULL, best = 1, which = NULL,
      hilite, hilite_penalty, main, sub, xlab = NULL, ylab,
      ann = par("ann"), axes = TRUE, col = c("gray40", "gray90"),
      lab = "lab", col_hilite = cbind("red", "pink"),
      lab_hilite = "lab", pad_size = 3, pad_best = 1,
      pad_which = 3, axis_pos = -4, axis_tck = -4,
      axis_lab = -10, ...)

\method{image}{lmSelect}(x, best = NULL, which = NULL, hilite,
      hilite_penalty, main, sub = NULL, xlab = NULL, ylab,
      ann = par("ann"), axes = TRUE, col = c("gray40", "gray90"),
      lab = "lab",  col_hilite = cbind("red", "pink"),
      lab_hilite = "lab", pad_best = 2, pad_which = 2,
      axis_pos = -4, axis_tck = -4, axis_lab = -10, ...)
}


\arguments{
  \item{x}{An object of class \code{"lmSubsets"} or \code{"lmSelect"}.}
  \item{main, sub, xlab, ylab}{Main, sub- and axis titles.}
  \item{size, best}{Submodels to be plotted.}
  \item{which}{Regressors to be plotted.}
  \item{col, lab}{Color and label style.}
  \item{hilite, hilite_penalty}{Submodels to be highlighted.}
  \item{col_hilite, lab_hilite}{Highlighting style.}
  \item{pad_size, pad_best, pad_which}{Padding.}
  \item{axis_pos, axis_tck, axis_lab}{Position of axes, tick length, and
    position of labels}
  \item{\dots}{Ignored.}
  \item{axes}{Plot axes.}
  \item{ann}{Annotate plot.}
}


\seealso{\code{\link{lmSubsets}}, \code{\link{lmSelect}}.}


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
