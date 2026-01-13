\name{AirPollution}


\docType{data}


\alias{AirPollution}


\title{Air pollution and mortality}


\description{

  Data relating air pollution and mortality, frequently used for
  illustrations in ridge regression and related tasks.

}


\usage{
data(AirPollution)
}


\format{

  A data frame containing 60 observations on 16 variables.

  \describe{

    \item{precipitation}{average annual precipitation in inches}

    \item{temperature1}{average January temperature in degrees
      Fahrenheit}

    \item{temperature7}{average July temperature in degrees
      Fahrenheit}

    \item{age}{percentage of 1960 SMSA population aged 65 or older}

    \item{household}{average household size}

    \item{education}{median school years completed by those over 22}

    \item{housing}{percentage of housing units which are sound and
      with all facilities}

    \item{population}{population per square mile in urbanized areas,
      1960}

    \item{noncauc}{percentage of non-Caucasian population in
      urbanized areas, 1960}

    \item{whitecollar}{percentage employed in white collar
      occupations}

    \item{income}{percentage of families with income < USD 3000}

    \item{hydrocarbon}{relative hydrocarbon pollution potential}

    \item{nox}{relative nitric oxides potential}

    \item{so2}{relative sulphur dioxide potential}

    \item{humidity}{annual average percentage of relative humidity at
      13:00}

    \item{mortality}{total age-adjusted mortality rate per 100,000}

  }

}


\source{

  \url{https://lib.stat.cmu.edu/datasets/pollution}

}


\references{

  McDonald GC, Schwing RC (1973).  Instabilities of regression estimates
  relating air pollution to mortality.  \emph{Technometrics},
  \strong{15}, 463--482.

  Miller AJ (2002). \emph{Subset selection in regression}.  New York:
  Chapman and Hall.

}


\examples{
## load data (with logs for relative potentials)
data("AirPollution", package = "lmSubsets")
for (i in 12:14)  AirPollution[[i]] <- log(AirPollution[[i]])

## fit subsets
lm_all <- lmSubsets(mortality ~ ., data = AirPollution)
plot(lm_all)

## refit best model
lm6 <- refit(lm_all, size = 6)
summary(lm6)
}


\keyword{datasets}
