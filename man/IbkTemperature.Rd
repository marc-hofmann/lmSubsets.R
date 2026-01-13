\name{IbkTemperature}


\alias{IbkTemperature}


\title{Temperature Observations and Numerical Weather Predictions for Innsbruck}


\description{
  00UTC temperature observations and corresponding 24-hour
  reforecast ensemble means from the Global Ensemble Forecast System
  (GEFS, Hamill et al. 2013) for SYNOP station Innsbruck Airport (11120;
  47.260, 11.357) from 2011-01-01 to 2015-12-31.
}


\usage{data("IbkTemperature")}


\format{
  A data frame containing 1824 daily observations/forecasts for 42
  variables.  The first column (\code{temp}) contains temperature
  observations at 00UTC (coordinated universal time), columns 2--37 are
  24-hour lead time GEFS reforecast ensemble means for different
  variables (see below).  Columns 38--42 are deterministic time
  trend/season patterns.

  \describe{
    \item{temp}{Observed temperature at Innsbruck Airport (deg \eqn{C}).} 
    \item{tp}{Total accumulated precipitation (\eqn{kg~m^{-2}}{kg m^-2}).} 
    \item{t2m}{Temperature at 2 meters (\eqn{K}).}
    \item{u10m}{U-component of wind at 10 meters (\eqn{m~s^{-1}}{m s^-1}).} 
    \item{v10m}{V-component of wind at 10 meters (\eqn{m~s^{-1}}{m s^-1}).} 
    \item{u80m}{U-component of wind at 80 meters (\eqn{m~s^{-1}}{m s^-1}).} 
    \item{v80m}{U-component of wind at 80 meters (\eqn{m~s^{-1}}{m s^-1}).} 
    \item{cape}{Convective available potential energy (\eqn{J~kg^{-1}}{J kg^-1}).} 
    \item{ci}{Convective inhibition (\eqn{J~kg^{-1}}{J kg^-1}).} 
    \item{sdlwrf}{Surface downward long-wave radiation flux (\eqn{W~m^{-2}}{W m^-2}).} 
    \item{sdswrf}{Surface downward short-wave radiation flux (\eqn{W~m^{-2}}{W m^-2}).} 
    \item{sulwrf}{Surface upward long-wave radiation flux (\eqn{W~m^{-2}}{W m^-2}).} 
    \item{suswrf}{Surface upward short-wave radiation flux (\eqn{W~m^{-2}}{W m^-2}).}
    \item{ghf}{Ground heat flux (\eqn{W~m^{-2}}{W m^-2}).} 
    \item{slhnf}{Surface latent heat net flux (\eqn{W~m^{-2}}{W m^-2}).} 
    \item{sshnf}{Surface sensible heat net flux (\eqn{W~m^{-2}}{W m^-2}).}
    \item{mslp}{Mean sea level pressure (\eqn{Pa}).} 
    \item{psfc}{Surface pressure (\eqn{Pa}).} 
    \item{pw}{Precipitable water (\eqn{kg~m^{-2}}{kg m^-2}).}
    \item{vsmc}{Volumetric soil moisture content (fraction).} 
    \item{sh2m}{Specific humidity at 2 meters (\eqn{kg~kg^{-1}}{kg kg^-1}).} 
    \item{tcc}{Total cloud cover (percent).} 
    \item{tcic}{Total column-integrated condensate (\eqn{kg~m^{-2}}{kg m^-2}).} 
    \item{tsfc}{Skin temperature (\eqn{K}).} 
    \item{tmax2m}{Maximum temperature (\eqn{K}).} 
    \item{tmin2m}{Minimum temperature (\eqn{K}).} 
    \item{st}{Soil temperature (0-10 cm below surface) (\eqn{K}).} 
    \item{ulwrf}{Upward long-wave radiation flux (\eqn{W~m^{-2}}{W m^-2}).}  
    \item{wr}{Water runoff (\eqn{kg~m^{-2}}{kg m^-2}).} 
    \item{we}{Water equivalent of accumulated snow depth (\eqn{kg~m^{-2}}{kg m^-2}).} 
    \item{wp}{Wind mixing energy (\eqn{J}).} 
    \item{w850}{Vertical velocity at 850 hPa surface (\eqn{Pa~s^{-1}}{Pa s^-1}).} 
    \item{t2pvu}{Temperature on 2 PVU surface (\eqn{K}).} 
    \item{p2pvu}{Pressure on 2 PVU surface (\eqn{Pa}).} 
    \item{u2pvu}{U-component of wind on 2 PVU surface (\eqn{m~s^{-1}}{m s^-1}).} 
    \item{v2pvu}{U-component of wind on 2 PVU surface (\eqn{m~s^{-1}}{m s^-1}).} 
    \item{pv}{Potential vorticity on 320 K isentrope (\eqn{K~m^{2}~kg^{-1}~s^{-1}}{K m^2 kg^-1 s^-1}).} 
    \item{time}{Time in years.}
    \item{sin, cos}{Sine and cosine component of annual harmonic pattern.}
    \item{sin2, cos2}{Sine and cosine component of bi-annual harmonic pattern.}
  }
}


\source{
  Observations: \url{http://www.ogimet.com/synops.phtml.en}

  Reforecasts: \url{http://www.esrl.noaa.gov/psd/forecasts/reforecast2/}
}


\references{
  Hamill TM, Bates GT, Whitaker JS, Murray DR, Fiorino M, Galarneau Jr. TJ, Zhu Y, Lapenta W (2013).
  NOAA's Second-Generation Global Medium-Range Ensemble Reforecast Data Set.
  \emph{Bulletin of the American Meteorological Society}, \bold{94}(10), 1553--1565.
  \doi{10.1175/BAMS-D-12-00014.1}
}


\examples{
## load data and omit missing values
data("IbkTemperature", package = "lmSubsets")
IbkTemperature <- na.omit(IbkTemperature)

## fit a simple climatological model for the temperature
## with a linear trend and annual/bi-annual harmonic seasonal pattern
CLIM <- lm(temp ~ time + sin + cos + sin2 + cos2,
  data = IbkTemperature)

## fit a simple MOS with 2-meter temperature forecast in addition
## to the climatological model
MOS0 <- lm(temp ~ t2m + time + sin + cos + sin2 + cos2,
  data = IbkTemperature)

## graphical comparison and MOS summary
plot(temp ~ time, data = IbkTemperature, type = "l", col = "darkgray")
lines(fitted(MOS0) ~ time, data = IbkTemperature, col = "darkred")
lines(fitted(CLIM) ~ time, data = IbkTemperature, lwd = 2)
MOS0

## best subset selection of remaining variables for the MOS
## (i.e., forcing the regressors of m1 into the model)
MOS1_all <- lmSubsets(temp ~ ., data = IbkTemperature,
  include = c("t2m", "time", "sin", "cos", "sin2", "cos2"))
plot(MOS1_all)
image(MOS1_all, size = 8:20)
## -> Note that soil temperature and maximum temperature are selected
## in addition to the 2-meter temperature

## best subset selection of all variables
MOS2_all <- lmSubsets(temp ~ ., data = IbkTemperature)
plot(MOS2_all)
image(MOS2_all, size = 2:20)
## -> Note that 2-meter temperature is not selected into the best
## BIC model but soil-temperature (and maximum temperature) are used instead

## refit the best BIC subset selections
MOS1 <- refit(lmSelect(MOS1_all))
MOS2 <- refit(lmSelect(MOS2_all))

## compare BIC
BIC(CLIM, MOS0, MOS1, MOS2)

## compare RMSE
sqrt(sapply(list(CLIM, MOS0, MOS1, MOS2), deviance)/
  nrow(IbkTemperature))

## compare coefficients
cf0 <- coef(CLIM)
cf1 <- coef(MOS0)
cf2 <- coef(MOS1)
cf3 <- coef(MOS2)
names(cf2) <- gsub("^x", "", names(coef(MOS1)))
names(cf3) <- gsub("^x", "", names(coef(MOS2)))
nam <- unique(c(names(cf0), names(cf1), names(cf2), names(cf3)))
cf <- matrix(NA, nrow = length(nam), ncol = 4,
  dimnames = list(nam, c("CLIM", "MOS0", "MOS1", "MOS2")))
cf[names(cf0), 1] <- cf0
cf[names(cf1), 2] <- cf1
cf[names(cf2), 3] <- cf2
cf[names(cf3), 4] <- cf3
print(round(cf, digits = 3), na.print = "")
}


\keyword{datasets}
