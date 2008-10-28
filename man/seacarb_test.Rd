\name{seacarb_test}
\docType{data}
\alias{seacarb_test}
\title{Test data file to test the use of the carb function}
\description{
The variables are:

\item Flag indicating which couple of variables is used
\item Value of the first variable in mol/kg, except for pH and for pCO2 in \eqn{\mu}atm
\item Value of the second variable in mol/kg, except for pH
\item Salinity
\item Temperature in degrees Celsius
\item Hydrostatic pressure in bar (surface = 0)
\item Value of the concentration of total phosphate in mol/kg
\item Value of the total silicate in mol/kg
}
\usage{seacarb_test}
\format{A data frame with 20 rows and 8 variables}
\source{None, these data were invented for this purpose. The input variables were chosen in order to check that the carbonate chemistry is identical for all flags.}
\keyword{datasets}
