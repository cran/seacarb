\encoding{latin1}
\name{Kf}
\alias{Kf}

\title{Stability constant of hydrogen fluoride (mol/kg)}
\description{Stability constant of hydrogen fluoride (mol/kg)}
\usage{
Kf(S = 35, T = 25, P = 0, kf = "pf", pHscale="T")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
  \item{kf}{"pf" for using Kf from Perez and Fraga (1987) and "dg" for using Kf from Dickson and Goyet (1979), default is "pf"}
  \item{pHscale}{choice of pH scale: "T" for the total scale, "F" for the free scale and "SWS" for using the seawater scale, default is "T" (total scale)}
}

\details{Note that the arguments can be given as a unique number or as vectors. If the lengths of the vectors are different, the longer vector is retained and only the first value of the other vectors is used. It can therefore be critical to use vectors of the same length.} 

\value{
  \item{Kf}{Stability constant of hydrogen fluoride (mol/kg)}

}
\references{The Perez and Fraga (1987) constant is recommended by Guide to Best Practices for Ocean CO2 Measurements (2007). The Dickson and Goyet (1979) constants is recommended by DOE (1994).

DOE 1994 \emph{Handbook of methods for the analysis of the various parameters of the carbon dioxide system in sea water}. ORNL/CDIAC-74. Oak Ridge,Tenn.: Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory.

Dickson, A. G. and Riley, J. P. 1979 The estimation of acid dissociation constants in seawater media from potentiometric titrations with strong base. I. The ionic product of water. \emph{Marine Chemistry} \bold{7}, 89-99.

Dickson A. G., Sabine C. L. and Christian J. R., 2007 Guide to best practices for ocean CO2 measurements. \emph{PICES Special Publication} \bold{3}, 1-191.

Perez, F. F. and Fraga, F. 1987 Association constant of fluoride and hydrogen ions in seawater. \emph{Marine Chemistry} \bold{21}, 161-168.

}

\author{
Héloïse Lavigne, Aurélien Proye and Jean-Pierre Gattuso \email{gattuso@obs-vlfr.fr}
}


\examples{
  Kf(S=35,T=25,P=0,kf="pf",pHscale="T")
}
\keyword{utilities}
