\encoding{latin1}
\name{pHinsi}
\alias{pHinsi}
\title{pH at in situ temperature}
\description{pH at in situ temperature}
\usage{
pHinsi(PH=8.2, ALK=2.4e-3, Tinsi=20, Tlab=25, S=35, Pt=0, Sit=0, k1k2 = "x", 
	kf = "x", pHscale = "T")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PH}{pH measured in the laboratory}
  \item{ALK}{ALK, total alkalinity (mol/kg)}
  \item{Tinsi}{In situ temperature in degrees Celsius}
  \item{Tlab}{Measurement temperature in degrees Celsius}
  \item{S}{Salinity}
  \item{Pt}{value of the concentration of total phosphate in mol/kg}
  \item{Sit}{the value of the total silicate in mol/kg}
  \item{k1k2}{"l" for using K1 and K2 from Lueker et al. (2000), "m" from Millero et al. (2006) and "r" from Roy et al. (1993). "x" is the default flag; the default value is then "l", except if T is outside the range 2 to 35oC and/or S is outside the range 19 to 43. In these cases, the default value is "m". }
	\item{kf}{"pf" for using Kf from Perez and Fraga (1987) and "dg" for using Kf from Dickson and Riley (1979). "x" is the default flag; the default value is then "pf", except if T is outside the range 9 to 33oC and/or S is outside the range 10 to 40. In these cases, the default is "dg".}
	\item{pHscale}{choice of pH scale: "T" for the total scale, "F" for the free scale and "SWS" for using the seawater scale, default is "T" (total scale)}}

\value{
  \item{pH}{pH at in situ temperature}
}

\details{The Lueker et al. (2000) constant is recommended by Guide to Best Practices for Ocean CO2 Measurements (2007). It is critical to consider that each formulation is valid in specific ranges of temperature and salinity:

\emph{For K1 and K2:}
\itemize{
\item Roy et al. (1993): S ranging between 0 and 45 and T ranging between 0 and 45oC.

\item Lueker et al. (2000): S ranging between 19 and 43 and T ranging between 2 and 35oC.

\item Millero et al. (2006): S ranging between 0 and 50 and T ranging between 0 and 50oC.
}

\emph{For Kh:}
\itemize{
\item Perez and Fraga (1987): S ranging between 10 and 40 and T ranging between 9 and 33oC.

\item Dickson and Goyet (1979): S ranging between 0 and 45 and T ranging between 0 and 45oC.
}
The arguments can be given as a unique number or as vectors. If the lengths of the vectors are different, the longer vector is retained and only the first value of the other vectors is used. It is recommended to use either vectors with the same dimension or one vector for one argument and numbers for the other arguments.
}

\references{
Hunter K. A., 1998. The temperature dependence of pH in surface seawater. \emph{Deep-Sea Research (Part I, Oceanographic Research Papers)} \bold{45}(11):1919-1930.

Millero, F. J. Graham, T. B. Huang, F. Bustos-Serrano, H. and Pierrot, D. 2006. Dissociation constants of carbonic acid in seawater as a function of salinity and temperature.  \emph{Marine Chemistry} \bold{100}, 80-84.
}
             
\author{
Jean-Pierre Gattuso, \email{gattuso@obs-vlfr.fr}
}

%\seealso{	\code{\link{carb}}}

\examples{
  pHinsi(8.2,2.4e-3,25,25,35,0,0)
}

\keyword{utilities}
