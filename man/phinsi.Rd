\encoding{latin1}
\name{pHinsi}
\alias{pHinsi}
\title{pH at in situ temperature}
\description{pH at in situ temperature}
\usage{
pHinsi(PH=8.2, ALK=2.4e-3, Tinsi=20, Tlab=25, S=35, Pt=0, Sit=0, k1k2 = "l", kf = "pf", pHscale = "T")
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
	\item{k1k2}{"l" for using K1 and K2 from Lueker et al. and "r" for using K1 and K2 from Roy et al. , default is "l"}
	\item{kf}{"pf" for using Kf from Perez and Fraga (1987) and "dg" for using Kf from Dickson and Goyet (1979), default is "pf"}
	\item{pHscale}{choice of pH scale: "T" for the total scale, "F" for the free scale and "SWS" for using the seawater scale, default is "T" (total scale)}}

\value{
  \item{pH}{pH at in situ temperature}
}

\references{
Hunter K. A., 1998. The temperature dependence of pH in surface seawater. \emph{Deep-Sea Research (Part I, Oceanographic Research Papers)} \bold{45}(11):1919-1930.
}
             
\author{
Jean-Pierre Gattuso, \email{gattuso@obs-vlfr.fr}
}

%\seealso{	\code{\link{carb}}}

\examples{
  pHinsi(8.2,2.4e-3,25,25,35,0,0)
}

\keyword{utilities}
