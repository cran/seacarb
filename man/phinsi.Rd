\encoding{latin1}
\name{phinsi}
\alias{phinsi}
\title{pH at in situ temperature}
\description{pH at in situ temperature}
\usage{
phinsi(PH=8.2, ALK=2.4e-3, Tinsi=20, Tlab=25, S=35, Pt=0, Sit=0)
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
}

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
  phinsi(8.2,2.4e-3,25,25,35,0,0)
}

\keyword{utilities}
