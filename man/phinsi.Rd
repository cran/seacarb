\name{phinsi}
\alias{phinsi}
\title{pH at in situ temperature}
\description{pH at in situ temperature}
\usage{
phinsi(PH=8.2, ALK=2.4e-3, Tinsi=20, Tlab=25, S=35, k1k2='r', phflag=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PH}{pH measured in the laboratory}
  \item{ALK}{ALK, total alkalinity (mol/kg)}
  \item{Tinsi}{In situ temperature in degrees Celsius}
  \item{Tlab}{Measurement temperature in degrees Celsius}
  \item{S}{Salinity}
  \item{k1k2}{"r" for using K1 and K2 from Roy et al. and "m" for using K1 and K2 from Mehrbach,}
  \item{phflag}{The pH scale is either the total scale (0) or the free scale (1), default is 0}
}

\value{
  \item{pH}{pH at in situ temperature}
}

\references{
Hunter K. A., 1998. The temperature dependence of pH in surface seawater. \emph{Deep-Sea Research (Part I, Oceanographic Research Papers)} \bold{45}(11):1919-1930.
}
             
\author{Jean-Pierre Gattuso}

%\seealso{	\code{\link{carb}}}

\examples{
  phinsi(8.2,2.4e-3,25,25,35,"r",0)
}

\keyword{utilities}
