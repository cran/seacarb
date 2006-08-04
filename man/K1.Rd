\name{K1}
\alias{K1}
\title{First dissociation constant of carbonic acid (mol/kg)}
\description{First dissociation constant of carbonic acid (mol/kg)}
\usage{
K1(S = 35, T = 25, P = 0, k1k2 = "r", phflag = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
  \item{k1k2}{"r" for using K1 and K2 from Roy et al. and "m" for using K1 and K2 from Mehrbach, default is "r"}
  \item{phflag}{The pH scale is either the total scale (0) or the free scale (1), default is 0}
}

\value{
  \item{K1}{First dissociation constant of carbonic acid (mol/kg)}
}

\references{The Roy et al. (1993) constant is recommended by DOE (1994). The Merbach et al. (1973) constants used are on the total scale from Lueker et al. (2000).

DOE 1994 \emph{Handbook of methods for the analysis of the various parameters of the carbon dioxide system in sea water}. ORNL/CDIAC-74. Oak Ridge,Tenn.: Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory.

Lueker, T. J. Dickson, A. G. and Keeling, C. D. 2000 Ocean pCO2 calculated from dissolved inorganic carbon, alkalinity, and equations for K1 and K2: validation based on laboratory measurements of CO2 in gas and seawater at equilibrium. \emph{Marine Chemistry} \bold{70} 105-119.

Mehrbach, C. Culberson, C. H. Hawley, J. E. and Pytkowicz, R. M. 1973 Measurement of the apparent dissociation constants of carbonic acid in seawater at atmospheric pressure. \emph{Limnology and Oceanography} \bold{18}, 897-907.

Roy, R. N. Roy, L. N. Vogel, K. M. Porter-Moore, C. Pearson, T. Good C. E., Millero F. J. and Campbell D. M., 1993. The dissociation constants of carbonic acid in seawater at salinities 5 to 45 and temperatures 0 to 45¡C. \emph{Marine Chemistry} \bold{44}, 249-267.
}
             
\author{Aurelien Proye and Jean-Pierre Gattuso}

\seealso{K2.rd}

\examples{
  K1(35,25,0)
}

\keyword{utilities}
