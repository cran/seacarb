# Copyright (C) 2010 Jean-Pierre Gattuso and Heloise Lavigne
#
# This file is part of seacarb.
#
# Seacarb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.
#
# Seacarb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with seacarb; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#

"Pcorrect" <-
function(Kvalue, Ktype, T=25, S=35, P=0, pHscale="T"){

nK <- max(length(Kvalue), length(Ktype), length(P), length(T), length(pHscale), length(S))

##-------- Creation de vecteur pour toutes les entrees (si vectorielles)

if(length(Kvalue)!=nK){Kvalue <- rep(Kvalue[1], nK)}
if(length(Ktype)!=nK){Ktype <- rep(Ktype[1], nK)}
if(length(T)!=nK){T <- rep(T[1], nK)}
if(length(S)!=nK){S <- rep(S[1], nK)}
if(length(P)!=nK){P <- rep(P[1], nK)}
if(length(pHscale)!=nK){pHscale <- rep(pHscale[1], nK)}

# Constants 
TK <- T + 273.15   # Temperature in Kelvin
R = 83.14472;             # mol bar deg-1
	
## loading of table with coefficients
if(!exists("Pcoeffs", where = .GlobalEnv)) data(Pcoeffs)

# ------------------- Pression effect --------------------------------
for(i in (1:nK)){
if (P[i] > 0.0) {

## 
##  here, pHscale given in argument is the pH scale of Kvalue also given in argument
##  if pHscale is different from seawater scale, it is necessary to convert K in seawater scale
##  before applying the pressure correction, except with Kf, Kspc or Kspa.
##  whith Kf, Kspc and Kspa pressure correction applied on free scale

# In case of pressure correction on the seawater scale
if (Ktype[i] %in% c("K1", "K2", "K1p", "K2p", "K3p", "Kb", "Khs", "Kn", "Ksi", "Kw")){
if (pHscale[i] == "T") { Kvalue[i] <- Kvalue[i]*kconv(S=S[i], T=T[i], P=0)$ktotal2SWS }
if (pHscale[i] == "F") { Kvalue[i] <- Kvalue[i]*kconv(S=S[i], T=T[i], P=0)$kfree2SWS }

    l <- which(Pcoeffs$K == Ktype[i])
    deltav  <-  Pcoeffs$a0[l] + Pcoeffs$a1[l] *T[i] + Pcoeffs$a2[l] *T[i]*T[i]
	  deltak  <-  Pcoeffs$b0[l]  + Pcoeffs$b1[l] *T[i] + Pcoeffs$b2[l] *T[i]*T[i]
	  lnkpok0 <-  -(deltav /(R*TK[i]))*P[i] + (0.5*deltak /(R*TK[i]))*P[i]*P[i];

	Kvalue[i] = Kvalue[i]*exp(lnkpok0);
	
## Then Kvalue is converted back in the pHscale ginven in argument
if (pHscale[i] == "T") { Kvalue[i] <- Kvalue[i]*kconv(S=S[i], T=T[i], P=300)$kSWS2total }
if (pHscale[i] == "F") { Kvalue[i] <- Kvalue[i]*kconv(S=S[i], T=T[i], P=300)$kSWS2free }

}

# In case of there is no pH scale choice 
# Ks -> "free scale"
# Kspa and Kspc do not need pH scale
if (Ktype[i] %in% c("Ks", "Kspa", "Kspc", "Kf")){

if(Ktype[i]=="Kf"){ #conversion on free scale
ST <- 0.14/96.062/1.80655*S[i]
FT <- 7e-5*(S[i]/35)
if(pHscale[i]=="T"){Kvalue[i] <- Kvalue[i] /(1+ST/Ks(S=S[i], T=T[i], P=0))}
if(pHscale[i]=="SWS"){Kvalue[i] <- Kvalue[i] / (1 + ST/Ks(S=S[i], T=T[i], P=0)+ FT/Kf(T=T[i], P=0, S=S[i], pHscale="F"))}
}
    l <- which(Pcoeffs$K == Ktype[i])
    deltav  <-  Pcoeffs$a0[l] + Pcoeffs$a1[l] *T[i] + Pcoeffs$a2[l] *T[i]*T[i]
    deltak  <-  Pcoeffs$b0[l]  + Pcoeffs$b1[l] *T[i] + Pcoeffs$b2[l] *T[i]*T[i]
    lnkpok0 <-  -(deltav /(R*TK[i]))*P[i] + (0.5*deltak /(R*TK[i]))*P[i]*P[i];
    
    Kvalue[i] = Kvalue[i]*exp(lnkpok0);
    if(Ktype[i] %in% c("Ks", "Kf")){pHscale[i] <- "F" }
}


}
}

phs <- rep(NA, nK)
for(i in 1:nK){
if(pHscale[i]=="SWS"){phs[i] <- "Seawater scale"}
if(pHscale[i]=="T"){phs[i] <- "Total scale"}
if(pHscale[i]=="F"){phs[i] <- "Free scale"}
if (Ktype[i] %in% c("Kspa", "Kspc")){phs[i] <- NA}
}



attr(Kvalue, "pH scale") <- phs
attr(Kvalue, "unit") <- "mol/kg-soln"
return(Kvalue)
}
