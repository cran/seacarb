# Copyright (C) 2008 Jean-Pierre Gattuso and Heloise Lavigne and Aurelien Proye
#
# This file is part of seacarb.
#
# Seacarb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.
#
# Seacarb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with seacarb; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#
"K2" <-
function(S=35,T=25,P=0,k1k2='x',pHscale="T"){

nK <- max(length(S), length(T), length(P), length(k1k2), length(pHscale))

##-------- Creation de vecteur pour toutes les entrees (si vectorielles)

if(length(S)!=nK){S <- rep(S[1], nK)}
if(length(T)!=nK){T <- rep(T[1], nK)}
if(length(P)!=nK){P <- rep(P[1], nK)}
if(length(k1k2)!=nK){k1k2 <- rep(k1k2[1], nK)}
if(length(pHscale)!=nK){pHscale <- rep(pHscale[1], nK)}

##----------  Initialisation of vectors
pHsc <- rep(NA,nK) # pHsc : this vector note the actual pHscale because it can change during processing
pHlabel <- rep(NA,nK)   
K2 <- rep(NA, nK)

##----------Check the validity of the method regarding the T/S range

for(i in 1:nK){
if(k1k2[i]=='x'){
k1k2[i] <- 'l'  ## luecker by default
if((T[i]>35)|(T[i]<2)|(S[i]<19)|(S[i]>43)){k1k2[i] <- 'm10' }
}
}

#-------Constantes----------------

#---- issues de equic----
tk = 273.15;           # [K] (for conversion [deg C] <-> [K])
TK = T + tk;           # TC [C]; T[K]

for(i in 1:nK){

# --------------------- K2 lueker et al. 2000 ------------------------------
#
#   second acidity constant:
#   [H^+] [CO_3^--] / [HCO_3^-] = K_2
#
#   Mehrbach et al. (1973) refit by Lueker et al. (2000).
#
#   (Lueker  et al., 2000 in Guide to the Best Practices for Ocean CO2 Measurements
#   Dickson, Sabin and Christian , 2007, Chapter 5, p. 14)
#
#   pH-scale: 'total'. mol/kg-soln

if(k1k2[i] == "l"){
logK2lue <- -471.78/TK[i] - 25.9290 + 3.16967*log(TK[i]) + 0.01781*S[i] - 0.0001122*S[i]*S[i]
K2[i] <- 10^(logK2lue)
pHsc[i] <- "T"
}

# --------------------- K2 Roy et al. 1993----------------------------------------
#
#   second acidity constant:
#   [H^+] [CO_3^--] / [HCO_3^-] = K_2
#
#   (Roy et al., 1993 in Dickson and Goyet, 1994, Chapter 5, p. 15)
#   pH-scale: 'total'. mol/kg-soln

if(k1k2[i] == "r"){	
tmp1 = -9.226508 - 3351.6106 / TK[i] - 0.2005743 * log(TK[i]);
tmp2 = (-0.106901773 - 23.9722 / TK[i]) * sqrt(S[i]);
tmp3 = 0.1130822 * S[i] - 0.00846934 * S[i]^1.5 + log(1 - 0.001005 * S[i]);
lnK2roy = tmp1 + tmp2 + tmp3;
K2[i] <- exp(lnK2roy)
pHsc[i] <- "T"
}

# --------------------- K2 ---------------------------------------
#   first acidity constant:
#   [H^+] [CO_3^--] / [HCO_3^-] = K_2
#
#   Millero et al. 2006 Marine Chemistry
#   pH-scale: 'SWS scale'. mol/kg-soln

if(k1k2[i] == "m06"){
pK2o <- -90.18333 + 5143.692/TK[i] + 14.613358*log(TK[i]) 
A2 <- 21.0894*S[i]^(0.5) + 0.1248*S[i] - (3.687e-4)*S[i]^2
B2 <- -772.483*S[i]^(0.5) - 20.051*S[i]
C2 <- -3.3336*S[i]^(0.5)
pK2 <- pK2o + A2 + B2/TK[i] + C2*log(TK[i])
K2[i] <- 10^(-pK2)
pHsc[i] <- "SWS"
}

# --------------------- K2 ---------------------------------------
#   first acidity constant:
#   [H^+] [CO_3^--] / [HCO_3^-] = K_2
#
#   Millero 2010 Marine and Fresh water research

if(k1k2[i]=="m10"){
pK2o <- 5143.692/TK[i] + 14.613358*log(TK[i]) -90.18333

#   pH-scale: 'SWS scale'. mol/kg-soln
if(pHscale[i]=="SWS" | P[i]>0){
A2 <- 21.3728*S[i]^(0.5) + 0.1218*S[i] - (3.688e-4)*S[i]^2
B2 <- -788.289*S[i]^(0.5) - 19.189*S[i]
C2 <- -3.374*S[i]^(0.5)
pK2 <- pK2o + A2 + B2/TK[i] + C2*log(TK[i])
K2[i] <- 10^(-pK2)
pHsc[i] <- "SWS"
}

#   pH-scale: 'Total scale'. mol/kg-soln
if(pHscale[i]=="T" & P[i]==0){
A2 <- 21.5724*S[i]^(0.5) + 0.1212*S[i] - (3.714e-4)*S[i]^2
B2 <- -798.292*S[i]^(0.5) - 18.951*S[i]
C2 <- -3.403*S[i]^(0.5)
pK2 <- pK2o + A2 + B2/TK[i] + C2*log(TK[i])
K2[i] <- 10^(-pK2)
pHsc[i] <- "T"
}

#   pH-scale: 'Free scale'. mol/kg-soln
if(pHscale[i]=="F" & P[i]==0){
A2 <- 11.0637*S[i]^(0.5) + 0.1379*S[i] - (3.688e-4)*S[i]^2
B2 <- -366.178*S[i]^(0.5) - 23.288*S[i]
C2 <- -1.810*S[i]^(0.5)
pK2 <- pK2o + A2 + B2/TK[i] + C2*log(TK[i])
K2[i] <- 10^(-pK2)
pHsc[i] <- "F"
}
}
}

# ------------------- Pression effect --------------------------------
for(i in (1:nK)){
if (P[i] > 0.0){    ## pressure correction
if (pHsc[i] == "T"){    ## conversion of K1 from Total scale to SWS scale in the case of Luecker or Roy calculation
factor <- kconv(S=S[i], T=T[i], P=0)$ktotal2SWS
K2[i] <- K2[i] * factor
pHsc[i] <- "SWS"
}
}
}

K2 <- Pcorrect(Kvalue=K2, Ktype="K2", S=S, T=T, P=P, pHscale=pHsc)

## ------------------- Last conversion in the require pHscale ------------------

for(i in (1:nK)){
## In case of pressure correction (K2 is in the SWS scale)
if (P[i] > 0.0){     
 if(pHscale[i]=="T"){factor <- kconv(S=S[i], T=T[i], P=P[i])$kSWS2total; pHlabel[i] <- "total scale"}
 if(pHscale[i]=="F"){factor <- kconv(S=S[i], T=T[i], P=P[i])$kSWS2free ; pHlabel[i] <- "free scale"}
 if(pHscale[i]=="SWS"){factor <- 1 ; pHlabel[i] <- "seawater scale"}
 }
## In case of no pressure correction (the pH scale of K2 depends on the calculation method)
 if(P[i]==0){ 
    if(pHsc[i] == "T"){  # with the Luecker, Roy and Millero 2010 formulations
      if(pHscale[i]=="T"){factor <- 1; pHlabel[i] <- "total scale"}
      if(pHscale[i]=="F"){factor <- kconv(S=S[i], T=T[i], P=P[i])$ktotal2free ; pHlabel[i] <- "free scale"}
      if(pHscale[i]=="SWS"){factor <- kconv(S=S[i], T=T[i], P=P[i])$ktotal2SWS  ; pHlabel[i] <- "seawater scale"}
    }
    if(pHsc[i] == "SWS"){ # whith the Millero 2006 and Millero 2010 formulations 
      if(pHscale[i]=="T"){factor <- kconv(S=S[i], T=T[i], P=P[i])$kSWS2total; pHlabel[i] <- "total scale"}
      if(pHscale[i]=="F"){factor <- kconv(S=S[i], T=T[i], P=P[i])$kSWS2free ; pHlabel[i] <- "free scale"}
      if(pHscale[i]=="SWS"){factor <- 1 ; pHlabel[i] <- "seawater scale"}
    }
    if(pHsc[i]=="F"){ # whith the Millero 2010 formulation
    factor <- 1
    pHlabel[i] <- "free scale"
    }
    }
K2[i] <- K2[i]*factor
}

##------------Warnings

for(i in 1:nK){
if((k1k2[i]=='l')&((T[i]>35)|(T[i]<2)|(S[i]<19)|(S[i]>43))){warning("S and/or T is outside the range of validity of the formulation used for K2.")}
if((k1k2[i]=='r')&((T[i]>45)|(S[i]>45))){warning("S and/or T is outside the range of validity of the formulation used for K2.")}
if((T[i]>50)|(S[i]>50)){warning("S and/or T is outside the range of validity of the formulations available for K2 in seacarb.")}
}

##----------Attributes

method <- c()
for(i in 1:nK){
m <- "Luecker et al. (2000)"
if(k1k2[i]=="m06"){m <- "Millero et al. (2006)"}
if(k1k2[i]=="m10"){m <- "Millero (2010)"}
if(k1k2[i]=="r"){m <- "Roy et al. (1993)"}
method <- c(method, m)
}

attr(K2,"unit")     = "mol/kg-soln"
attr(K2,"pH scale") = pHlabel
attr(K2,"method") = method
return(K2)
}



	
