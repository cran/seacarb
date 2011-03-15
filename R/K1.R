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
"K1" <-
function(S=35,T=25,P=0,k1k2='x',pHscale="T")
{

nK <- max(length(S), length(T), length(P), length(k1k2), length(pHscale))

##-------- Creation de vecteur pour toutes les entrees (si vectorielles)

if(length(S)!=nK){S <- rep(S[1], nK)}
if(length(T)!=nK){T <- rep(T[1], nK)}
if(length(P)!=nK){P <- rep(P[1], nK)}
if(length(k1k2)!=nK){k1k2 <- rep(k1k2[1], nK)}
if(length(pHscale)!=nK){pHscale <- rep(pHscale[1], nK)}

##---------- pHsc : this vector note the actual pHscale because it can change during processing
pHsc <- rep(NA,nK)
pHlabel <- rep(NA,nK)
K1 <- rep(NA, nK)

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

# --------------------- K1 ---------------------------------------
#   first acidity constant:
#   [H^+] [HCO_3^-] / [CO2] = K_1
#
#     Mehrbach et al (1973) refit by Lueker et al. (2000).
#
#(Lueker  et al., 2000 in Guide to the Best Practices for Ocean CO2 Measurements
#   Dickson, Sabine and Christian , 2007, Chapter 5, p. 13)
#
#   pH-scale: 'total'. mol/kg-soln
if(k1k2[i] == "l"){
logK1lue <- (-3633.86)/TK[i] + 61.2172 - 9.67770*log(TK[i]) + 0.011555*S[i] - 0.0001152*S[i]*S[i]
K1[i]<- 10^logK1lue
pHsc[i] <- "T"
}

# --------------------- K1 ---------------------------------------
#   first acidity constant:
#   [H^+] [HCO_3^-] / [CO2] = K_1
#
#   (Roy et al., 1993 in Dickson and Goyet, 1994, Chapter 5, p. 14)
#   pH-scale: 'total'. mol/kg-soln
if(k1k2[i] == "r"){
tmp1  <- 2.83655 - 2307.1266 / TK[i] - 1.5529413 * log(TK[i])
tmp2  <- - (0.20760841 + 4.0484 / TK[i]) * sqrt(S[i])
tmp3 <- 0.08468345 * S[i] - 0.00654208 * S[i] * sqrt(S[i])   
tmp4 <- log(1 - 0.001005 * S[i])
lnK1roy <- tmp1 + tmp2 + tmp3 + tmp4
K1[i] <- exp(lnK1roy)
pHsc[i] <- "T"
}
        
# --------------------- K1 ---------------------------------------
#   first acidity constant:
#   [H^+] [HCO_3^-] / [CO2] = K_1
#
#   Millero et al. 2006 Marine Chemistry
#   pH-scale: 'SWS scale'. mol/kg-soln
if(k1k2[i] == "m06"){
pK1o <- 6320.813/TK[i] + 19.568224*log(TK[i]) -126.34048
A1 <- 13.4191*S[i]^(0.5) + 0.0331*S[i] - (5.33e-5)*S[i]^2
B1 <- -530.123*S[i]^(0.5) - 6.103*S[i]
C1 <- -2.06950*S[i]^(0.5)
pK1 <- pK1o + A1 + B1/TK[i] + C1*log(TK[i])
K1[i] <- 10^(-pK1)
pHsc[i] <- "SWS"
}

# --------------------- K1 ---------------------------------------
#   first acidity constant:
#   [H^+] [HCO_3^-] / [CO2] = K_1
#
#   Millero 2010 Marine and Fresh water research

if(k1k2[i] == "m10"){
pK1o <- 6320.813/TK[i] + 19.568224*log(TK[i]) -126.34048

#   pH-scale: 'SWS scale'. mol/kg-soln
if(pHscale[i]=="SWS" | P[i]>0){
A1 <- 13.4038*S[i]^(0.5) + 0.03206*S[i] - (5.242e-5)*S[i]^2
B1 <- -530.659*S[i]^(0.5) - 5.8210*S[i]
C1 <- -2.0664*S[i]^(0.5)
pK1 <- pK1o + A1 + B1/TK[i] + C1*log(TK[i])
K1[i] <- 10^(-pK1)   # K1 according to Millero et al. 2010 at Seawater scale
pHsc[i] <- "SWS" 
}

#   pH-scale: 'Total scale'. mol/kg-soln
if(pHscale[i]=="T" & P[i]==0){
A1 <- 13.4051*S[i]^(0.5) + 0.03185*S[i] - (5.218e-5)*S[i]^2
B1 <- -531.095*S[i]^(0.5) - 5.7789*S[i]
C1 <- -2.0663*S[i]^(0.5)
pK1 <- pK1o + A1 + B1/TK[i] + C1*log(TK[i])
K1[i] <- 10^(-pK1)   # K1 according to Millero et al. 2010 at Total scale 
pHsc[i] <- "T"
}

#   pH-scale: 'Free scale'. mol/kg-soln
if(pHscale[i]=="F" & P[i]==0){
A1 <- 5.09247*S[i]^(0.5) + 0.05574*S[i] - (9.279e-5)*S[i]^2
B1 <- -189.879*S[i]^(0.5) - 11.3108*S[i]
C1 <- -0.8080*S[i]^(0.5)
pK1 <- pK1o + A1 + B1/TK[i] + C1*log(TK[i])
K1[i] <- 10^(-pK1)  # K1 according to Millero et al. 2010 at Total scale
pHsc[i] <- "F"
}
}
}



# ------------------- Pression effect --------------------------------
for(i in (1:nK)){
if (P[i] > 0.0){    ## pressure correction
if (pHsc[i] == "T"){    ## conversion of K1 from Total scale to SWS scale in the case of Luecker or Roy calculation
factor <- kconv(S=S[i], T=T[i], P=0)$ktotal2SWS
K1[i] <- K1[i] * factor
pHsc[i] <- "SWS"
}
}
}

K1 <- Pcorrect(Kvalue=K1, Ktype="K1", T=T, S=S, P=P, pHscale=pHsc)

## --------------- Last conversion in the require pHscale ----------------

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
K1[i] <- K1[i]*factor
}


##---------------Attributes
method <- c()
for(i in 1:nK){
m <- "Luecker et al. (2000)"
if(k1k2[i]=="m06"){m <- "Millero et al. (2006)"}
if(k1k2[i]=="m10"){m <- "Millero (2010)"}
if(k1k2[i]=="r"){m <- "Roy et al. (1993)"}
method <- c(method, m)
}

##------------Warnings

for(i in 1:nK){
if((k1k2[i]=='l')&((T[i]>35)|(T[i]<2)|(S[i]<19)|(S[i]>43))){warning("S and/or T is outside the range of validity of the formulation chosen for K1.")}
if((k1k2[i]=='r')&((T[i]>45)|(S[i]>45))){warning("S and/or T is outside the range of validity of the formulation chosen for K1.")}
if((T[i]>50)|(S[i]>50)){warning("S and/or T is outside the range of validity of the formulations available for K1 in seacarb.")}
}

attr(K1,"unit") <- "mol/kg-soln"
attr(K1,"pH scale") <- pHlabel
attr(K1, "method") <- method
return(K1)
}
