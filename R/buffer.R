# Copyright (C) 2008 Jean-Pierre Gattuso and Heloise Lavigne and Aurelien Proye
# with a most valuable contribution of Bernard Gentili <gentili@obs-vlfr.fr>
# and valuable suggestions from Jean-Marie Epitalon <epitalon@lsce.saclay.cea.fr>
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

buffer <- 
function(flag, var1, var2, S=35, T=25, P=0, Pt=0, Sit=0, k1k2='l', kf='pf', pHscale="T"){
Carb <- carb(flag=flag, var1=var1, var2=var2, S=S, T=T, P=P, Pt=Pt, Sit=Sit, k1k2=k1k2, kf=kf, pHscale=pHscale)
RES <- data.frame()

n <- nrow(Carb)

k1k2.V <- k1k2
if(length(k1k2.V)!=n){k1k2.V <- rep(k1k2.V[1], n)}
kf.V <- kf
if(length(kf.V)!=n){kf.V <- rep(kf.V[1], n)}
pHscale.V <- pHscale
if(length(pHscale.V)!=n){pHscale.V <- rep(pHscale.V[1], n)}

	for (i in (1:nrow(Carb))){

	S <- Carb[i,2]
	T <- Carb[i,3]
	P <- Carb[i,4]
	PH <- Carb[i,5]
	h <- 10^(-PH)
	CO2 <- Carb[i,6]
	pCO2 <- Carb[i,7]
	fCO2 <- Carb[i,8]
	HCO3 <- Carb[i,9]
	CO3 <- Carb[i,10]
	DIC <- Carb[i,11]
	ALK <- Carb[i,12]
	Oa <- Carb[i,13]
	Oc <- Carb[i,14]

	k1k2 <- k1k2.V[i]
	kf <- kf.V[i]
	pHscale <- pHscale.V[i]

#-------Constantes----------------

tk = 273.15;           # [K] (for conversion [deg C] <-> [K])

# JME: moved following code block here, after reading imput file

TK = T + tk;           # TK [K]; T[C]

#---- issues de equic----
Cl = S / 1.80655;      # Cl = chlorinity; S = salinity (per mille)
cl3 = Cl^(1/3);
ION = 0.00147 + 0.03592 * Cl + 0.000068 * Cl * Cl;   # ionic strength
iom0 = 19.924*S/(1000-1.005*S);
ST = 0.14/96.062/1.80655*S;   # (mol/kg soln) total sulfate
bor = (416.*(S/35.))* 1e-6;   # (mol/kg), DOE94 boron total
fluo = (7*(S/35))*1e-5        # (mol/kg), DOE94 fluoride total

#---------------------------------------------------------------------
#--------------------- calcul des K ----------------------------------
#---------------------------------------------------------------------

K1 <- K1(S=S, T=T, P=P, pHscale=pHscale, k1k2=k1k2)   
K2 <- K2(S=S, T=T, P=P, pHscale=pHscale, k1k2=k1k2)
Ks <- Ks(S=S, T=T, P=P)
Kf <- Kf(S=S, T=T, P=P, pHscale=pHscale, kf=kf)
Kw <- Kw(S=S, T=T, P=P, pHscale=pHscale)
Kh <- Kh(S=S, T=T, P=P)
Kb <- Kb(S=S, T=T, P=P, pHscale=pHscale)
K1p <- K1p(S=S, T=T, P=P, pHscale=pHscale)
K2p <- K2p(S=S, T=T, P=P, pHscale=pHscale)
K3p <- K3p(S=S, T=T, P=P, pHscale=pHscale)
Ksi <- Ksi(S=S, T=T, P=P, pHscale=pHscale)
Kspa <- Kspa(S=S, T=T, P=P)
Kspc <- Kspc(S=S, T=T, P=P)
	
rho <- rho(S=S,T=T,P=P)

#---------------------------------------------------------------------
#--------------------    buffer effects    ---------------------------
#---------------------------------------------------------------------

	DD=-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DD)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DD)/((h+2*K2)*(h+2*K2));
	PhiD=-1/(h*log(10) * ( B+A+C ) );
	BetaD=-h*log(10)*DIC/CO2*B*PhiD;


	Q=(h+2*K2);
	V=(Kb*bor)/((h+Kb)*(h+Kb)) + Kw/(h*h)+1;

	DB=(( K2*(2*CO3+HCO3)+ Q*V *(h+K2)+(h/K1)*( (2*CO3+HCO3)*Q+2*K2*(2*CO3+HCO3)+h*Q*V))/Q)*(1/(Q-(h+K2+h*h/K1)))-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DB)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DB)/((h+2*K2)*(h+2*K2));
	PhiB=-1/(h*log(10) * ( B+A+C ) );
	BetaB=-h*log(10)*DIC/CO2*B*PhiB;

	DC=2*(( K2*(2*CO3+HCO3)+ Q*V *(h+K2)+(h/K1)*( (2*CO3+HCO3)*Q+2*K2*(2*CO3+HCO3)+h*Q*V))/Q)*(1/(Q-2*(h+K2+h*h/K1)))-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DC)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DC)/((h+2*K2)*(h+2*K2));
	PhiC=-1/(h*log(10) * ( B+A+C ) );
	BetaC=-h*log(10)*DIC/CO2*B*PhiC;

	D1=(K1*(K1*K2-h*h)*DIC)   /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	D2=(-K1*K2*(2*h+K1)*DIC)  /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	D=D1+2*D2;
	PhiH=1/ (h*log(10)* (D +(-Kb*bor/((h+Kb)*(h+Kb)))  + (-Kw/(h*h))-1))  ; 

	Pi=(h*K1*(h+2*K2)*DIC)  /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	PiH=((-h/Kh)*log(10)*Pi)*PhiH;
	PiB=CO2/(Kh*DIC)*BetaB;
	PiD=CO2/(Kh*DIC)*BetaD;
	PiC=CO2/(Kh*DIC)*BetaC;
	
	col <- c("PhiD", "BetaD", "PiD", "PhiB", "BetaB", "PiB", "PhiC", "BetaC", "PiC", "PhiH", "PiH")
	res <- data.frame(PhiD,BetaD,PiD,PhiB,BetaB,PiB,PhiC,BetaC,PiC,PhiH,PiH)
	RES <- rbind(RES, res)
	}
return(RES)
}


